#!/usr/bin/env python3
"""Wright–Fisher neutral simulations for treatment-specific Ne trajectories.

This script samples initial allele frequencies from the empirical G1 spectrum
and propagates them for six generations using the provided Ne values.
It produces null distributions for |ΔAF| and FST-like divergence per treatment
and writes summary tables and optional full simulation traces.

Example:
    python3 scripts/simulate_neutral.py \
        --ne inputs/raw/Ne.csv \
        --freq inputs/raw/allele_frequencies_and_coverage.txt \
        --out_dir outputs/neutral_sims \
        --replicates 10000 \
        --seed 20240816
"""

import argparse
import csv
import math
import os
import random
from collections import defaultdict, namedtuple
from pathlib import Path

import numpy as np

SimulationConfig = namedtuple(
    "SimulationConfig",
    ["ne_path", "freq_path", "out_dir", "replicates", "seed", "full_trace", "chunk_size"],
)


def parse_args() -> SimulationConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ne", default="inputs/raw/Ne.csv",
                        help="Path to treatment × generation Ne table (wide format).")
    parser.add_argument("--freq", default="inputs/raw/allele_frequencies_and_coverage.txt",
                        help="Path to allele frequency table with AF_pool_G1 column.")
    parser.add_argument("--out_dir", default="outputs/neutral_sims",
                        help="Directory to write simulation outputs.")
    parser.add_argument("--replicates", type=int, default=10000,
                        help="Number of SNP replicates to simulate per treatment.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    parser.add_argument("--full_trace", action="store_true",
                        help="Write per-generation trajectories for each replicate (large files).")
    parser.add_argument("--chunk_size", type=int, default=1000,
                        help="Batch size when writing full traces or sampling p0.")
    args = parser.parse_args()
    return SimulationConfig(
        ne_path=args.ne,
        freq_path=args.freq,
        out_dir=args.out_dir,
        replicates=args.replicates,
        seed=args.seed,
        full_trace=args.full_trace,
        chunk_size=args.chunk_size,
    )


def read_ne_table(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Ne table not found at {path}")
    treatments = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            treat = row.get("TREAT")
            if not treat:
                raise ValueError("Ne table requires a 'TREAT' column")
            gens = []
            for key in sorted(row.keys()):
                if key.upper().startswith("G"):
                    val = row[key]
                    if val in ("", None):
                        raise ValueError(f"Missing Ne for {treat} {key}")
                    ne_val = float(val)
                    if ne_val <= 0:
                        raise ValueError(f"Non-positive Ne for {treat} {key}: {ne_val}")
                    gens.append(ne_val)
            if not gens:
                raise ValueError(f"No generation columns found for treatment {treat}")
            treatments[treat] = gens
    return treatments


def sample_p0(freq_path: str, chunk_size: int):
    if not os.path.exists(freq_path):
        raise FileNotFoundError(f"Allele frequency table not found at {freq_path}")

    p0_values = []
    with open(freq_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "AF_pool_G1" not in reader.fieldnames:
            raise ValueError("Allele frequency table must contain 'AF_pool_G1'")
        for row in reader:
            val = row.get("AF_pool_G1")
            if val in (None, "", "NA"):
                continue
            p0 = float(val)
            if 0.0 <= p0 <= 1.0:
                p0_values.append(p0)
            if len(p0_values) >= chunk_size:
                yield from p0_values
                p0_values.clear()
    if p0_values:
        yield from p0_values


def drift_step(p, ne):
    """One Wright–Fisher generation with effective size ne."""
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    count = np.random.binomial(2 * int(round(ne)), p)
    return count / (2 * int(round(ne)))


def simulate_replicate(p0, ne_values):
    traj = [p0]
    p = p0
    for ne in ne_values:
        p = drift_step(p, ne)
        traj.append(p)
    return traj


def weir_cockerham_fst(p0, pt, ne_mean):
    """Approximate FST using temporal variance."""
    # Temporal FST under neutral drift: Var(p_t - p_0) / (p_bar(1-p_bar))
    p_bar = 0.5 * (p0 + pt)
    denom = p_bar * (1 - p_bar)
    if denom <= 0:
        return 0.0
    return ((pt - p0) ** 2 - (p0 * (1 - p0) / (2 * ne_mean))) / denom


def main():
    cfg = parse_args()

    np.random.seed(cfg.seed)
    random.seed(cfg.seed)

    ne_map = read_ne_table(cfg.ne_path)
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    trace_dir = out_dir / "traces"

    p0_pool = list(sample_p0(cfg.freq_path, cfg.chunk_size))
    if not p0_pool:
        raise RuntimeError("No baseline allele frequencies found")

    for treat, ne_values in ne_map.items():
        ne_array = np.asarray(ne_values, dtype=float)
        ne_mean = float(np.mean(ne_array))
        abs_deltas = []
        fst_values = []

        if cfg.full_trace:
            trace_dir.mkdir(parents=True, exist_ok=True)
            trace_file = trace_dir / f"{treat}_traces.csv"
            trace_handle = trace_file.open("w")
            trace_writer = csv.writer(trace_handle)
            header = ["replicate", "generation", "allele_frequency"]
            trace_writer.writerow(header)
        else:
            trace_handle = None
            trace_writer = None

        for rep in range(cfg.replicates):
            p0 = random.choice(p0_pool)
            traj = simulate_replicate(p0, ne_array)
            pt = traj[-1]
            abs_deltas.append(abs(pt - p0))
            fst_values.append(weir_cockerham_fst(p0, pt, ne_mean))

            if trace_writer:
                for gen_idx, freq in enumerate(traj):
                    trace_writer.writerow([rep, gen_idx, f"{freq:.6f}"])

        if trace_handle:
            trace_handle.close()

        abs_array = np.asarray(abs_deltas)
        fst_array = np.asarray(fst_values)

        summary_rows.append({
            "pool": treat,
            "replicates": cfg.replicates,
            "mean_abs_delta": float(abs_array.mean()),
            "p95_abs_delta": float(np.quantile(abs_array, 0.95)),
            "p99_abs_delta": float(np.quantile(abs_array, 0.99)),
            "mean_fst": float(fst_array.mean()),
            "p95_fst": float(np.quantile(fst_array, 0.95)),
            "p99_fst": float(np.quantile(fst_array, 0.99)),
        })

        np.save(out_dir / f"{treat}_abs_delta.npy", abs_array)
        np.save(out_dir / f"{treat}_fst.npy", fst_array)

    summary_path = out_dir / "neutral_sim_summary.csv"
    with summary_path.open("w", newline="") as f:
        fieldnames = [
            "pool", "replicates", "mean_abs_delta", "p95_abs_delta", "p99_abs_delta",
            "mean_fst", "p95_fst", "p99_fst"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Simulation complete. Results saved to {out_dir}")


if __name__ == "__main__":
    main()
