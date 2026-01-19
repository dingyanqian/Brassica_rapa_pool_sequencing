#!/usr/bin/env python3
"""Flag SNPs exceeding temporal FST thresholds per pool."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freq", default="inputs/raw/allele_frequencies_and_coverage.txt",
                        help="Allele frequency table with AF_pool_G1 and AF_pool_<POOL> columns.")
    parser.add_argument("--ne", default="inputs/raw/Ne.csv",
                        help="CSV listing per-pool Ne across generations (wide format).")
    parser.add_argument("--thresholds", default="outputs/final_analysis/final_thresholds.csv",
                        help="Final thresholds CSV with neutral_fst_q95/neutral_fst_q99 columns.")
    parser.add_argument("--out", default="outputs/final_analysis/fst_threshold_flags.csv",
                        help="Output TSV with per-pool FST exceedance flags.")
    return parser.parse_args()


def read_ne(path: str) -> dict[str, float]:
    ne_map: dict[str, float] = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pool = row.get("TREAT")
            if not pool:
                continue
            gens = [float(val) for key, val in row.items() if key.upper().startswith("G") and val not in ("", None)]
            if not gens:
                raise ValueError(f"No generation columns found for pool {pool}")
            ne_map[pool] = sum(gens) / len(gens)
    return ne_map


def read_thresholds(path: str) -> dict[str, dict[str, float]]:
    thresholds: dict[str, dict[str, float]] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pool = row["pool"]
            thresholds[pool] = {
                "q95": float(row["neutral_fst_q95"]),
                "q99": float(row["neutral_fst_q99"]),
            }
    return thresholds


def temporal_fst(p0: float, pt: float, ne_mean: float) -> float:
    p_bar = 0.5 * (p0 + pt)
    denom = p_bar * (1 - p_bar)
    if denom <= 0:
        return 0.0
    sampling = p0 * (1 - p0) / (2 * ne_mean)
    fst = ((pt - p0) ** 2 - sampling) / denom
    return max(fst, 0.0)


def main():
    args = parse_args()
    ne_map = read_ne(args.ne)
    thresholds = read_thresholds(args.thresholds)

    pools = sorted(set(ne_map.keys()) & set(thresholds.keys()))
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.freq, newline="") as handle, out_path.open("w", newline="") as out_handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "AF_pool_G1" not in fieldnames:
            raise RuntimeError("AF_pool_G1 column missing in frequency table")
        pool_cols = {pool: f"AF_pool_{pool}" for pool in pools}
        missing = [col for col in pool_cols.values() if col not in fieldnames]
        if missing:
            raise RuntimeError(f"Missing pool columns: {missing}")

        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["CHROM", "POS", "pool", "fst", "threshold_q95", "threshold_q99",
                         "pass_q95", "pass_q99"])

        for row in reader:
            chrom = row.get("CHROM")
            pos = row.get("POS")
            p0_raw = row.get("AF_pool_G1")
            if p0_raw in (None, "", "NA"):
                continue
            p0 = float(p0_raw)
            for pool in pools:
                pt_raw = row.get(pool_cols[pool])
                if pt_raw in (None, "", "NA"):
                    continue
                pt = float(pt_raw)
                ne_mean = ne_map[pool]
                fst_val = temporal_fst(p0, pt, ne_mean)
                thr_q95 = thresholds[pool]["q95"]
                thr_q99 = thresholds[pool]["q99"]
                writer.writerow([
                    chrom,
                    pos,
                    pool,
                    f"{fst_val:.6f}",
                    f"{thr_q95:.6f}",
                    f"{thr_q99:.6f}",
                    int(fst_val > thr_q95),
                    int(fst_val > thr_q99),
                ])

    print("FST flags written to", out_path)


if __name__ == "__main__":
    main()
