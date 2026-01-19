
# compuatation of the autocorrelation
autocor=function(x){abs(cor(x[-1],x[-length(x)]))} 


### computation of the lindley process from scores

lindley=function(scores){
  L=length(scores)
	sl=rep(0,L+1)
	for (i in 1:L){
		sl[i+1]= max(0, (sl[i]+scores[i]))
		}
	return(sl[-1])
	}


# computation of the local score from a lindley process and the region under 
# the maximum peak

scorelocal=function(lind){
	nC=max(lind[,1])
	sloc=matrix(0, nrow=nC, ncol=3)
	colnames(sloc)=c('chr','beg','end')
	for (i in 1:nC){
		tmp=which(lind[,1]==i)		
		list=lind[tmp,3]
		M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))>0){
			m_loc=max(which(list[1:M_loc]==0))+1			
			}else{m_loc=1}
	sloc[i,]=c(i,lind[tmp[m_loc],2],lind[tmp[M_loc],2])	
	}
	return(data.frame(zone=seq(1,nrow(sloc)),sloc))
}

#computation of the significance threshold if the distribution of the p-values is uniform

thresUnif=function(L, cor, xi, alpha = 0.05){
  if (xi==1){ 
    a=log(L) - 5.5*cor^3 + 6.76*cor^2 - 5.66*cor -2.51
    b = -1.22*cor^2 + 3.17*cor - 1.99
      }else{if (xi==2){ 
      a   = log(L) + 2.47*cor^3 - 4.16*cor^2 - 1.82*cor - 4.58  
      b = 0.37*cor^2 + 2.14*cor - 2.35}else{print('xi not equal to 1 or 2')}
      }       
  #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
}

# computation of the significative regions from a lindley process given a significance threshold

sig_sl=function(lind,pos, th){
	zones=c(0,0,0)	
	list=lind
	auxpos=pos
	while(max(list)>=th){
	  M_loc=which.max(list)
		if(length(which(list[1:M_loc]==0))==0){ #the peak is at the beginning of the chrom
			m_loc=1
			zones=rbind(zones, c(auxpos[m_loc],auxpos[M_loc],max(list)))
			tmp=which.min[which(list[M_loc+1:length(list)]==0)] #first 0 score after peak
			list=list[tmp:length(list)]
			auxpos=pos[tmp:length(list)]
			}else{	
				m_loc=max(which(list[1:M_loc]==0))			
				max=max(list)
				zones=rbind(zones, c(auxpos[m_loc+1],auxpos[M_loc],max))
				tmp=which(list[M_loc:length(list)]==0) #first 0 score after peak
				if (length(tmp)>0){
				  auxpos=auxpos[c(1:m_loc,(min(tmp)+M_loc):length(list))]
				  list=list[c(1:m_loc, (min(tmp)+M_loc):length(list))]
				  }else{ #if the peak is at the end of the chromosome
				    auxpos=auxpos[1:m_loc]
				    list=list[1:m_loc]
				    }				
				}
	  }
	zones=matrix(zones, ncol=3)
	zones=data.table(beg=zones[,1],end=zones[,2],peak=zones[,3])
	if (nrow(zones)>1){zones=zones[-1,]}
	return(zones)
	}

### Estimation of Gumble coeficients.

# Estimation of the coeficients of the Gumbel distributions

lineGumb=function(x){
  x1tmp=x
  if (length(table(x1tmp)) >  5){
    Fn = ecdf(x1tmp)
    Fnx=Fn(x1tmp)
    filtre = ( Fnx < 1 - 10**(-1) )  & ( x1tmp > 1 ) #sinon pbs numq
    lm0 = lm(log(-log(Fnx[filtre])) ~ x1tmp[filtre])$coefficients
  }else{lm0=c(0,0)}  
  return(lm0)
}

# Estimation of the coeficients of the polynomes to compute the Gumbel coefficients 
# depending on the length of the chromosomes and the chromosome autocorrelation

coefsGumb=function(mydata, Ls=seq(10000,80000,10000), nSeq=5000){
  bins=seq(0.025,9.975,0.05) 
  coefs=array(0, dim=c(length(bins),5, length(Ls)))  
  coefs[,5,]=matrix(rep(Ls,length(bins)), ncol=length(Ls),byrow=TRUE)
  as=matrix(0, ncol=length(Ls), nrow=length(bins))
  bs=matrix(0, ncol=length(Ls), nrow=length(bins))
  xx=seq(0,max(mydata$lindley), 1)
  for (j in 1:length(Ls)){
    coefs[,1,j]=bins  
    len=Ls[j]
    tmp=sample(nrow(mydata)-len,nSeq,replace=F) #samples the sequences'beginnings
    DTL=data.table(LocScore = vector('numeric', length=nSeq),cors= vector('numeric', length=nSeq)) 
    for (l in 1:nSeq){
      DTL[l]=mydata[seq(tmp[l],length.out=len),.(max(lindley(score)),autocor(pval))]
    }
    DTL[,bin:=which.min(abs(bins-cors)),cors]
    binNE=sort((unique(DTL$bin)))
    coefs[binNE,4,j]=table(DTL$bin)
    DTCoef=DTL[,.(lineGumb(LocScore)),bin]
    coefs[unique(DTL$bin),c(2,3),j]=matrix(DTCoef$V1, ncol=2, byrow=T)
    ys=coefs[,c(2,3),j]%*%rbind(rep(1,length(xx)), xx)
    pdf(paste0('GumbelinesAllxi',xi,'L',len,'.pdf'))
    for (i in binNE){
      Fn=ecdf(DTL[bin==i,LocScore]) 
      tmp=DTL[bin==i,.(LocScore,Fn=Fn(LocScore)),][(Fn< 1 - 10**(-5)) & (LocScore > 1)]
      if(i==min(binNE)){
        plot(tmp$LocScore, log(-log(tmp$Fn)), col=i,main=paste("L = ",len,sep=""), xlim=quantile(DTL$LocScore,c(0,0.9)),ylim=range(ys[binNE,xx<quantile(DTL$LocScore,0.8)], na.rm=T), pch=20, xlab='Local Score', ylab='log(-log(Fn(LS)))')
        }else{points(tmp$LocScore, log(-log(tmp$Fn)), col=i, pch=20)
      }
      if (sum(ys[i,]!=0)){lines(xx, ys[i,], col=i, lwd=1.5)}
    }
    dev.off()
    as[,j]=as.numeric(coefs[,2,j])
    bs[,j]=as.numeric(coefs[,3,j])		
  }
  
  auxWhich=which(coefs[,4,]>100)
  rhos=coefs[,1,][auxWhich]
  auxAs=as[auxWhich]-log(coefs[,5,][auxWhich])
  auxBs=bs[auxWhich]
  
  fitA=lm((auxAs)~I(rhos)+I(rhos^2))  
  fitB=lm((auxBs)~I(rhos)+I(rhos^2))  
  
  pdf(paste0('FitAandB.pdf'))
  par(mfrow=c(1,2))
  plot(rhos, auxAs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='a-log(L)')
  xslm=seq(0,1,0.01)
  lines(xslm, fitA$coefficients%*%rbind(rep(1,length(xslm)), xslm, xslm^2))
  plot(rhos, auxBs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='b')
  lines(xslm, fitB$coefficients%*%rbind(rep(1,length(xslm)), xslm, xslm^2))
  legend(min(rhos), max(auxBs), legend=sort(unique(coefs[,5,][auxWhich])), col=unique(coefs[,5,][auxWhich]/10000), pch=16)
  
  dev.off()
  
  return(list(aCoef=fitA$coefficients, bCoef=fitB$coefficients))
  }

# Computation of the significance threshold given the computed polynomes for computing
# the Gumbel coefficients depending of length and autocorrelation

threshold=function(L, cor, aCoef, bCoef, alpha = 0.05){
  degA=length(aCoef)
  degB=length(bCoef)
  a=log(L)
  b=0
  for (i in 1:degA){
    a=a+aCoef[i]*cor^(i-1)  
    }
  for (i in 1:degB){
    b=b+bCoef[i]*cor^(i-1)  
  }
  #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
  }


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Function to create a Manhattan plot with chromosome-based coloring
create_manhattan_plot <- function(data, 
                                  title = "Manhattan Plot",
                                  threshold = -log10(0.05/nrow(data)),
                                  chr_column = "CHR",
                                  bp_column = "BP",
                                  p_column = "P",
                                  snp_column = NULL) {
  
  # Check if the required columns exist in the data
  required_cols <- c(chr_column, bp_column, p_column)
  if (!all(required_cols %in% colnames(data))) {
    stop("Required columns not found in data. Please ensure your data contains chromosome, position, and p-value columns.")
  }
  
  # Convert data for plotting
  data_plot <- data %>%
    # Ensure chromosome is a factor for proper ordering
    mutate(CHR = factor(.data[[chr_column]], levels = sort(unique(as.factor(.data[[chr_column]])))))
  
  # Calculate cumulative position for x-axis
  data_plot <- data_plot %>%
    group_by(CHR) %>%
    arrange(CHR, .data[[bp_column]]) %>%
    ungroup()
  
  # Compute cumulative positions for plotting
  chr_lengths <- data_plot %>%
    group_by(CHR) %>%
    summarize(max_pos = max(.data[[bp_column]]), .groups = "drop")
  
  chr_cumsum <- cumsum(c(0, chr_lengths$max_pos[-nrow(chr_lengths)]))
  names(chr_cumsum) <- chr_lengths$CHR
  
  data_plot$pos_cum <- as.numeric(data_plot[[bp_column]]) + chr_cumsum[as.character(data_plot$CHR)]
  
  # Calculate midpoints for chromosome labels
  axis_df <- data_plot %>%
    group_by(CHR) %>%
    summarize(center = mean(pos_cum), .groups = "drop")
  
  # Create the plot
  p <- ggplot(data_plot, aes(x = pos_cum, y = -log10(.data[[p_column]]), color = CHR)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = rep(c("#0066CC", "#009900", "#CC0000", "#FF9900", "#663300", "#000000", "#9900CC", "#00CCCC"), 
                                    length.out = length(unique(data_plot$CHR)))) +
    scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "gray") +
    labs(title = title,
         x = "Chromosome",
         y = expression(-log[10](p))) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5)
    )
  
  # Add labels for significant points if SNP column is provided and there are significant points
  if (!is.null(snp_column) && any(-log10(data_plot[[p_column]]) > threshold)) {
    signif_points <- data_plot %>% 
      filter(-log10(.data[[p_column]]) > threshold) %>%
      arrange(desc(-log10(.data[[p_column]]))) %>%
      head(10)  # Only label top 10 most significant points
    
    p <- p + geom_text(data = signif_points, 
                       aes(label = .data[[snp_column]]), 
                       vjust = -1, size = 3)
  }
  
  return(p)
}
