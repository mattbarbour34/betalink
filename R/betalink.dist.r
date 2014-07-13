# Modified so that betalink.dist returns distance matrices for upper and lower trophic levels in the bipartite matrix.  
# added an option to accommodate quantitative metrics
# Also added an option for rectangular output of distance matrices

betalink.dist = function(W,triangular=F,bf="jaccard"){
	dWN = matrix(NA,ncol=length(W),nrow=length(W))
	colnames(dWN) = names(W)
	rownames(dWN) = names(W)
	diag(dWN)=0
	dOS = dWN
	dS_all.taxa = dWN
	dS_upper = dWN
	dS_lower = dWN
	dST = dWN
	dContrib = dWN
	for(i in 1:(length(W)-1)){
		for(j in (i+1):length(W)){
		  index=cbind(c(j,i),c(i,j))
			partition = betalink(W[[i]],W[[j]],bf)
			dWN[index]		= partition$WN
			dOS[index]		= partition$OS
			dS_all.taxa[index]  = partition$S
			dS_upper[index]   = partition$U
			dS_lower[index]    = partition$L
			dST[index]		= partition$ST
			dContrib[index]	= partition$contrib
		}
	}
	distances = list(WN=dWN, OS=dOS, S_all.taxa=dS_all.taxa, S_upper=dS_upper ,S_lower=dS_lower, ST=dST, contrib=dContrib)
	if(triangular==T) distances = lapply(distances,as.dist)
	return(distances)
}