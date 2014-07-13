# now amenable to quantitative metrics

beta.os_prime = function(W, bf="jaccard"){
	metaweb = metaweb(W)$web
	os_prime = NULL
	for(w in W){
		partitions = betalink(as.table(w), as.table(metaweb[rownames(w),colnames(w)]), bf)
		os_prime = c(os_prime,partitions$OS)
	}
	return(os_prime)
}