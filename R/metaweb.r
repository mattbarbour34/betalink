# Paul Rabie modified to make this function amenable to quantitative distance metrics.  
# metaweb is now the sum of observations over all realizations, 
# co.oc is still the number of realizations for which each species pair co-occurs, and 
# null.template is now the expected number of observed interactions per co-occurence.

metaweb = function(W, binary = F){
	Lo = unique(unlist(lapply(W,colnames)))
	Up = unique(unlist(lapply(W,rownames)))
	meta = matrix(0,ncol=length(Lo),nrow=length(Up))
	colnames(meta) = Lo
	rownames(meta) = Up
	co.oc = meta
	for(w in W){
	  if(binary==T)  w[w>0] = 1   ##binarization is now optional.
		meta[rownames(w),colnames(w)] = meta[rownames(w),colnames(w)] + w
		co.oc[rownames(w),colnames(w)] = co.oc[rownames(w),colnames(w)] + 1
	}
	null.template = meta/co.oc
	null.template[is.nan(null.template)] = 0
	# meta[meta>0] = 1 # removed binarization
	return(list(web=meta,template=null.template, cooc = co.oc))
}