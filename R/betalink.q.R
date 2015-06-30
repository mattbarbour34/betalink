### Quantitative version of the betalink function.

# Paul Rabie adapted this code to accomodate quantitative data, and Matt Barbour debugged Paul's code to make it work.

# Treatment of beta_WN, beta_ST and contrib is the same as in betalink.b  See notes associated with that function. 
# Jan 12: commented out "if common species > 0 code", because it doesn't make sense to have it this way.
betalink.q = function(w1,w2,bf="jaccard"){
  w1 <- as.matrix(w1) # need to be converted into matrices for "as.table" function to work for quantitative data.
  w2 <- as.matrix(w2) # need to be converted into matrices for "as.table" function to work for quantitative data.
  
  sp1 = list(top=rowSums(w1),bottom=colSums(w1),all=c(rowSums(w1),colSums(w1)[!names(colSums(w1))%in%names(rowSums(w1))])) # note that rowSums correspond to "top" and colSums correspond to "bottom". This is not that intuitive
  sp2 = list(top=rowSums(w2),bottom=colSums(w2),all=c(rowSums(w2),colSums(w2)[!names(colSums(w2))%in%names(rowSums(w2))]))
  topCom=vec2data.frame(sp1$top,sp2$top,all=T)
  bottomCom=vec2data.frame(sp1$bottom,sp2$bottom,all=T)
  allCom=vec2data.frame(sp1$all,sp2$all,all=T)  
  topCom[is.na(topCom)]=0
  bottomCom[is.na(bottomCom)]=0
  allCom[is.na(allCom)]=0
  beta_U = vegdist(topCom,method=bf,diag=T,upper=T) # based on nodes in the ROWS
  beta_L = vegdist(bottomCom,method=bf,diag=T,upper=T) # based on nodes in the COLUMNS
  beta_S = vegdist(allCom,method=bf,diag=T,upper=T)
  
  #browser()
  # Common species
  Csp = names(sp1$all)[names(sp1$all) %in% names(sp2$all)]
  CUsp = names(sp1$top)[names(sp1$top) %in% names(sp2$top)]
  CLsp =names(sp1$bottom)[names(sp1$bottom) %in% names(sp2$bottom)]

  # I commented out the 'if' statement below, because it doesn't seem necessary to me to have this as a requirement for running this function.
  if((length(CUsp)>0) & (length(CLsp)>0))
  {
    w1ConInts = as.data.frame(as.table(w1[CUsp,CLsp]))$Freq 
    names(w1ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
    w2ConInts= as.data.frame(as.table(w2[CUsp,CLsp]))$Freq 
    names(w2ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
    commonTaxaWeb=vec2data.frame(w1ConInts,w2ConInts,all=T) 
    # I don't understand why this code is necessary. Again, it makes the code inconsistent with the full range of dissimilarity measures.
    #if(all(rowSums(commonTaxaWeb)>0)){ #  formerly there
    beta_OS = vegdist(commonTaxaWeb,method=bf,diag=T,upper=T)#} else   {
    # beta_OS = 1
    # } 
  } else {
    beta_OS = 0
  }
   
  w1=as.data.frame(as.table(w1)) # update from Matt Barbour to Paul Rabie's code
  w1Ints = w1$Freq
  names(w1Ints) = paste(w1$Var1,w1$Var2,sep="-")
  w2=as.data.frame(as.table(w2)) # update from Matt Barbour to Paul Rabie's code
  w2Ints = w2$Freq
  names(w2Ints) = paste(w2$Var1,w2$Var2,sep="-")
  ints=vec2data.frame(w1Ints,w2Ints,all=T)
  ints[is.na(ints)]=0  
  beta_WN = vegdist(ints,method=bf,diag=T,upper=T)
  
  #browser()
    


    beta_ST = beta_WN - beta_OS
    if(beta_WN > 0){
      b_contrib = beta_ST / beta_WN
    } else {
      b_contrib = 0
    }
  # I have commented out everything below, because if someone uses Euclidean distance, then beta_WN and beta_ST can be greater than 1.
  #} else {
  #  beta_WN = 1
  #  beta_OS = 0
   # beta_ST = 1
   # b_contrib = 1 
  #}
  return(list(U = beta_U[1], L = beta_L[1], S = beta_S[1], OS = beta_OS[1], WN = beta_WN[1], ST = beta_ST[1], contrib = b_contrib[1]))
}
