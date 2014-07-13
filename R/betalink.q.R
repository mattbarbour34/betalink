### Quantitative version of the betalink function.

# Paul Rabie adapted this code to accomodate quantitative data, and Matt Barbour debugged Paul's code to make it work.

# Treatment of beta_WN, beta_ST and contrib is the same as in betalink.b  See notes associated with that function.  
betalink.q = function(w1,w2,bf="jaccard"){
  sp1 = list(top=rowSums(w1),bottom=colSums(w1),all=c(rowSums(w1),colSums(w1)[!names(colSums(w1))%in%names(rowSums(w1))]))
  sp2 = list(top=rowSums(w2),bottom=colSums(w2),all=c(rowSums(w2),colSums(w2)[!names(colSums(w2))%in%names(rowSums(w2))]))
  topCom=vec2data.frame(sp1$top,sp2$top,all=T)
  bottomCom=vec2data.frame(sp1$bottom,sp2$bottom,all=T)
  allCom=vec2data.frame(sp1$all,sp2$all,all=T)  
  topCom[is.na(topCom)]=0
  bottomCom[is.na(bottomCom)]=0
  allCom[is.na(allCom)]=0
  beta_U = vegdist(topCom,method=bf,diag=T,upper=T)
  beta_L = vegdist(bottomCom,method=bf,diag=T,upper=T)
  beta_S = vegdist(allCom,method=bf,diag=T,upper=T)
  # Common species
  Csp = names(sp1$all)[names(sp1$all) %in% names(sp2$all)]
  CUsp = names(sp1$top)[names(sp1$top) %in% names(sp2$top)]
  CLsp =names(sp1$bottom)[names(sp1$bottom) %in% names(sp2$bottom)]
  if((length(CUsp)>0) & (length(CLsp)>0))
  {
    w1ConInts = as.data.frame(as.table(w1[CUsp,CLsp]))$Freq 
    names(w1ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
    w2ConInts= as.data.frame(as.table(w2[CUsp,CLsp]))$Freq 
    names(w2ConInts) = as.vector(outer(CUsp,CLsp,FUN="paste",sep="-"))
    commonTaxaWeb=vec2data.frame(w1ConInts,w2ConInts,all=T)   
    w1=as.data.frame(as.table(w1)) # update from Matt Barbour to Paul Rabie's code
    w1Ints = w1$Freq
    names(w1Ints) = paste(w1$Var1,w1$Var2,sep="-")
    w2=as.data.frame(as.table(w2)) # update from Matt Barbour to Paul Rabie's code
    w2Ints = w2$Freq
    names(w2Ints) = paste(w2$Var1,w2$Var2,sep="-")
    ints=vec2data.frame(w1Ints,w2Ints,all=T)
    ints[is.na(ints)]=0  
    beta_WN = vegdist(ints,method=bf,diag=T,upper=T)
    if(all(rowSums(commonTaxaWeb)>0)){
      beta_OS = vegdist(commonTaxaWeb,method=bf,diag=T,upper=T)} else   {
        beta_OS = 1
      }
    beta_ST = beta_WN - beta_OS
    if(beta_WN > 0){
      b_contrib = beta_ST / beta_WN
    } else {
      b_contrib = 0
    }
  } else {
    beta_WN = 1
    beta_OS = 0
    beta_ST = 1
    b_contrib = 1 
  }
  return(list(U = beta_U[1], L = beta_L[1], S = beta_S[1], OS = beta_OS[1], WN = beta_WN[1], ST = beta_ST[1], contrib = b_contrib[1]))
}
