# if one data frame is a subset of a second data frame, the merge function will return the larger data frame.  This function forces the merge function to preserve all rows from both data sets

vec2data.frame=function(x,y,all=T){
  temp=merge(t(x),t(y),all=all)
  if(nrow(temp)<sum(min(nrow(x),1),min(nrow(y),1))){
    x=c(a = -1,x)
    y=c(a = -2,y)
    temp=merge(t(x),t(y),all=all)
    n=names(temp)
    temp=data.frame(temp[,-1])
    names(temp)=n[-1]}
  return(temp)}