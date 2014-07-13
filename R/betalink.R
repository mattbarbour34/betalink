# betalink is now a dispatch function that calls its binary version if the bf argument is itself a function, and calls its quantitative version if the bf argument is a character string

betalink = function(w1, w2, bf="jaccard"){
  if(is.character(bf)) betalink.q(w1,w2,bf=bf) else betalink.b(w1,w2,bf=bf)
}
