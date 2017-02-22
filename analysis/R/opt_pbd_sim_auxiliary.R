opt1 = function(contr){ # speciation-initiation of good species
  parent = as.numeric(ifelse(length(contr$sg)>1,
                             sample(contr$sg,1),
                             contr$sg))
  contr$id = contr$id + 1
  contr$L = rbind(contr$L, c(contr$id, parent, contr$t, -1, -1,
                             contr$L[abs(parent) - contr$id1, 6]))
  # cat("\n\n")
  # cat(c(contr$id, parent, contr$t, -1, -1, contr$L[abs(parent) - contr$id1, 6]))
  # cat("\n")
  # cat(abs(parent), contr$id1, contr$L[abs(parent) - contr$id1, 6])
  # cat("\n")
  contr$si = c(contr$si, -contr$id)
  contr$Ni = contr$Ni + 1
  
  contr$parent = parent
  
  return(contr)
}
opt2 = function(contr){ # extinction of good species
  todie = as.numeric(ifelse(length(contr$sg)>1,
                            sample(contr$sg,1),
                            contr$sg))
  contr$L[todie - contr$id1, 5] = contr$t
  contr$sg = contr$sg[-which(contr$sg == todie)]
  contr$Ng = contr$Ng - 1
  
  contr$todie = todie
  
  return(contr)
}
opt3 = function(contr){ # speciation-completion
  tocomplete = abs(as.numeric(ifelse(length(contr$si)>1,
                                     sample(contr$si,1),
                                     contr$si)))
  contr$L[tocomplete - contr$id1, 4] = contr$t
  contr$Sid = contr$Sid + 1
  contr$L[tocomplete - contr$id1, 6] = contr$Sid
  contr$sg = c(contr$sg, tocomplete)
  contr$si = contr$si[-which(abs(contr$si) == tocomplete)]
  contr$Ng = contr$Ng + 1
  contr$Ni = contr$Ni - 1
  
  contr$tocomplete = tocomplete
  
  return(contr)
}
opt4 = function(contr){ # speciation-initiation of incipient species
  parent = as.numeric(ifelse(length(contr$si)>1,
                             sample(contr$si,1),
                             contr$si))
  contr$id = contr$id + 1
  contr$L = rbind(contr$L, c(contr$id, parent, contr$t, -1, -1,
                             contr$L[abs(parent) - contr$id1, 6]))
  cat(c(contr$id, parent, contr$t, -1, -1,contr$L[abs(parent) - contr$id1, 6]))
  contr$si = c(contr$si, -contr$id)
  contr$Ni = contr$Ni + 1
  
  contr$parent = parent
  
  return(contr)
}
opt5 = function(contr){ # extinction of incipient species 
  todie = abs(as.numeric(ifelse(length(contr$si)>1,
                                sample(contr$si,1), 
                                contr$si)))
  contr$L[todie - contr$id1, 5] = contr$t
  contr$si = contr$si[-which(abs(contr$si) == todie)]
  contr$Ni = contr$Ni - 1
  
  contr$todie = todie
  
  return(contr)
}




