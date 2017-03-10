

# load("analysis/R/2017_03_01_expand_rcpp/rcpp/pbd_sim_cpp0.3_0.3_0_0_0.1.RData")
# L = out[[9]]$L
# age = out[[9]]$L[1, 3]
# my = get_phylo(L, age)
# plot(my)
# axisPhylo()

# set.seed(69)
# age = 10
# aux = pbd_sim(c(0.7, 0.3, 0, 0, 0.5), age = age, soc = 1)
# (L = aux$L)
# 
# plot(aux$stree_oldest)
# aux$stree_oldest$edge
# aux$stree_oldest$edge.length
# 
# aaa = get_phylo(L, age)
# plot(aaa, show.node.label = 1)
# axisPhylo()
# 
# all.equal(aaa$edge.length, aux$stree_oldest$edge.length)
# cophyloplot(aaa, aux$stree_oldest)

# test = list(L, age)

library(PBD)

get_phylo = function(L, age){
  sp = unique(L[ , 6])
  Nsp = length(sp)
  good = get_good(L, sp)
  good = cbind(good, NA)
  good[2, 7] = Nsp + 1
  daugthers = lapply(good[ , 1], function(x) which(good[ , 2] == x))
  pp = list("edge" = matrix(ncol = 2), "Nnode" = Nsp+1, "mat" = good)
  ppp = build_edge(pp, daugthers, Nsp+1, 1)
  edgll = build_edge.length(ppp$edge, ppp$mat, Nsp)
  phy = list("edge" = ppp$edge, "Nnode" = ppp$Nnode, "edge.length" = edgll)
  phy = clean_phylo(phy, Nsp, good[ , 6])
  
  return(phy)
}

clean_phylo = function(ttree, ns, taxa){
  # get the length of the root
  ttree$root.edge = ttree$edge.length[1]
  
  # removes root from 'edge.length' and 'edge'
  ttree$edge.length = ttree$edge.length[-1]
  ttree$edge = ttree$edge[-c(1:2), ]
  
  # re-calculates the node numbers
  ttree$edge[ , 1] = ttree$edge[ , 1] - 1
  ttree$edge[ , 2] = sapply(ttree$edge[ , 2], function(x) ifelse(x <= ns, x, x - 1))
  
  # adds 'Nnode' to the phylogeny
  ttree$Nnode = ns - 1
  
  # creates 'tip.label'
  ttree$tip.label = paste0("sp.", taxa)
  
  # assigns the class of the object
  class(ttree) = "phylo"
  
  return(ttree)
}

#edge=ppp$edge; mat=ppp$mat; tips=Nsp
build_edge.length = function(edge, mat, tips){
  out = c()
  edge = edge[-1, ]
  for(w in 1:nrow(edge)){
    leaf = edge[w, 2]
    stem = edge[w, 1]
    if(leaf <= tips){
      out[w] = mat[which(mat[ , 7] == stem), 3] - mat[which(mat[ , 6] == leaf), 5]
    } else{
      out[w] = mat[which(mat[ , 7] == leaf), 4]
    }
  }
  
  return(out)
}

# filog=pp;family =daugthers;nn= Nsp+1;id =1
build_edge = function(filog, family, nn, id){
  # get the list of cladogenetic events
  fam = family[[id]]
  
  bb = length(fam) > 0
  if(bb){ # if there's at least one cladogenesis along this branch
    # creates the list of nodes
    nodes = filog$Nnode+(1:length(fam))
    nodes.old = c(nn, nodes[-length(nodes)])
    # for each cladogenetic event; ie, node
    # adds event/node to the 'edge'
    filog$edge = rbind(filog$edge, matrix(c(nodes.old, nodes), ncol = 2))
    
    # gets info about the node names
    filog$mat[fam, 7] = nodes
    
    # updates the Nnode
    filog$Nnode = rev(nodes)[1]
    nn = filog$Nnode
  } 
  # add the tip to the phylogeny
  filog$edge = rbind(filog$edge, c(nn, id))
  
  if(bb){ # now it opens the left-hand side of the phylogeny
    for(k in rev(fam)){ # goes backwards checking for more cladogenetic events
      filog = build_edge(filog, family, nn, k)
      nn = nn-1
    }
  }
  
  return(filog)
}

#mat=L;taxa=sp
get_good = function(mat, taxa){
  incip = lapply(taxa, FUN = function(s) mat[which(mat[ , 6] == s), ])
  orig = sapply(incip, function(x) ifelse(class(x) == "numeric", x[3], max(x[ , 3])))
  ext = sapply(incip, function(x) ifelse(class(x) == "numeric", x[5], x[ , 5]))
  ext[ext < 0] = 0
  out = lapply(incip, function(x) if(class(x) == "numeric"){x}else{x[1, ]})
  out = do.call(rbind, out)
  out[ , 3] = orig
  out[ , 5] = ext
  out[ , 4] = get_length(out)
  
  return(out)
}

#mat=out
get_length = function(mat){
  out = 0
  parents = unique(mat[ , 2])
  id.par = match(parents, mat[ , 1])
  all.id = c()
  for(i in 2:length(parents)){
    id = which(mat[ , 2] == parents[i])
    all.id = c(all.id, id)
    if(length(id) == 1){
      out = c(out, mat[id.par[i], 3] - mat[id, 3])
    } else{
      aux = c(mat[id.par[i], 3], mat[id, 3])
      aux = abs(diff(aux))
      out = c(out, aux)
    }
  }
  
  return(out[order(c(1, all.id))])
}