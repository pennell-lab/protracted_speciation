

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

require(PBD)

get2phylo = function(L, age){
  ll = length(unique(L[ , 6]))
  if(ll >= 2){
    out = get_phylo(L, age)
  } else{
    edge = matrix(ncol = 2)
    edge.length = NA
    Nnode = 0
    root.edge = max(L[ , 3]) - max(min(L[ , 5]), 0)
    tip.label = "sp.1"
    # } else{
    #   edge = matrix(c(3, 3, 1, 2), ncol = 2)
    #   edge.length = rep(max(L[which(L[ , 6] == 2), 3]), 2)
    #   edge.length = edge.length - sapply(1:2, function(x) max(min(L[which(L[ , 6] == x), 5]), 0))
    #   Nnode = 1
    #   root.edge = L[1, 1] - max(L[which(L[ , 6] == 2), 3])
    #   tip.label = c("sp.1", "sp.2")
    # }
    
    out = list("edge" = edge, "edge.length" = edge.length, "Nnode" = Nnode, "tip.label" = tip.label, "root.edge" = root.edge)
  }
  
  return(out)
}

get_phylo = function(L, age){
  sp = unique(L[ , 6])
  Nsp = length(sp)
  good = get_good(L, sp)
  good = cbind(good, NA)
  good[2, 7] = Nsp + 1
  daugthers = lapply(good[ , 1], function(x) which(good[ , 2] == x))
  names(daugthers) = sp
  pp = list("edge" = matrix(ncol = 2), "Nnode" = Nsp+1, "mat" = good)
  ppp = build_edge(pp, daugthers, Nsp+1, 1, names(daugthers))
  edgll = build_edge.length(ppp$edge, ppp$mat, Nsp)
  phy = list("edge" = ppp$edge, "Nnode" = ppp$Nnode, "edge.length" = edgll)
  phy = clean_phylo(phy, Nsp, 1:Nsp)#na.exclude(phy$edge[phy$edge[ , 2] <= Nsp, 2]))
  
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
  #mat[ , 5] = mat[ , 5] - mat[1, 3]
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
build_edge = function(filog, family, nn, id, taxa){
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
  filog$edge = rbind(filog$edge, c(nn, as.numeric(taxa[id])))
  
  if(bb){ # now it opens the left-hand side of the phylogeny
    for(k in rev(fam)){ # goes backwards checking for more cladogenetic events
      filog = build_edge(filog, family, nn, k, taxa)
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





############ Funcions for 2 separated analysis
# For "opt_pbd_sim_taxa2tts.cpp" --> pbdLoop_taxa2tts

#phy1=aa[[1]];phy2=aa[[2]]
combineL = function(phy1, phy2){
  if(phy1$Nnode == 0){
    Nnode = phy2$Nnode + 1
    tips = c(1, Ntip(phy2))
    edge = phy2$edge
    edge[ , 1] = edge[ , 1] + 2
    edge[ , 2] = sapply(edge[ , 2], function(x) ifelse(x > tips[2], x + 2, x))
    edge = rbind(c(sum(tips) + 1, tips[2] + 1),
                 c(sum(tips) + 1, sum(tips) + 2),
                 edge) # add node that conects the 2 phylos
    # creates the edge.length
    edge.length = c(phy1$root.edge, phy2$root.edge, phy2$edge.length)
  } else if(phy2$Nnode == 0){
    Nnode = phy1$Nnode + 1
    tips = c(Ntip(phy1), 1)
    edge = phy1$edge
    edge[ , 1] = edge[ , 1] + 2
    edge[ , 2] = sapply(edge[ , 2], function(x) ifelse(x > tips[1], x + 2, x))
    edge = rbind(c(sum(tips) + 1, tips[1] + 1),
                 c(sum(tips) + 1, sum(tips) + 2),
                 edge) # add node that conects the 2 phylos
    # creates the edge.length
    edge.length = c(phy2$root.edge, phy1$root.edge, phy1$edge.length)
  } else{
    Nnode = phy1$Nnode + phy2$Nnode + 1
    tips = sapply(list(phy1, phy2), function(x) sum(!x$edge[ , 2] %in% x$edge[ , 1]))
    # use the 1st phylo as reference
    edge = phy1$edge
    
    # modify "edge" to include the nodes and tips from the other phylo
    edge[ , 1] = edge[ , 1] + tips[2] + 1
    edge[ , 2] = sapply(edge[ , 2], function(x) ifelse(x > tips[1], x + tips[2] + 1, x))
    edge = rbind(c(sum(tips) + 1, sum(tips) + 2), edge) # add node that conects the 2 phylos
    # modify the edge from phylo2 to include the nodes and tips from the other phylo
    ee = phy2$edge
    ee[ , 1] = ee[ , 1] + max(edge) - tips[2]
    ee[ , 2] = sapply(ee[ , 2], function(x) ifelse(x > tips[2], x + max(edge) - tips[2], x + tips[1]))
    edge = rbind(edge, c(sum(tips) + 1, max(edge) + 1), ee) # paste the 2 edges, adding the connection to the second phylo
    
    # creates the edge.length
    edge.length = c(phy1$root.edge, phy1$edge.length, phy2$root.edge, phy2$edge.length)
  }
  
  out = list("edge" = edge, "edge.length" = edge.length, "Nnode" = Nnode, "tip.label" = paste0("sp.", 1:sum(tips)))
  class(out) = "phylo"
  
  return(out)
}
#aa=twoL(try[[1]][[1]])
twoL = function(L){
  out = list(NA, NA)
  for(i in 1:2){
    ind = sort(open_tree(L, i))
    if(length(ind) > 0){
      out[[i]] = L[c(i, ind), ]
      sp = length(unique(out[[i]][ , 6]))
      or = order(out[[i]][ , 6])
      mul = table(out[[i]][ , 6])
      out[[i]][or, 6] = rep(1:sp, mul)
      if(i == 2){
        out[[2]][1, 1] = 1
        out[[2]][ , 2] = sapply(out[[2]][ , 2], function(x) ifelse(x == 2, 1, x))
      }
    } else{
      out[[i]] = matrix(L[i, ], nrow = 1)
    }
  }
  
  names(out) = c("L1", "L2")
  
  return(out)
}
#mat=L;ind=1
open_tree = function(mat, ind){
  out = which(mat[ , 2] == mat[ind, 1])
  if(length(out) > 0){
    for(i in 1:length(out)){
      aux2 = open_tree(out[i], mat = mat)
      out = c(out, aux2)
    }
  }
  return(out)
}
