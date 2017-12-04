
get_SpecExt4birthdeath = function(turnover, divers = NULL){
  # turnover = d/b
  # divers = b-d
  
  if(is.null(divers)){
    divers = turnover[2]
    turnover = turnover[1]
  }
  
  birth = divers / (1 - turnover)
  death = birth - divers
  
  return(structure(c(birth, death), names = c("birth", "death")))
}

unlist_branches = function(bran){
  out = list()
  ll = length(bran)
  ll1 = c(1, cumsum(sapply(bran, length))[-ll] + 1)
  ll2 = cumsum(sapply(bran, length))
  for(i in 1:ll){
    out[ (ll1[i]):(ll2[i]) ] = bran[[i]]
  }
  return(out)
}

table_prop_phylo = function(vec, categories = NULL){
  test = names(attributes(vec))
  if("harshness" %in% test){
    aux = attr(vec, "harshness")
  } else if ("latitude" %in% test){
    aux = attr(vec, "latitude")
  } else{
    stop("'vec' must have an attribute called either 'harshness' or 'latitude'.")
  }
  
  ll = length(vec)
  if(is.null(categories)){
    out = table(aux) / ll
  } else{
    out = sapply(categories, function(x) sum(aux == x) / ll)
  }
  
  return(out)
}

bin_harshness = function(x, bins){
  # Function to create a discrete categorization of 'ENV.HARSHNESS' from Botero14
  
  if(bins == 5){
    if(x < -1.5){return("ExtStable")}
    if(x < -0.5){return("Stable")}
    if(x < 0.5){return("Regular")}
    if(x < 1.5){return("Variable")}
    return("ExtVariable")
  } else if(bins == 3){
    if(x < -1.0){return("Stable")}
    if(x < 1.0){return("Regular")}
    return("Variable")
  } else{
    stop("'bins' must be 3 or 5.")
  }
}

get_Botero_spp = function(these, col.sp = "SPECIES", col.subsp = "SUBSPECIES", ref = birds){
  # Gets th number of subspecies per species (default set to Botero14)
  #   - these, character vector with the species of interest
  #   - col.sp, the column witht the species names
  #   - col.subsp, the column with subspecies info
  #   - ref, the reference dataset. It MUST be a data.frame
  if(class(these) == "phylo"){
    these = these$tip.label
  }
  if(class(these) != "character"){
    stop("'these' must be either a character vector or a phylogeny with the species names.")
  }
  if(class(ref) != "data.frame"){
    stop("'ref' must be a data.frame")
  }
  
  ind = match(these, ref[[col.sp]])
  if(any(is.na(ind))){
    warning("Some species names were not found in the 'ref' and will be excluded from analysis!")
    ind = na.exclude(ind)
  }
  
  out = ref[ind, col.subsp]
  
  return(out)
}

get_prior_exp = function(subspp, output = "function", log = TRUE){
  #   Given a vector with the number of subspecies ('subspp') it estimates the 
  # best rate parameter, with ML. It then returns either the fitted model
  # (output = "fit"), the rate parameter (output = "coef"), or the prior function 
  # for the b1 parameter to be used in pbd_Bayes (output = "function"; the default).
  LLexp = function(r) {
    out = suppressWarnings(dexp(x = subspp, rate = r, log = TRUE))
    -sum(out)
  }
  par = stats4:::mle(LLexp, start = list(r = .1))
  if(output == "fit"){
    return(par)
  }
  if(output == "coef"){
    return(par@coef)
  }
  priorb1exp = function(x) dexp(x, rate = par@coef, log = log)
  return(priorb1exp)
}

# tree = phy = rcoal(10);plot(phy)
# tipLabels = sample(c("Temperate", "Tropical", "Mixed"), 10, TRUE)
get_matching_clades = function(phy, clade.size, tipLabels, minimum = TRUE, mixed = 1, multiLat = 0){
  # Finds all subclades from the phylogenies that meet certain criteria:
  #   - phy, phylogeny that'll serve as base
  #   - clade.size, minimum (or maximum) clade size
  #   - tipLabels, ordered (to match 'phy' tip.label) character vector indicating the latitudinal distribution of the species, one of: "Temperate", "Tropical", "Mixed"
  #   - minimum = TRUE, is 'clade.size' the minimum species number desired?
  #   - mixed = 1, proportion of "Mixed" species permited in the subclades
  #   - multiLat = 0, proportion of 'Tropical' species allowed in 'Temperate' subclades (and vice-versa)
  #
  # Returns one of the three: 
  # - NULL, if there is NO clade that meets 'clade.size' (and 'minimum') parameters
  # - NA, if there is NO clade that meets 'mixed' and 'multiLat' parameters
  # - a list of numeric vectors indicating the tips that meet the criteria. Each vector has 2 attributes: "species", with the corresponding tip.label; and "latitude", with the latitudinal distribution ("Temperate", "Tropical", "Mixed")
  
  if(mixed < 0 | mixed > 1){
    stop("'mixed' indicates proportion of 'Mixed' species in the subclades,\nso it must be between 0 and 1!")
  }
  if(multiLat < 0 | multiLat > 1){
    stop("'multiLat' indicates proportion of 'Tropical' species allowed in 'Temperate' subclades (and vice-versa),\nso it must be between 0 and 1!")
  }
  
  # get species names from phylogeny
  real.nms = phy$tip.label
  # assigns the latitudinal distribution to the tip.labels
  phy$tip.label = tipLabels
  
  # function to get all subclades from a phylogeny
  foo = ape:::prop.part(phy)
  # 'prop.part' returns numeric vectors with the tip numbers and tip.labels as an attribute
  nms = attr(foo, which = "labels")
  
  # filters and tests if any subclade meets 'clade.size' criteria
  if(minimum){
    foo2 = foo[sapply(foo, length) >= clade.size]
  } else{
    foo2 = foo[sapply(foo, length) <= clade.size]
  }
  if(length(foo2) == 0){
    return(NULL)
  }
  
  # gets the latitudinal distribution for all subclades of interest
  these.nms = lapply(foo2, function(x) nms[x])
  # gets the species names for all subclades of interest
  sp.nms = lapply(foo2, function(x) real.nms[x])
  # assigns "species" and "latitude" as attributes
  for(i in 1:length(foo2)){
    attr(foo2[[i]], which = "latitude") = these.nms[[i]]
    attr(foo2[[i]], which = "species") = sp.nms[[i]]
  }
  
  # creates a function to evaluate the subclades according to 'mixed' and 'multiLat' parameters
  if(multiLat == 0){
    fun = function(x){
      bool1 = ("Tropical" %in% x) != ("Temperate" %in% x)
      bool2 = (sum(x == "Mixed") / length(x)) <= mixed
      return(bool1 & bool2)
    }
  } else if(multiLat == 1){
    fun = function(x){
      bool = (sum(x == "Mixed") / length(x)) <= mixed
      return(bool)
    }
  } else{
    fun = function(x){
      p = c(sum(x == "Tropical"),
            sum(x == "Temperate"),
            sum(x == "Mixed")
      ) / length(x)
      bool1 = (p[1] <= multiLat) | (p[2] <= multiLat)
      bool2 = p[3] <= mixed
      return(bool1 & bool2)
    }
  }
  
  # filters and tests if any subclade follows the 'mixed' and 'multiLat' parameters
  bool = sapply(these.nms, FUN = fun)
  if(sum(bool) == 0){
    return(NA)
  }
  out = foo2[bool]
  
  return(out)
}


# tree = phy = rcoal(20);plot(phy)
# (tipLabels = sample(c("Variable", "Regular", "Stable"), 20, TRUE))
# clade.size = 3; mixed = 0.7
get_matching_cladesHarsh = function(phy, clade.size, tipLabels, minimum = TRUE, mixed = 1){
  # Finds all subclades from the phylogenies that meet certain criteria:
  #   - phy, phylogeny that'll serve as base
  #   - clade.size, minimum (or maximum) clade size
  #   - tipLabels, ordered (to match 'phy' tip.label) character vector indicating the latitudinal distribution of the species, one of: "Variable", "Regular", "Stable"
  #   - minimum = TRUE, is 'clade.size' the minimum species number desired?
  #   - mixed = 1, proportion of "Variable", "Regular" species allowed in "Stable" subclades (and all other combinations)
  #
  # Returns one of the three: 
  # - NULL, if there is NO clade that meets 'clade.size' (and 'minimum') parameters
  # - NA, if there is NO clade that meets 'mixed' parameter
  # - a list of numeric vectors indicating the tips that meet the criteria. Each vector has 2 attributes: "species", with the corresponding tip.label; and "harshness", with the harshness classification ("Variable", "Regular", "Stable")
  
  if(mixed < 0 | mixed > 1){
    stop("'mixed' indicates proportion of mix between different categories in the subclades,\nso it must be between 0 and 1!")
  }
  
  # get species names from phylogeny
  real.nms = phy$tip.label
  # assigns the latitudinal distribution to the tip.labels
  phy$tip.label = tipLabels
  
  # function to get all subclades from a phylogeny
  foo = ape:::prop.part(phy)
  # 'prop.part' returns numeric vectors with the tip numbers and tip.labels as an attribute
  nms = attr(foo, which = "labels")
  
  # filters and tests if any subclade meets 'clade.size' criteria
  if(minimum){
    foo2 = foo[sapply(foo, length) >= clade.size]
  } else{
    foo2 = foo[sapply(foo, length) <= clade.size]
  }
  if(length(foo2) == 0){
    return(NULL)
  }
  
  # gets the latitudinal distribution for all subclades of interest
  these.nms = lapply(foo2, function(x) nms[x])
  # gets the species names for all subclades of interest
  sp.nms = lapply(foo2, function(x) real.nms[x])
  # assigns "species" and "latitude" as attributes
  for(i in 1:length(foo2)){
    attr(foo2[[i]], which = "harshness") = these.nms[[i]]
    attr(foo2[[i]], which = "species") = sp.nms[[i]]
  }
  
  # creates a function to evaluate the subclades according to 'mixed' parameter
  if(mixed == 0){
    fun = function(x){
      bool1 = length(unique(x)) == 1
      return(bool1)
    }
  } else{
    fun = function(x){
      p = sapply(unique(x), FUN = function(y) sum(x == y) / length(x))
      bool1 = (max(p) + mixed) >= 1
      return(bool1)
    }
  }
  
  # filters and tests if any subclade follows the 'mixed' parameters
  bool = sapply(these.nms, FUN = fun)
  if(sum(bool) == 0){
    return(NA)
  }
  out = foo2[bool]
  
  return(out)
}

##### Functions to analyse how the different taxonomic categories are distributed
get_prop = function(column, type = NULL, ref = birds){
  if(is.null(type)){
    stop("'type' must be either 'harshness' or 'latitude'.")
  }
  
  if(type == "harshness"){
    out = get_prop_harshness(column, ref)
  } else if(type == "latitude"){
    out = get_prop_latitude(column, ref)
  } else{
    stop("'type' must be either 'harshness' or 'latitude'.")
  }
  
  return(out)
}
plot_prop = function(out, type = NULL, minimum = 1){
  if(is.null(type) | !(type %in% c('harshness', 'latitude'))){
    stop("'type' must be either 'harshness' or 'latitude'.")
  }
  
  if(type == "harshness"){
    gg = plot_prop_harshness(out, minimum)
  } else{ # if(type == "latitude")
    gg = plot_prop_latitude(out, minimum)
  }
  
  return(gg)
}
get_prop_harshness = function(column, ref = birds){
  # estimates the number of species that are: "ExtVariable", "Variable", "Regular", "Stable", or "ExtStable"
  #   - column, numeric or character indicating the column of the taxonomy level to be used
  #   - ref, data.frame with the information about harshness MUST have a column named "harshness"
  
  require(dplyr)
  ref = ref[!is.na(ref[[column]]), ]
  groups = unique(ref[[column]])
  if(class(column) == "numeric"){
    column = colnames(ref)[column]
  }
  
  out = data.frame(matrix(data = NA, nrow = length(groups), ncol = 6,
                          dimnames = list(groups, c("ExtVariable", "Variable", "Regular", "Stable", "ExtStable", "N"))))
  
  for(i in 1:length(groups)){
    this_row = ref %>% 
      filter_(paste(column, "==", shQuote(groups[i]))) %>%
      summarise(ExtVariable = sum(harshness == "ExtVariable"),
                Variable = sum(harshness == "Variable"),
                Regular = sum(harshness == "Regular"),
                Stable = sum(harshness == "Stable"),
                ExtStable = sum(harshness == "ExtStable"),
                N = length(harshness))
    if(sum(this_row[1:5]) == this_row[6]){
      out[i, ] = as.numeric(this_row)
    } else{
      warning(paste0("Problems talling the latitudinal distribution for group ", groups[i]))
    }
  }
  
  return(out)
}
plot_prop_harshness = function(out, minimum = 1){
  # creates a stacked bar plot from the output of "get_prop_latitude"
  #   - minimum, is the minimum number of species a group must have to be plotted
  require(ggplot2)
  require(reshape2)
  out$taxa = rownames(out)
  
  gg = out %>% 
    filter(N >= minimum) %>%
    mutate(ExtVariable = 100 * (ExtVariable / N),
           Variable = 100 * (Variable / N),
           Regular = 100 * (Regular / N),
           Stable = 100 * (Stable / N),
           ExtStable = 100 * (ExtStable / N)) %>%
    arrange(desc(N)) %>%
    mutate(taxa = factor(taxa, levels = taxa)) %>%
    melt(id.var = 6:7) %>%
    ggplot(aes(x = taxa, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(y = -1, label = N), size = 3)
  
  return(gg)
}
get_prop_latitude = function(column, ref = birds){
  # estimates the number of species that are "Temperate", "Tropical", or "Mixed"
  #   - column, numeric or character indicating the column of the taxonomy level to be used
  #   - ref, data.frame with the information about latitudinal distribution. MUST have a column named "LAT.RANGE"
  
  require(dplyr)
  ref = ref[!is.na(ref[[column]]), ]
  groups = unique(ref[[column]])
  
  out = data.frame(matrix(data = NA, nrow = length(groups), ncol = 4,
                          dimnames = list(groups, c("Temperate", "Tropical", "Mixed", "N"))))
  
  if(class(column) == "numeric"){
    column = colnames(ref)[column]
  }
  
  for(i in 1:length(groups)){
    this_row = ref %>% 
      filter_(paste(column, "==", shQuote(groups[i]))) %>%
      summarise(Temp = sum(LAT.RANGE == "Temperate"),
                Trop = sum(LAT.RANGE == "Tropical"),
                Mixed = sum(LAT.RANGE == "Mixed"),
                N = length(LAT.RANGE))
    if(sum(this_row[1:3]) == this_row[4]){
      out[i, ] = as.numeric(this_row)
    } else{
      warning(paste0("Problems talling the latitudinal distribution for group ", groups[i]))
    }
  }
  
  return(out)
}
plot_prop_latitude = function(out, minimum = 1){
  # creates a stacked bar plot from the output of "get_prop_latitude"
  #   - minimum, is the minimum number of species a group must have to be plotted
  require(ggplot2)
  require(reshape2)
  out$taxa = rownames(out)
  
  gg = out %>% 
    filter(N >= minimum) %>%
    mutate(Temperate = 100 * (Temperate / N),
           Tropical = 100 * (Tropical / N),
           Mixed = 100  * (Mixed / N) ) %>%
    arrange(desc(N)) %>%
    mutate(taxa = factor(taxa, levels = taxa)) %>%
    melt(id.var = 4:5) %>%
    ggplot(aes(x = taxa, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(y = -1, label = N), size = 3)
  
  return(gg)
}


correct_name = function(vec){
  # transforms the given names first letter into uppercase and all other into lowercase
  # just because.....
  
  out = strsplit(vec, split = "")[[1]]
  paste0(toupper(out[1]), paste0(tolower(out[-1]), collapse = ""))
}