#########################
# Import and clean tree #
#########################

GetSingleGenusSpecies <- function(x) {
  return(paste(strsplit(x, " |_")[[1]][1:2], collapse=" "))
}

GetAllGenusSpecies <- function(x) {
  sapply(x, GetSingleGenusSpecies)
}

GetTreeWithNameProcessing <- function(treefile) {
  phy <- read.nexus(treefile)
  phy$tip.label <- unname(GetAllGenusSpecies(phy$tip.label))
  phy$tip.label <- make.unique(phy$tip.label)
  phy <- multi2di(phy)
  phy$edge.length[phy$edge.length==0]<-1e-16
  return(phy)
}

#############################################################
# Import data and remove blank values of column of interest #
#############################################################

find_unique <-  function(diet_data){
  x <- unique(diet_data)
  return(x)
} 

clean_up_diet <- function(unique_data){
  unique_data$diet_binary <- ifelse(unique_data$diet == "insectivore, omnivore" | unique_data$diet == "insectivore, herbivore" |
                                      unique_data$diet == "insectivore, frugivore" | unique_data$diet =="insectivore, carnivore" |
                                      unique_data$diet == "insectivore" | unique_data$diet == "omnivore", 1, 0)
  #rownames(unique_data) <- unique_data[, "GenusSpecies"]
  return(unique_data)
}

DataImport <- function(data1, data2, col1="", col2="", CleanColumn = ""){
  import1 <- read.csv(data1, strip.white=TRUE) # Import continuous data
  import2 <- read.csv(data2, strip.white=TRUE) # Import discrete data
  subset <- import2[,c(col1, col2)] # Subset discrete data
  UniqueIm2 <- find_unique(subset)
  CleanIm2 <- clean_up_diet(UniqueIm2)
  x <- import1[!(import1[,CleanColumn] ==""),]
  y <- unite(x, GenusSpecies, c(1:2), sep=" ", remove=FALSE) # Merge genus and species name into one col of continuous data
  merge <- CleanIm2 %>% inner_join(., y, by = col1)
  row.names(merge) <- merge$GenusSpecies
  return(merge)
}


##############################
# Combine Tree and Dataframe #
##############################
CleanData <- function(phy, data) {
  clean <- geiger::treedata(phy, data, sort = TRUE, warnings=TRUE)
  return(clean)
}

################
# Making trees #
################
plot_tree <- function(tree, file, label) {
  pdf(file=file)
  plot(tree, type = "fan", show.tip.label=label, edge.width = 0.1, cex = .5)
  dev.off()
}

VisualizeData <- function(phy, data, file_1, file_2) {
  
  plot(phy, show.tip.label=TRUE, edge.width = 0.1, cex = .25)
  newdata <- data
  print(newdata)
  
  plot_tree(phy, file_1, TRUE)
  write.csv(newdata, file = file_2)
  
}

#################################
# Building and Comparing Models #
#################################
EvoMod <- function(datatree){
  BM1 <- geiger::fitContinuous(datatree$phy, datatree$data, model="BM", ncores = 1)
  OU  <- geiger::fitContinuous(datatree$phy, datatree$data, model="OU", ncores = 1)
  EB  <- geiger::fitContinuous(datatree$phy, datatree$data, model="EB", ncores = 1)
  aic.vals <-setNames(c(BM1$opt$aicc,
                        OU$opt$aicc,
                        EB$opt$aicc),
                      c("Brownian Motion",
                        "Ornstein-Uhlenbeck",
                        "Early-Bust"))
  table <- aic.w(aic.vals)
  return(list("tables" = table, "aic" = aic.vals, "BM" = BM1, "OU" = OU, "EB" = EB))
}

BestMod <- function(mods){
  if(names(which.max(mods$table)) == "Brownian Motion"){
    print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$BM$opt$sigsq, "."))
  } else if(names(which.max(mods$table)) == "Ornstein-Uhlenbeck"){print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$OU$opt$sigsq, "."))
  } else if(names(which.max(mods$table)) == "Early-Bust"){print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$EB$opt$sigsq, "."))
  }
}

TreeCompare <- function(datatree, mods){
  par(mfcol = (c(1,2)))
  plot(datatree$phy, show.tip.label=FALSE, main = "Original")
  if(names(which.max(mods$table)) == "Brownian Motion"){
    transformed.tree <- geiger::rescale(datatree$phy, model="BM", sigsq = mods$BM$opt$sigsq)
  } else if(names(which.max(mods$table)) == "Ornstein-Uhlenbeck"){
    transformed.tree <- geiger::rescale(datatree$phy, model="OU", alpha = mods$OU$opt$alpha)
  } else if(names(which.max(mods$table)) == "Early-Bust"){
    transformed.tree <- geiger::rescale(datatree$phy, model="EB", a = mods$EB$opt$a)
  }
  plot(transformed.tree, show.tip.label=FALSE, main = names(which.max(mods$table)))
}

##################
# OUwie Analysis #
##################
RegimeAssignment <- function(discdatatree, datatree){
  #label tree with best states
  reconstruction.info <- ace(discdatatree$data, datatree$phy, type="discrete", method="ML", CI=TRUE)
  best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]
  labeledtree <- datatree$phy
  labeledtree$node.label <- best.states
  
  #Continued analysis
  data <- as.data.frame(discdatatree$data)
  names<- as.data.frame(tibble::rownames_to_column(data, "SPECIES"))
  colnames(names)[2] <- "Species"
  regime = as.data.frame(datatree$data)
  names2<- as.data.frame(tibble::rownames_to_column(regime, "SPECIES"))
  colnames(names2)[2] <- "Species"
  temp2 <- data.frame("species"= names[,1], "regime"=names[,2],
                      "continuous"=names2[,2])
  return(list("temp2" = temp2, "labeledtree" = labeledtree))
}
 # I wanted to place this in my plan, but drake kept skipping it and 
                                                           # moving to the next step, which caused a fault.
AllOuwie <- function(phy, model, data){
  results <- lapply(model, OUwie, phy = phy, data = data)
  AICc.values<-sapply(results, "[[", "AICc")
  names(AICc.values)<-model
  AICc.values<-AICc.values-min(AICc.values)
  best<-results[[which.min(AICc.values)]] #store for later
  return(list("AIC_values" = AICc.values, "best" = best)) #prints info on best model
}

