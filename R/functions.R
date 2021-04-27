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
  raw <- readLines(treefile)
  raw <- gsub(" CSM", "", raw)
  raw <- gsub(" ", "_", raw)
  phy <- ape::read.tree(text=raw)
  phy$tip.label <- unname(GetAllGenusSpecies(phy$tip.label))
  phy <- root(phy, outgroup = "Notharctus tenebrosus")
  phy$edge.length = phy$edge
  phy <- multi2di(phy)
  return(phy)
}

#############################################################
# Import data and remove blank values of column of interest #
#############################################################
DataImport <- function(x, CleanColumn = NA){
  import <- read.csv(x, strip.white=TRUE)
  x <- import[!(import[,CleanColumn] ==""),]
  y <- unite(x, GenusSpecies, c(1:2), sep=" ", remove=FALSE)
  row.names(y) <- y$GenusSpecies
  return(y)
}

#For continuous data
ContinData <- function(dataimport){
  con <- subset(dataimport, select = c(11))
  con[,1] <- as.numeric(as.character(con[,1]))
  return(con)
}

#Ford disctrete data
DiscData <- function(dataimport){
  disc <- subset(dataimport, select = c(5))
  return(disc)
}

##############################
# Combine Tree and Dataframe #
##############################
CleanData <- function(phy, data) {
  DNIT <- geiger::name.check(phy, data)
  newphy <- geiger::drop.tip(phy, tip = DNIT$tree_not_data)
  #clean <- geiger::treedata(phy, data, warnings=TRUE)
  return(n)
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
                      c("Browning Motion",
                        "Ornstein-Uhlenbeck",
                        "Early-Bust"))
  table <- aic.w(aic.vals)
  return(list("tables" = table, "aic" = aic.vals, "BM" = BM1, "OU" = OU, "EB" = EB))
}

BestMod <- function(mods){
  if(names(which.max(mods$table)) == "Browning Motion"){
    print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$BM$opt$sigsq, "."))
  } else if(names(which.max(mods$table)) == "Ornstein-Uhlenbeck"){print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$OU$opt$sigsq, "."))
  } else if(names(which.max(mods$table)) == "Early-Bust"){print(paste0("The best model is ", names(which.max(mods$table)), " with its weighted aic value of ", max(mods$table), ". This model has an evolution rate of ", mods$EB$opt$sigsq, "."))
  }
}

TreeCompare <- function(datatree, mods){
  par(mfcol = (c(1,2)))
  plot(datatree$phy, show.tip.label=FALSE, main = "Original")
  if(names(which.max(mods$table)) == "Browning Motion"){
    bm.tree <- rescale(datatree$phy, model="BM")
    plot(ou.tree,show.tip.label=FALSE, main = "Browning Motion")
    } else if(names(which.max(mods$table)) == "Ornstein-Uhlenbeck"){
      ou.tree <- rescale(datatree$phy, model="OU", alpha = mods$OU$opt$alpha)
      plot(ou.tree, show.tip.label=FALSE, main = "Ornstein-Uhlenbeck")
      } else if(names(which.max(mods$table)) == "Early-Bust"){
        eb.tree <- rescale(datatree$phy, model="EB", alpha = mods$EB$opt$alpha)
        plot(ou.tree, show.tip.label=FALSE, main = "Early-Bust")
        }
  }







