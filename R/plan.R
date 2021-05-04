plan <- drake_plan(
  tree_import = GetTreeWithNameProcessing("R/data/consensusTree_10kTrees_Primates_Version3.nex"),
  data_import = DataImport(data1 ="R/data/continuous_primate_data.csv",
                           data2 = "R/data/primate_data.csv",
                           col1 = "GenusSpecies",
                           col2 = "diet",
                           CleanColumn = 10 
                           ), # Data of Ape biomass and diet
  DataTree = CleanData(tree_import, data_import["BodySize_g"]/1000), #Brian said divide by 1000 to get OUwie to work...and that worked...so here we are
  discDataTree = CleanData(tree_import, data_import["diet_binary"]),
  Mods = EvoMod(DataTree), #This has a list of the three modles, the aic output, and weight aic output
  Best = BestMod(Mods),
  CompareTrees = TreeCompare(DataTree, Mods),
  Contin = RegimeAssignment(discdatatree = discDataTree, datatree = DataTree),
  ouwie_analysis = OUwie(Contin$labeledtree, Contin$temp2, model="OUMV", simmap.tree=FALSE, diagn=FALSE),
  # Based on the values of the Rates, sigma (BM) seems to have the highest values. Further, it seems that
  #there is more selection for an herbavor lifestyle (0) than an ominavor one (1). 
  models = c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),
  FullOuwieAnalysis = AllOuwie(phy = Contin$labeledtree,  data = Contin$temp2, model = models)
  )
