plan <- drake_plan(
  tree_import = GetTreeWithNameProcessing("R/data/Primate_tree.tre"),
  data_import = DataImport("R/data/continuous_primate_data.csv",CleanColumn = 10), # Data of Ape biomass
  con_data = ContinData(data_import),
  disc_data = DiscData(data_import),
  DataTree = CleanData(tree_import, con_data),
  #pl = VisualizeData(DataTree$phy, DataTree$data, "out_data/continuous_tree.tree", "out_data/continuous_taxa.csv"),
  Mods = EvoMod(DataTree), #This has a list of the three modles, the aic output, and weight aic output
  Best = BestMod(Mods),
  CompareTrees = TreeCompare(DataTree, Mods)
  )




