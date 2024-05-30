multigroupSOM <- function(env, group.label.list) {
  require(oposSOM)
  # group.label.list - new group labels for SOM as list
  # env - SOM environment
  
  group.label.list <- lapply(group.label.list, function(x) 
  {names(x) <- colnames(env$indata) 
  return(x)
  } )
  
  if (is.null(names(group.label.list)) ) {
    names(group.label.list) <- LETTERS[length(group.label.list)]
  }
  
  for(i in 1:length(group.label.list)) {
    env$preferences$dataset.name = names(group.label.list)[i]
    env$files.name = names(group.label.list)[i]
    
    env$output.paths <- c(CSV = "CSV Sheets", `Summary Sheets Samples` = "Summary Sheets - Samples")
    
    dir.create(paste(env$files.name, "- Results"), showWarnings = FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), 
               showWarnings = FALSE)
    setwd(paste(env$files.name, "- Results"))
    env$group.labels <- group.label.list[[i]]
    
    
    o <- order(env$group.labels)
    env$pat.labels <- env$pat.labels[o]
    env$group.labels <- env$group.labels[o]
    env$p.g.m <- env$p.g.m[, o]
    env$fdr.g.m <- env$fdr.g.m[, o]
    env$n.0.m <- env$n.0.m[o]
    env$perc.DE.m <- env$perc.DE.m[o]
    env$p.m <- env$p.m[, o]
    env$indata <- env$indata[, o]
    env$metadata <- env$metadata[, o]
    env$indata.sample.mean <- env$indata.sample.mean[o]
    
    env$group.colors <- rep("#000000", ncol(env$indata))
    
    for (i in seq_along(unique(env$group.labels)))
    {
      env$group.colors[which(env$group.labels == unique(env$group.labels)[i])] <- color.palette.discrete(length(unique(env$group.labels)))[i]
    }
    
    if (length(unique(substr(env$group.colors, 1, 1)) > 1) || unique(substr(env$group.colors, 1, 1))[1] != "#")
    {
      env$group.colors <- apply(col2rgb(env$group.colors), 2, function(x) { rgb(x[1]/255, x[2]/255, x[3]/255) })
    }
    names(env$group.colors) <- colnames(env$indata)
    
    env$groupwise.group.colors <- env$group.colors[match(unique(env$group.labels), env$group.labels)]
    names(env$groupwise.group.colors) <- unique(env$group.labels)
    
    env$preferences$activated.modules$primary.analysis <- FALSE
    
    
    
    # filename <- paste(env$files.name, "pre.RData")
    # util.info("Saving environment image:", filename)
    # save(env, file = filename)
    util.info("Processing Differential Expression Statistics")
    # env <- pipeline.diffExpressionStatistics(env)
    util.info("Detecting Spots")
    #env <- pipeline.detectSpotsSamples(env) obsolete
    env <- pipeline.detectSpotsModules(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
    
    #Needed for saving the environment  
    # filename <- paste(env$files.name, ".RData", sep = "")
    # util.info("Saving environment image:", filename)
    # save(env, file = filename)
    
    # if (file.exists(paste(env$files.name, "pre.RData")) && 
    #     file.exists(filename)) {
    #   file.remove(paste(env$files.name, "pre.RData"))
    # }
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      util.info("Calculating Geneset Enrichment")
      env <- pipeline.genesetStatisticSamples(env)
      env <- pipeline.genesetStatisticModules(env)
    }
    
    util.info("Plotting Supporting Information")
    # pipeline.supportingMaps(env)
    # pipeline.entropyProfiles(env)
    # pipeline.topologyProfiles(env)
    
    if (ncol(env$indata) < 1000) {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    }
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      dir.create("Geneset Analysis", showWarnings=FALSE)
      
      util.info("Plotting Geneset Enrichment Heatmaps")
      pipeline.genesetOverviews(env)
      
      util.info("Plotting Geneset Profiles and Maps")
      pipeline.genesetProfilesAndMaps(env)
      
      util.info("Calculating Cancer Hallmark Enrichment")
      pipeline.cancerHallmarks(env)
    }
    
    util.info("Writing Gene Lists")
    pipeline.geneLists(env)
    util.info("Plotting Summary Sheets (Samples)")
    pipeline.summarySheetsSamples(env)
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    pipeline.summarySheetsPATs(env)
    
    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      pipeline.groupAnalysis(env)
    }
    
    util.info("Generating HTML Report")
    pipeline.htmlSampleSummary(env)
    pipeline.htmlModuleSummary(env)
    pipeline.htmlGenesetAnalysis(env)
    pipeline.htmlPsfAnalysis(env)
    pipeline.htmlSummary(env)
    setwd("..")
    
  }
}