rm(list=ls())

## Source all classes and packages ####

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Project Title: Expression Analysis for Landscape paper

## Instantiate a new Object of type ProjectSetUp ####
rnaseqProject <- ProjectSetUp$new(
  
  date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
  time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
  projectName             = "Glioma.data",
  annotationRDS           = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  pcRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  tfRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  ewsr1Fli1RDS            = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/EWSR1_FL1_DownstreamTargets.RDS",
  pax3Foxo1RDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",
  
  BrainExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Brain.v2.RDS",
  HeartExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Heart.v2.RDS", 
  KidneyExpRDS            = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Kidney.v2.RDS", 
  LiverExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Liver.v2.RDS", 
  LungExpRDS              = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Lung.v2.RDS", 
  
  outputPrefix            = "GliomaProject",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "MutationStatus.Wu",
  metadataFileRefCol      = "Sample.Biowulf.ID.GeneExp",
  metaDataFileName        = "MetadataMapper.v2.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  DiffGeneExpAnaDir       = "DiffExpResults",
  DiffGeneExpRDS          = "DiffGeneExpRDSOutput",
  ## Keep only Ribozero
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"),
  #                          "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                              "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  ## Keep only PolyA
  ## factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LibraryPrep" = "Ribozero"))
  ## Remove Celllines
  ## factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"))
  # factorsToExclude          = list('None'=list("LIBRARY_TYPE"=""))
)

## Add utility functions to the project ####
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

# Make a Tree Map ####
StatsFinal <-  rnaseqProject$metaDataDF %>% group_by_(.dots= c(rnaseqProject$factorName, "DIAGNOSIS.SubMut.Wu", "Color", 
                                                               "MutationStatus.Wu.color","LIBRARY_TYPE.TreeMap") ) %>% 
  count_(var=as.name("Sample.ID.Alias")) %>% dplyr::summarise(Count=n()) %>% 
  dplyr::group_by_(.dots= c(rnaseqProject$factorName)) %>%  
  dplyr::mutate( SampleSum := sum(Count)) %>% 
  spread_("LIBRARY_TYPE.TreeMap", "Count") %>% 
  mutate_( .dots = setNames( list( interp(~paste(rnaseqProject$factorName ,"(", Sum , ")"), 
                                          factorName=as.name(rnaseqProject$factorName), Sum=as.name("Sample") ) ), "LegendSampleSum") ) %>% 
  data.frame() %>% distinct(MutationStatus.Wu, DIAGNOSIS.SubMut.Wu,SampleSum, .keep_all = TRUE)

StatsFinal[,"LegendSampleSum"] <- paste(StatsFinal[,"MutationStatus.Wu"],"( ",StatsFinal[,"Sample"], " )",sep="")
StatsFinal

## Make the plot
pdf(file=paste(
    paste(rnaseqProject$workDir,rnaseqProject$projectName, rnaseqProject$plotsDir,"Mutation.Wu_Tree_Map_All", sep = "/"),
    "pdf",sep="."), height=8, width= 10)

#ggplot(StatsFinal, aes(area = Sample, fill = Color, label=MutationStatus.Wu, subgroup=DIAGNOSIS.SubMut.Wu)) +
ggplot(StatsFinal, aes(area = Sample, fill = MutationStatus.Wu.color, label=DIAGNOSIS.SubMut.Wu, subgroup=MutationStatus.Wu)) +
  geom_treemap(show.legend = TRUE) +
  geom_treemap_subgroup_border(colour="#4d4d4d") +
  geom_treemap_text(fontface = "bold.italic",
                    colour = "white",
                    place = "topleft",
                    padding.x = grid::unit(1.5, "mm"),
                    padding.y = grid::unit(1.5, "mm"),
                    grow = F,
                    reflow=T) +
  geom_treemap_subgroup_text(place = "bottom",
                             grow = T,
                             alpha = 0.5,
                             colour = "#FAFAFA",
                             min.size = 0,
                             padding.x = grid::unit(3, "mm"),
                             padding.y = grid::unit(1.5, "mm")) +
  scale_fill_identity()
dev.off()

## Generate expression matrix ####
rm(mergeObjectsNoDup)
# mergeObjectsNoDup <- corUtilsFuncs$getMergedMatrix(dir               = "TPM_Genes",
#                                                    fileFormat        = "txt",
#                                                    colNameSelect     = "expected_count",
#                                                    isRowNames        = TRUE,
#                                                    rowNamesColInFile = 1,
#                                                    fileSuffix        = ".rsem_ENS.genes.results",
#                                                    primaryID         = "gene_id",
#                                                    metadata          = rnaseqProject$metaDataDF,
#                                                    metadataFileRefCol= rnaseqProject$metadataFileRefCol )
# 
# saveRDS(mergeObjectsNoDup, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RawCount",paste0("mergeObjectsNoDup.rawcounts.",rnaseqProject$date,".rds"),sep="/") )

### Read DataSets ####
mergeObjectsNoDup_data <- readRDS(paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RawCount","mergeObjectsNoDup.rawcounts.2020-02-10.rds",sep="/"))

### Filter specific Histology samples ####
to_filter_by_histology = FALSE
if(to_filter_by_histology==TRUE){
  
} else {
  design <- rnaseqProject$metaDataDF
}

# ## Rearrange the design matrix and data matrix based on original Diagnosis
# design$DIAGNOSIS <- factor(design$DIAGNOSIS, levels =c("Astrocytoma", "Oligodendroglioma", "Glioblastoma", "Ependymoma",
#                                                        "Pleomorphic Xanthoastrocytoma", "Glioma", "Gliosarcoma",
#                                                        "Reccurent glioneuronal tumor"), ordered = TRUE)

# ## Rearrange the design matrix and data matrix based on Re-Assigned alias Diagnosis
design$MutationStatus.Wu <- factor(design$MutationStatus.Wu, levels = c("IDH_WT", "IDH1_Mutant", "IDH2_Mutant"), ordered = TRUE)
design$DIAGNOSIS <- factor(design$DIAGNOSIS.SubMut.Wu, levels =c("Anaplastic.PXA_III" , "Ependymoma_II" , "Ependymoma_Myxo" ,
                              "GBM_IV" , "GBM_IV_HMP" , "GN_II" , "GS_IV" , "PA" , "A_II" , "AA_III" , 
                              "AO_III" , "AO_III_HMP" , "O_II"), ordered = TRUE)

design %<>% arrange( MutationStatus.Wu,DIAGNOSIS )
print("design")
print(dim(design))
mergeObjectsNoDup <- mergeObjectsNoDup_data %>% dplyr::select(one_of(as.character(design[,rnaseqProject$metadataFileRefCol]))); 
print("data dim")
print(dim(mergeObjectsNoDup))

## Check if designmatrix and count matrix have same order of columns
View(data.frame(count_names=colnames(mergeObjectsNoDup), design_names=design[,rnaseqProject$metadataFileRefCol]))

## Evaluate presence of duplicate features (genes) and consolidate them ####
setDT(mergeObjectsNoDup, keep.rownames = TRUE)
mergeObjectsNoDup.pre <- mergeObjectsNoDup          %>% 
  dplyr::rename(GeneName = rn) 
mergeObjectsNoDup.pre <- dplyr::left_join(rnaseqProject$annotationDF[,c("GeneID", "GeneName")], mergeObjectsNoDup.pre, by="GeneName") %>% 
  data.table()
mergeObjectsConso     <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup.pre[,-c("GeneID")], funcName = "max", featureName = "GeneName")
mergeObjectsConso     <- dplyr::full_join(mergeObjectsConso, rnaseqProject$annotationDF[,c("GeneID", "GeneName")], by="GeneName") %>%  
  data.table()
mergeObjectsConso     <- subset(mergeObjectsConso,!duplicated(mergeObjectsConso$GeneName))
mergeObjectsConso     <- mergeObjectsConso[complete.cases(mergeObjectsConso), ]; dim(mergeObjectsConso)
mergeObjectsConso     <- mergeObjectsConso[,-c("GeneName")]         %>% 
  data.frame()                               %>% 
  tibble::column_to_rownames(var = "GeneID") %>% 
  as.matrix() ; dim(mergeObjectsConso)
## matching above data frame with the annotationDF
rnaseqProject$annotationDF <- rnaseqProject$annotationDF %>% dplyr::filter(GeneID %in% rownames(mergeObjectsConso)); dim(rnaseqProject$annotationDF)

# ## If TPM values needed then annotate the file and save file here
# mergeObjectsConsoTPM <- mergeObjectsConso %>% data.frame() %>%  tibble::rownames_to_column(var = "GeneID")
# mergeObjectsConsoTPM.Annot <- dplyr::left_join(rnaseqProject$annotationDF, mergeObjectsConsoTPM, by="GeneID")
# write.table(mergeObjectsConsoTPM.Annot, paste("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/",
#                                           paste0("mergeObjectsConsoTPM.Annot.TPM.annot",rnaseqProject$date,".txt"),sep="/"),
#             sep="\t", row.names = FALSE, quote = FALSE)

## Subset metaDataDF by the number of samples in the folder ####
colnamesDF           <- data.frame( "Sample.Biowulf.ID.GeneExp"= colnames(mergeObjectsConso))
corUtilsFuncs$subsetMetaData(colnamesDF=colnamesDF)

## Instantiate a new Object of type GeneExpNormalization ####
expressionObj        <- GeneExpNormalization$new(
  
  countObj          = as.matrix(mergeObjectsConso), 
  featureType       = "Gene", 
  packageRNAseq     = "edgeR", 
  annotationDF      = rnaseqProject$annotationDF, 
  design            = design[,rnaseqProject$factorName], 
  #design           = newMetaDataDF[,rnaseqProject$factorName],
  proteinCodingOnly = FALSE,
  corUtilsFuncs     = corUtilsFuncs
)

## Get expression in desired units ####
### RawCounts
#expressionTMM.Counts          = expressionObj$edgeRMethod("RawCounts")
## Normalised counts
#expressionTMM.NormDF         = expressionObj$edgeRMethod("NormFactorDF")

### RPKM
expressionTMM.RPKM            = expressionObj$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = FALSE)
designMatrix                  <- corUtilsFuncs$validfMatrix(df = design)
### Zscore ###
expressionTMM.RPKM.zscore <- expressionObj$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = TRUE)
## Replace NA with 0
expressionTMM.RPKM.zscore[is.na(expressionTMM.RPKM.zscore)] <- 0

## Arrange data by histology and Library type
expressionTMM.RPKM.arr <- expressionTMM.RPKM %>% dplyr::select(one_of("Chr","Start","End","Strand","GeneID","GeneName","Length",
                                                                      as.character(gsub("-",".",designMatrix[,rnaseqProject$metadataFileRefCol]))))

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM.arr)), 
                                                   designMatrix[,c(rnaseqProject$metadataFileRefCol,rnaseqProject$factorName,"Sample.ID.Alias", 
                                                                   "Sample.Data.ID", "MutationStatus.Wu", "DIAGNOSIS.SubMut.Wu", "Database.name")] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),"Database.name"])); 
View(AliasNames_df)

## Check if designmatrix and count matrix have same order of columns
View(data.frame(count_names=colnames(expressionTMM.RPKM.arr)[-c(1:7)], 
                design_names=designMatrix[,rnaseqProject$metadataFileRefCol],
                AliasColnames = as.character(AliasNames_df[-c(1:7),"Database.name"])))

## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM.arr)) == length(AliasColnames) )
colnames(expressionTMM.RPKM.arr)  <- AliasColnames
# 
### Save expression (TMM-RPKM/whatwever asked for in the above step) to a file ####
# write.table(expressionTMM.RPKM.arr, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
#                                           paste0("RPKM_Data_Filt_Consolidated.GeneNames.log2.",rnaseqProject$date,".txt"),sep="/"),
#             sep="\t", row.names = FALSE, quote = FALSE)

## Arrange data by histology and Library type (Zscore)
expressionTMM.RPKM.arr.zscore <- expressionTMM.RPKM.zscore  %>% dplyr::select(one_of("Chr","Start","End","Strand","GeneID","GeneName","Length",
                                                                  as.character(gsub("-",".",designMatrix[,rnaseqProject$metadataFileRefCol]))))

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM.arr.zscore)), 
                                                   designMatrix[,c(rnaseqProject$metadataFileRefCol,rnaseqProject$factorName,"Sample.ID.Alias", 
                                                                   "Sample.Data.ID","MutationStatus.Wu", "DIAGNOSIS.SubMut.Wu",  "Database.name")] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),"Database.name"])); AliasColnames

## Check if designmatrix and count matrix have same order of columns
View(data.frame(count_names=colnames(expressionTMM.RPKM.arr.zscore)[-c(1:7)], 
                design_names=designMatrix[,rnaseqProject$metadataFileRefCol],
                AliasColnames = as.character(AliasNames_df[-c(1:7),"Database.name"])))

## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM.arr.zscore)) == length(AliasColnames) )
colnames(expressionTMM.RPKM.arr.zscore)  <- AliasColnames

### Save expression (TMM-RPKM/whatwever asked for in the above step) to a file ####
#rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
# write.table(expressionTMM.RPKM.arr.zscore, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
#                                                  paste0("RPKM_Data_Filt_Consolidated.GeneNames.log2.zscore.",rnaseqProject$date,".txt"),sep="/"),
#             sep="\t", row.names = FALSE, quote = FALSE)
# 

### Perform hirarchial clustering of all the samples & Remove frozen samples

## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    code <- substr(label, 1, 1)
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=customColorsVector[label])
  }
  return(x)
}

## set colorCodes
hcColorPalette <- unique(as.character( rnaseqProject$metaDataDF[, "Color"] ))

## Consider all samples
RPKM_Data_Filt.zscore.hclust <- expressionTMM.RPKM.arr.zscore[, -c(1:7)];

## Remove frozen samples
# frozen_samples <- c( "CL0045_T2R.IDH2.AA",
#                      "CL0046_T2R.IDH1_HMP.AA",
#                      "CL0051_T2R.IDH1.DAC",
#                      "CL0052_T2R.WT.PA",
#                      "CL0033_T3R_T2.IDH1.AA",
#                      "CL0048_T2R.WT.GBM",
#                      "CL0053_T2R.WT.GBM",
#                      "CL0047_T2R.IDH1.AOG",
#                      "CL0049_T2R.IDH1_HMP.AOG")
RPKM_Data_Filt.zscore.hclust <- expressionTMM.RPKM.arr.zscore[, -c(1:7)]; dim(RPKM_Data_Filt.zscore.hclust)
# RPKM_Data_Filt.zscore.hclust <- RPKM_Data_Filt.zscore.hclust[, !(names(RPKM_Data_Filt.zscore.hclust) %in% frozen_samples)]; dim(RPKM_Data_Filt.zscore.hclust)
  
rownames(RPKM_Data_Filt.zscore.hclust) <-expressionTMM.RPKM.arr.zscore[, 6]; 
RPKM_Data_Filt.zscore.hclust.t <- t(RPKM_Data_Filt.zscore.hclust)

hc<-hclust(dist(RPKM_Data_Filt.zscore.hclust.t,"euclidean"),"ward.D")
clus2=cutree(hc,4)
#dend <- dendrapply(as.dendrogram(hc), labelCol)
op = par(bg="#FFFFFF")

## Make plot 
pdf(file=paste(
  paste(rnaseqProject$workDir,rnaseqProject$projectName, rnaseqProject$plotsDir,paste0("hirarchial_clustering_nofrozen.",
                                                                                       rnaseqProject$date), sep = "/"),
  "pdf",sep="."), height=15, width= 15)
plot(as.phylo(hc),type = "fan", tip.color=hcColorPalette[clus2],label.offset=1,col="red",cex=1)
dev.off()

### Performing ssGSEA output analysis. ( Plotting the scores across histology ) ##########

### Prepare input for ssGSEA broad gene pattern

expressionTMM.RPKM.GSEA.Input <- expressionTMM.RPKM.arr[, -c(1:7)];
rownames(expressionTMM.RPKM.GSEA.Input) <-expressionTMM.RPKM.arr[, 6]; 
expressionTMM.RPKM.GSEA.print = corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.GSEA.Input)

# # ## Save input for ssGSEA 
# write.table(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,"GSEA/rnk",
#                                                  paste0("RPKM_Data_Filt_Consolidated.GeneNames.log2.ssGSEA",rnaseqProject$date,".txt"),sep="/"),
#             sep="\t", row.names = FALSE, quote = FALSE)
# saveRDS(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.",rnaseqProject$date,".rds"),sep="/"))

## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../Glioma.data/GSEA/output/RPKM_Data_Filt_Consolidated.GeneNames.log2.ssGSEA2020-02-10.PROJ.gct")

## Add custom expression like cytolytic scre and HLA gene expression to the ssGSEA Outpuut file.
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
HLA_cytolyticScore      <- rbind(expressionTMM.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore)
View(data.frame(colnames(HLA_cytolyticScore), colnames(ssGSEAScores)))
HLA_cytolyticScore.ordered      <- HLA_cytolyticScore %>% tibble::rownames_to_column() %>% 
                                                    dplyr::select(one_of(c("rowname",colnames(ssGSEAScores))))
HLA_cytolyticScore.ordered %<>% tibble::column_to_rownames(var = "rowname")
ssGSEAScores.HLA.Cyto   <- rbind(ssGSEAScores,HLA_cytolyticScore.ordered)
dim(ssGSEAScores.HLA.Cyto)

## Plot the one variable plot
## Sanity Check: Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto) == length(as.character(rnaseqProject$validMetaDataDF$Database.name)) )

## Filter specified Diagnosis
factorsToExclude              = c("None")
selected.metadata              <- rnaseqProject$validMetaDataDF  %>% 
  filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
  dplyr::select_( .dots=c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName, "Database.name", "MutationStatus.Wu", "HMP_profile" ) )
dim(selected.metadata)

ssGSEAScores.HLA.Cyto.Selected <- ssGSEAScores.HLA.Cyto %>% 
                        dplyr::select(one_of(as.character(selected.metadata$Database.name)))
dim(ssGSEAScores.HLA.Cyto.Selected)

## sanity check Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto.Selected) == length(as.character(selected.metadata$Database.name)) )

## Preparing the expression matrix for string plot, by appending metadata
Scores <- cbind(t(ssGSEAScores.HLA.Cyto.Selected), selected.metadata[,c("MutationStatus.Wu","HMP_profile"), drop=FALSE]) %>% 
  dplyr::rename_(.dots = setNames( list("MutationStatus.Wu"), list("Diagnosis") )) #%>%
#dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))

### Setting up variables for  string plot
## Set the order of Diagnosis to appear
orderOfFactor    <- unique(Scores$Diagnosis)
## Set the order of signature to appear
orderOfSignature <- colnames(Scores)[1:(ncol(Scores)-2)]
## Total list of signatures
colList          <- c(1:(ncol(Scores)-2))
## Generate custom colors
customColorDF <- data.frame(Color=unique(as.character(rnaseqProject$metaDataDF$MutationStatus.Wu.color)) , 
                                     Diagnosis= unique(as.character(rnaseqProject$metaDataDF$MutationStatus.Wu)))
customColorDF

## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE,
                                                      yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF,
                                                      sizeOfDots=1.5, extradimension_colname="HMP_profile", extradimension_val="HMP")
plotLists[[1]]
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,
                               paste0("TMM-RPKM.ssGSEA.enrichmentScores.byVariantProfile.bean.log.", rnaseqProject$date, ".pdf"),sep="/")
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 10, height = 10)




