rm(list=ls())

library(grid)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
setwd("T:/Sivasish_Sindiri/Collaboration/Jack")

## Memo sort
memoSort = function(M = NA) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

## Load variant dataset
variantData <- read.csv("./COG.variant.data.txt", sep="\t", header = T,row.names = 1 )
variantData[is.na(variantData)] <- 0;  dim(variantData)
variantDataE <- variantData %>% dplyr::select(-one_of("PKN1","IGF2","ALK","SOS1","ROBO1","PDGFRA","GAB1","CCND1","CCND2","SOS2","AKT1","SMARCA4")); dim(variantDataE)


## Load metadataset 
metadata <- read.csv("./COG_metadata.txt", sep="\t", header = T, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE); head(metadata)

# View(data.frame(data.col=colnames(variantData_memo), metadata.col=rownames(metadata_reorder)))

######################################################### Version 1 ###########################################
variantData_memo <- memoSort(t(variantDataE));
metadata_reorder <- metadata[colnames(variantData_memo),];

## Adding annotation
histology_col = structure(names = c("ARMS", "ERMS","MIXED A & E RMS", "PRMS", "RMS NOS", "SC-RMS"),
                          c("#6699ff", "#ff944d", "#3cb44b", "#66ffff", "#cc00ff", "#ffff66"))
EnrollmentAge_col = structure(names = c("< 1yr", "1-5yr","5-10yr", "10-15yr", "> 15yr", "unknown"),
                          c("#6699ff", "#33cc33", "#ff6600", "#ffbf00", "#cc66ff", "#e6e6e6"))
Sex_col = structure(names = c("Female","Male"),
                     c("#6699ff", "#e6e6e6"))
Cohort_col = structure(names = c("COG", "UK"),
                    c("#e6e6e6","#ff6600"))
FusionStatus_col  = structure(names= c("FN", "FP"), c("#e6e6e6","#6699ff"))
Anatomic_col = structure(names=
                         c("Bladder_Prostate","Extremity", "Female GU", "Head_Neck", "Orbital", "Other", "Parameningeal",  
                         "Paratesticular", "RPT", "Unknown"),
                         c("#911eb4", "#ffe119", "#ffd8b1", "#ff66b3", "#42d4f4", "#a9a9a9", "#F8636B", "#6699ff", 
                                    "#3cb44b", "#e6e6e6"))
Riskgroup_col = structure(names=c("High", "Intermediate", "Low", "Unknown"),c("#F8636B", "#F8D663", "#6699ff","#e6e6e6"))


top_ha = HeatmapAnnotation(
  Cohort = metadata_reorder[[1]],
  Sex = metadata_reorder[[3]],
  EnrollmentAge = metadata_reorder[[10]],
  Anatomic = metadata_reorder[[8]],
  Riskgroup = metadata_reorder[[9]],
  Histology = metadata_reorder[[5]],
  FusionStatus = metadata_reorder[[7]],
  col = list(
    Histology  = histology_col,
    FusionStatus = FusionStatus_col,
    Cohort = Cohort_col,
    Sex = Sex_col,
    EnrollmentAge = EnrollmentAge_col,
    Anatomic = Anatomic_col,
    Riskgroup = Riskgroup_col
  ),
  show_legend = c(TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Histology  = list(title = "Histology"),
    FusionStatus = list(title = "Fusion Status"),
    Cohort = list(title= "Cohort"),
    Sex = list(title="Sex"),
    EnrollmentAge = list(title="Age at enrollment"),
    Anatomic = list(title="Anatomic Group"),
    Riskgroup = list(title="Risk Group")
  )
)


pdf("plot5u.pdf", height = 8, width = 19)
## Displayinh Heatmap

color_list = c('#F4F4F4', "#0073e6")
#color_list = c('#F4F4F4', "#F86262")
Heatmap(variantData_memo, 
        name = "mat", 
        col = color_list, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        rect_gp = gpar(col= "#e6e6e6"),
        #rect_gp = gpar(col= NA),
        show_column_names = FALSE,
        top_annotation = top_ha,
        row_names_gp =  gpar(fontsize = 10),
        row_names_side = c("left"),
        show_heatmap_legend = FALSE,
        width = unit(390, "mm"),
    )

annotation_titles = c(Cohort = "Cohort",
                      Sex = "Sex",
                      EnrollmentAge = "Age",
                      Histology = "Histology",
                      FusionStatus = "Fusion Status",
                      Anatomic = "Anatomic Group",
                      Riskgroup = "Risk Group"
                    )

for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

decorate_annotation("Cohort", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Sex", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("EnrollmentAge", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Histology", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("FusionStatus", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Anatomic", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Riskgroup", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})

dev.off()
# 


######################################################### Version 2 ###########################################

## Adding annotation
histology_col = structure(names = c("ARMS", "ERMS","MIXED A & E RMS", "PRMS", "RMS NOS", "SC-RMS"),
                          c("#6699ff", "#ff944d", "#3cb44b", "#66ffff", "#cc00ff", "#ffff66"))
EnrollmentAge_col = structure(names = c("< 1yr", "1-5yr","5-10yr", "10-15yr", "> 15yr", "unknown"),
                              c("#6699ff", "#33cc33", "#ff6600", "#ffbf00", "#cc66ff", "#e6e6e6"))
Sex_col = structure(names = c("Female","Male"),
                    c("#6699ff", "#e6e6e6"))
Cohort_col = structure(names = c("COG", "UK"),
                       c("#e6e6e6","#ff6600"))
FusionStatus_col  = structure(names= c("FN", "FP"), c("#e6e6e6","#6699ff"))
Anatomic_col = structure(names=
                           c("Bladder_Prostate","Extremity", "Female GU", "Head_Neck", "Orbital", "Other", "Parameningeal",  
                             "Paratesticular", "RPT", "Unknown"),
                         c("#911eb4", "#ffe119", "#ffd8b1", "#ff66b3", "#42d4f4", "#a9a9a9", "#F8636B", "#6699ff", 
                           "#3cb44b", "#e6e6e6"))
Riskgroup_col = structure(names=c("High", "Intermediate", "Low", "Unknown"),c("#F8636B", "#F8D663", "#6699ff","#e6e6e6"))
color_list = c('#F4F4F4', "#0073e6")
#color_list = c('#F4F4F4', "#F86262")

### Split the data
## COG & FP
metadata_COG_FP     <- metadata %>% 
  tibble::rownames_to_column(var = 'SampleID') %>% 
  dplyr::filter( Cohort == "COG" & `Grouping FP or FN` == "FP") %>% 
  tibble::column_to_rownames(var='SampleID');
variantDataE_COG_FP_memo <- variantDataE[rownames(metadata_COG_FP), ] %>% t() %>% memoSort()
metadata_COG_FP_reorder     <- metadata_COG_FP[colnames(variantDataE_COG_FP_memo), ]
if(!identical(colnames(variantDataE_COG_FP_memo), rownames(metadata_COG_FP_reorder))){
  print( "Please check the metadata and variant data")
}


top_ha_COG_FP = HeatmapAnnotation(
  Cohort = metadata_COG_FP_reorder[[1]],
  Sex = metadata_COG_FP_reorder[[3]],
  EnrollmentAge = metadata_COG_FP_reorder[[10]],
  Anatomic = metadata_COG_FP_reorder[[8]],
  Riskgroup = metadata_COG_FP_reorder[[9]],
  Histology = metadata_COG_FP_reorder[[5]],
  FusionStatus = metadata_COG_FP_reorder[[7]],
  col = list(
    Histology  = histology_col,
    FusionStatus = FusionStatus_col,
    Cohort = Cohort_col,
    Sex = Sex_col,
    EnrollmentAge = EnrollmentAge_col,
    Anatomic = Anatomic_col,
    Riskgroup = Riskgroup_col
  ),
  show_legend = c(TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Histology  = list(title = "Histology"),
    FusionStatus = list(title = "Fusion Status"),
    Cohort = list(title= "Cohort"),
    Sex = list(title="Sex"),
    EnrollmentAge = list(title="Age at enrollment"),
    Anatomic = list(title="Anatomic Group"),
    Riskgroup = list(title="Risk Group")
  )
)

hm_COG_FP = Heatmap(variantDataE_COG_FP_memo, 
                    name = "mat", 
                    col = color_list, 
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE, 
                    rect_gp = gpar(col= "#e6e6e6"),
                    #rect_gp = gpar(col= NA),
                    show_column_names = FALSE,
                    top_annotation = top_ha_COG_FP,
                    row_names_gp =  gpar(fontsize = 10),
                    row_names_side = c("left"),
                    show_heatmap_legend = FALSE,
                    width = unit(100, "mm"),
)


## COG & FN
metadata_COG_FN     <- metadata %>% 
  tibble::rownames_to_column(var = 'SampleID') %>% 
  dplyr::filter( Cohort == "COG" & `Grouping FP or FN` == "FN") %>% 
  tibble::column_to_rownames(var='SampleID');
variantDataE_COG_FN_memo <- variantDataE[rownames(metadata_COG_FN), ] %>% t() %>% memoSort()
metadata_COG_FN_reorder     <- metadata_COG_FN[colnames(variantDataE_COG_FN_memo), ]
if(!identical(colnames(variantDataE_COG_FN_memo), rownames(metadata_COG_FN_reorder))){
  print( "Please check the metadata and variant data")
}


top_ha_COG_FP = HeatmapAnnotation(
  Cohort = metadata_COG_FN_reorder[[1]],
  Sex = metadata_COG_FN_reorder[[3]],
  EnrollmentAge = metadata_COG_FN_reorder[[10]],
  Anatomic = metadata_COG_FN_reorder[[8]],
  Riskgroup = metadata_COG_FN_reorder[[9]],
  Histology = metadata_COG_FN_reorder[[5]],
  FusionStatus = metadata_COG_FN_reorder[[7]],
  col = list(
    Histology  = histology_col,
    FusionStatus = FusionStatus_col,
    Cohort = Cohort_col,
    Sex = Sex_col,
    EnrollmentAge = EnrollmentAge_col,
    Anatomic = Anatomic_col,
    Riskgroup = Riskgroup_col
  ),
  show_legend = c(TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Histology  = list(title = "Histology"),
    FusionStatus = list(title = "Fusion Status"),
    Cohort = list(title= "Cohort"),
    Sex = list(title="Sex"),
    EnrollmentAge = list(title="Age at enrollment"),
    Anatomic = list(title="Anatomic Group"),
    Riskgroup = list(title="Risk Group")
  )
)

hm_COG_FN= Heatmap(variantDataE_COG_FN_memo, 
                   name = "mat", 
                   col = color_list, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   rect_gp = gpar(col= "#e6e6e6"),
                   #rect_gp = gpar(col= NA),
                   show_column_names = FALSE,
                   top_annotation = top_ha_COG_FP,
                   row_names_gp =  gpar(fontsize = 10),
                   row_names_side = c("left"),
                   show_heatmap_legend = FALSE,
                   width = unit(390, "mm"),
)

## UK & FP
metadata_UK_FP     <- metadata %>% 
  tibble::rownames_to_column(var = 'SampleID') %>% 
  dplyr::filter( Cohort == "UK" & `Grouping FP or FN` == "FP") %>% 
  tibble::column_to_rownames(var='SampleID');
variantDataE_UK_FP_memo <- variantDataE[rownames(metadata_UK_FP), ] %>% t() %>% memoSort()
metadata_UK_FP_reorder     <- metadata_UK_FP[colnames(variantDataE_UK_FP_memo), ]
if(!identical(colnames(variantDataE_UK_FP_memo), rownames(metadata_UK_FP_reorder))){
  print( "Please check the metadata and variant data")
}

top_ha_UK_FP = HeatmapAnnotation(
  Cohort = metadata_UK_FP_reorder[[1]],
  Sex = metadata_UK_FP_reorder[[3]],
  EnrollmentAge = metadata_UK_FP_reorder[[10]],
  Anatomic = metadata_UK_FP_reorder[[8]],
  Riskgroup = metadata_UK_FP_reorder[[9]],
  Histology = metadata_UK_FP_reorder[[5]],
  FusionStatus = metadata_UK_FP_reorder[[7]],
  col = list(
    Histology  = histology_col,
    FusionStatus = FusionStatus_col,
    Cohort = Cohort_col,
    Sex = Sex_col,
    EnrollmentAge = EnrollmentAge_col,
    Anatomic = Anatomic_col,
    Riskgroup = Riskgroup_col
  ),
  show_legend = c(TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Histology  = list(title = "Histology"),
    FusionStatus = list(title = "Fusion Status"),
    Cohort = list(title= "Cohort"),
    Sex = list(title="Sex"),
    EnrollmentAge = list(title="Age at enrollment"),
    Anatomic = list(title="Anatomic Group"),
    Riskgroup = list(title="Risk Group")
  )
)

hm_UK_FP= Heatmap(variantDataE_UK_FP_memo, 
                   name = "mat", 
                   col = color_list, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   rect_gp = gpar(col= "#e6e6e6"),
                   #rect_gp = gpar(col= NA),
                   show_column_names = FALSE,
                   top_annotation = top_ha_UK_FP,
                   row_names_gp =  gpar(fontsize = 10),
                   row_names_side = c("left"),
                   show_heatmap_legend = FALSE,
                   width = unit(390, "mm"),
)


## UK & FN
metadata_UK_FN     <- metadata %>% 
  tibble::rownames_to_column(var = 'SampleID') %>% 
  dplyr::filter( Cohort == "UK" & `Grouping FP or FN` == "FN") %>% 
  tibble::column_to_rownames(var='SampleID');
variantDataE_UK_FN_memo <- variantDataE[rownames(metadata_UK_FN), ] %>% t() %>% memoSort()
metadata_UK_FN_reorder  <- metadata_UK_FN[colnames(variantDataE_UK_FN_memo), ]
if(!identical(colnames(variantDataE_UK_FN_memo), rownames(metadata_UK_FN_reorder))){
  print( "Please check the metadata and variant data")
}

top_ha_UK_FN = HeatmapAnnotation(
  Cohort = metadata_UK_FN_reorder[[1]],
  Sex = metadata_UK_FN_reorder[[3]],
  EnrollmentAge = metadata_UK_FN_reorder[[10]],
  Anatomic = metadata_UK_FN_reorder[[8]],
  Riskgroup = metadata_UK_FN_reorder[[9]],
  Histology = metadata_UK_FN_reorder[[5]],
  FusionStatus = metadata_UK_FN_reorder[[7]],
  col = list(
    Histology  = histology_col,
    FusionStatus = FusionStatus_col,
    Cohort = Cohort_col,
    Sex = Sex_col,
    EnrollmentAge = EnrollmentAge_col,
    Anatomic = Anatomic_col,
    Riskgroup = Riskgroup_col
  ),
  show_legend = c(TRUE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Histology  = list(title = "Histology"),
    FusionStatus = list(title = "Fusion Status"),
    Cohort = list(title= "Cohort"),
    Sex = list(title="Sex"),
    EnrollmentAge = list(title="Age at enrollment"),
    Anatomic = list(title="Anatomic Group"),
    Riskgroup = list(title="Risk Group")
  )
)

hm_UK_FN= Heatmap(variantDataE_UK_FN_memo, 
                   name = "mat", 
                   col = color_list, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   rect_gp = gpar(col= "#e6e6e6"),
                   #rect_gp = gpar(col= NA),
                   show_column_names = FALSE,
                   top_annotation = top_ha_UK_FN,
                   row_names_gp =  gpar(fontsize = 10),
                   row_names_side = c("left"),
                   show_heatmap_legend = FALSE,
                   width = unit(390, "mm"),
)

### Concatenate heatmaps
pdf("plot_1.hm_UK_FN.pdf", height = 8, width = 19)

## Displayinh Heatmap
#hm_COG_FP
#hm_COG_FN
#hm_UK_FP
hm_UK_FN

### Annotation 
annotation_titles = c(Cohort = "Cohort",
                      Sex = "Sex",
                      EnrollmentAge = "Age",
                      Histology = "Histology",
                      FusionStatus = "Fusion Status",
                      Anatomic = "Anatomic Group",
                      Riskgroup = "Risk Group"
)

for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

decorate_annotation("Cohort", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Sex", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("EnrollmentAge", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Histology", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("FusionStatus", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Anatomic", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Riskgroup", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})

dev.off()

########## All together ###############
ht_list = hm_COG_FP + hm_COG_FN
pdf("plot_al_Four.1c.pdf", height = 8, width = 30)
draw(ht_list)
### Annotation 
annotation_titles = c(Cohort = "Cohort",
                      Sex = "Sex",
                      EnrollmentAge = "Age",
                      Histology = "Histology",
                      FusionStatus = "Fusion Status",
                      Anatomic = "Anatomic Group",
                      Riskgroup = "Risk Group"
)

for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

decorate_annotation("Cohort", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Sex", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("EnrollmentAge", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Histology", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("FusionStatus", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Anatomic", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("Riskgroup", {
  grid.lines(unit(c(-35, 0), "mm"), unit(c(1, 1), "npc"))
})

dev.off()



























