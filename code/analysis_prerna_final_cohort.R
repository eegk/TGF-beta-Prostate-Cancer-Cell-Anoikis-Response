#############################################################################################################################################
### libraries and functions
#############################################################################################################################################
#############################################################################################################################################
### Load packages for analysis
#############################################################################################################################################
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "annotables","Rsamtools","GenomicFeatures","ggarrange","pheatmap","rms",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","org.Hs.eg.db","treemapify",'ggpubr',
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "colorspace","Vennerable","enrichR","cowplot","data.table",
        "ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","SCENIC","AUCell","RcisTarget","plyr",
        "tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")
lapply(libs, require, character.only = TRUE) ; rm(libs)
#if (!requireNamespace("BiocManager", quietly = TRUE)) 
###
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "ggarrange","pheatmap","rms", "ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace","reshape2",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges","ComplexHeatmap",
        "scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "colorspace","Vennerable","enrichR","cowplot","tidyverse","ggalluvial")
lapply(libs, require, character.only = TRUE) ; rm(libs)
#if (!requireNamespace("BiocManager", quietly = TRUE)) ; install.packages("BiocManager") ; BiocManager::install("survminer")
#############################################################################################################################################
### functions
#############################################################################################################################################
my_filter_list_by_row <- function(x){
  ii <- NULL
  for (i in 1:length(x)) { ii[[i]]<- nrow(x[[i]]) > 0}
  ii <- unlist(ii)
  return(x[ii]) }
#############################################################################################################################################
run_enrichr_on_list <- function(y) {
  
  enrichr_results <- NULL
  generic_var_x <- NULL
  dbs <- listEnrichrDbs()
  dbs <- dbs$libraryName
  
  for (i in 1:length(y)) {
    my_genes <- y[[i]]  
    generic_var_x <- enrichr( genes=my_genes, databases = dbs )
    generic_var_x <- my_filter_list_by_row(generic_var_x)
    
    for (x in 1:length(generic_var_x)) { generic_var_x[[x]]$database <- names(generic_var_x)[x] }
    enrichr_results[[i]] <- do.call(rbind, generic_var_x)
    enrichr_results[[i]]$Comparison <- names(y)[i]
  }
  return(enrichr_results) }
###
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#############################################################################################################################################

#############################################################################################################################################
### get log (in minerva)
#############################################################################################################################################
#temp = list.files(pattern="^Log.final.out$",recursive = TRUE,full.name = TRUE)
#temp <- temp[-grep('_STARpass1',temp)]
#my_log_files = lapply(temp, read.table, sep="\t",header=TRUE,fill=TRUE) 
#names(my_log_files) <- temp
#save.image("logs.RData")
#############################################################################################################################################
### get reads from STAR
#############################################################################################################################################
#temp = list.files(pattern="^ReadsPerGene.out.tab$",recursive = TRUE,full.name = TRUE)
#my_count_files = lapply(temp, read.table, sep="\t",header=TRUE,fill=TRUE) 
#names(my_count_files) <- temp
#save.image("count_files.RData")
#############################################################################################################################################

#############################################################################################################################################
### SET DIRECTORY
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_dogra/projects_kyprianou/prerna_natasha/")
#############################################################################################################################################

#############################################################################################################################################
### stats and data prep
#############################################################################################################################################
load(file = 'counts_RData_2024/count_files.RData')
load(file = 'counts_RData_2024/logs.RData')
#############################################################################################################################################

#############################################################################################################################################
### LOGs
#############################################################################################################################################
names(my_log_files) <- gsub('./align_star/','',names(my_log_files))
names(my_log_files) <- gsub('.fastq.trimmed_fastq/Log.final.out','',names(my_log_files))
names(my_log_files) <- gsub('\\.\\/','',names(my_log_files))
names(my_log_files) <- make.names(names(my_log_files))
names(my_log_files) <- gsub('\\.','_',names(my_log_files))

### add labels
for (i in 1:length(my_log_files)) { my_log_files[[i]]$sample <- names(my_log_files)[i] ; colnames(my_log_files[[i]]) <- c("stat","value","sample")}
### to character
for (i in 1:length(my_log_files)) { my_log_files[[i]] <- my_log_files[[i]] %>% mutate_at(c(1:3), as.character) }
###
for (i in 1:length(my_log_files)) { names(my_log_files[[i]]) <- names(my_log_files[[i]]) <- paste(names(my_log_files[[i]]), names(my_log_files)[i],sep="___") }

### extract
alignment_stats <- do.call(cbind,my_log_files)
alignment_stats <- cbind(alignment_stats[,1],alignment_stats[,grep('value',colnames(alignment_stats))])
alignment_stats <- as.data.frame(alignment_stats)
### sort.process.filter.rename
rownames(alignment_stats) <- alignment_stats$`alignment_stats[, 1]`
alignment_stats <- alignment_stats[,-1]
###
rownames(alignment_stats) <- gsub("%","Percentage",rownames(alignment_stats))
rownames(alignment_stats) <- make.names(rownames(alignment_stats))
rownames(alignment_stats) <- gsub("  ","",rownames(alignment_stats)) # nrepeat
rownames(alignment_stats) <- gsub("[[:punct:]]","_",rownames(alignment_stats))
rownames(alignment_stats) <- gsub("^ ","",rownames(alignment_stats))
rownames(alignment_stats) <- gsub(" $","",rownames(alignment_stats))
rownames(alignment_stats) <- gsub("__","_",rownames(alignment_stats))
rownames(alignment_stats) <- gsub("__","_",rownames(alignment_stats))
alignment_stats$Parameter <- rownames(alignment_stats)
###
alignment_stats <- alignment_stats[grep('percent',rownames(alignment_stats),ignore.case = T),]
alignment_stats$Parameter <- NULL
colnames(alignment_stats) <- tidyr::separate(data.frame(colnames(alignment_stats)),1, sep="___", c("a","b",'c','d'))$b
rownames(alignment_stats) <- gsub("X_","",rownames(alignment_stats))
rownames(alignment_stats) <- gsub("_$","",rownames(alignment_stats))
###
for(i in 1:ncol(alignment_stats)){ alignment_stats[,i] <- gsub('%','',alignment_stats[,i]) }
alignment_stats <- type_convert(alignment_stats)
alignment_stats <- alignment_stats[,-grep('Undetermined',colnames(alignment_stats),ignore.case = T)]
rownames(alignment_stats) <- c("Unique","Mismatch","Multiple_loci",
                               "Many_loci","Unmapped_mismatches","Unmapped_too_short",
                               "Unmapped_other","Chimeric_reads")
temp_params <- melt(as.matrix(alignment_stats[-c(4,8),]))
###
pdf(file="figures_TD006643/alignment_statistics.pdf",width = 5,height = 5)
ggplot(data=temp_params) + theme_classic() +
  aes(y=Var2,x=value,color=Var1,fill=Var1) + 
  geom_bar(stat="identity") +
  geom_text(data=temp_params[temp_params$Var1 %in% "Unique",],aes(label=value), vjust=0.5, color="black") +
  theme(legend.position = "bottom") + 
  labs(x="Percent",y="Samples",title = "Alignment Statistics")
dev.off()
### clean
results_analysis <- list()
results_analysis$alignment_stats <- alignment_stats
rm(temp_params,my_log_files,alignment_stats,i,temp)
#############################################################################################################################################

#############################################################################################################################################
### Counts
#############################################################################################################################################

names(my_count_files) <- gsub('.fastq.trimmed_fastq/ReadsPerGene.out.tab','',names(my_count_files))
names(my_count_files) <- gsub('./','',names(my_count_files))

third_columns <- lapply(my_count_files, function(df) df[,3])
my_counts <- do.call(cbind,third_columns)
rownames(my_counts) <- my_count_files[[1]][,1]
results_analysis$my_counts <- my_counts

results_analysis$my_counts <- results_analysis$my_counts[,-13]

library(org.Hs.eg.db)
results_analysis$my_genes <- data.frame( ensembl=rownames(results_analysis$my_counts),
                                         symbol=as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$my_counts)), column="SYMBOL", keytype="ENSEMBL", multiVals="first")) )
rm(third_columns,my_count_files,my_counts)

results_analysis$alignment_stats <- t(results_analysis$alignment_stats)
results_analysis$alignment_stats <- as.data.frame(results_analysis$alignment_stats)
results_analysis$alignment_stats$sampleID <- colnames(results_analysis$my_counts)
colnames(results_analysis$my_counts) <- rownames(results_analysis$alignment_stats)

results_analysis$alignment_stats$cell_line <- "TRII"
results_analysis$alignment_stats$cell_line[grep('VCaP',rownames(results_analysis$alignment_stats))] <- 'VCaP'
results_analysis$alignment_stats$treatment <- c('DZ50','DZ50','DZ50',
                                                'TGFB','TGFB',
                                                'TGFB_DZ50','TGFB_DZ50','TGFB_DZ50',
                                                'TGFB',
                                                'Untreated','Untreated','Untreated',
                                                'Untreated','Untreated','Untreated',
                                                'DZ50','DZ50','DZ50',
                                                'TGFB','TGFB',
                                                'TGFB_DZ50','TGFB_DZ50','TGFB_DZ50',
                                                'TGFB')
###
results_analysis$alignment_stats$individualID <- paste(results_analysis$alignment_stats$cell_line,results_analysis$alignment_stats$treatment)
results_analysis$alignment_stats$sampleID <- rownames(results_analysis$alignment_stats)
#############################################################################################################################################

#############################################################################################################################################
### checkpoint
#############################################################################################################################################
save.image(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################

#############################################################################################################################################
### REVIEW
#############################################################################################################################################
mydata <- prcomp(t( cpm(results_analysis$my_counts[rowSums(results_analysis$my_counts)>0,]) ), scale=TRUE,center = TRUE) 
pca_mts<-as.data.frame(mydata$x)
pca_mts$sampleID <- rownames(pca_mts)
pca1 <- format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2)
pca3 <- format(round(get_eigenvalue(mydata)$variance.percent[3], 2), nsmall = 2)
pca_mts <- merge(pca_mts,results_analysis$alignment_stats,by='sampleID')

pdf(file="figures_TD006643/pca_rnaseq_screeplot_samples_TD006643_rawdata.pdf",width = 3,height = 2)
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_classic() + ylim(0,25)
dev.off()

pdf(file="figures_TD006643/pca_rnaseq_PC1_2_samples_TD006643_rawdata.pdf",width = 7.5,height = 7)
ggplot(pca_mts, aes(PC2, PC1, color = paste(cell_line,treatment) )) + 
    geom_point(show.legend = TRUE) + theme_classic() + 
    geom_text_repel(aes(label=individualID),max.overlaps = 100,show.legend = F)+
    labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Type / Origin") 
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### normalization
#############################################################################################################################################
###
colSums(results_analysis$my_counts)
grep("TRII-TGFB",colnames(results_analysis$my_counts),value = F)
ix <- colnames(results_analysis$my_counts)[c(4,5,9)]
### investigate
table(rowSums(cpm(results_analysis$my_counts)>1)>=1)
table(rowSums(cpm(results_analysis$my_counts[,ix])>1)>=1)
table(rowSums(cpm(results_analysis$my_counts)>=0.2)>=1)
table(rowSums(cpm(results_analysis$my_counts)>=0.2)>=3)
#table(rowSums(cpm(results_analysis$my_counts)>5)>=3)
#table(rowSums(cpm(results_analysis$my_counts)>5)>=5)
#table(rowSums(cpm(results_analysis$my_counts)>9)>=5)
###
keep.counts <- which(rowSums(cpm(results_analysis$my_counts)>=0.2)>=3)
###
deglist_obj <- DGEList(counts = results_analysis$my_counts[keep.counts,],
                       genes = results_analysis$my_genes[keep.counts,])
###
deglist_obj <- calcNormFactors(deglist_obj, method = "TMM") 
deglist_obj <- estimateCommonDisp(deglist_obj, verbose=TRUE)

### dream

### formula
form = ~ cell_line + treatment + (1|individualID)
###
results_analysis$voom_TD006643 <- voomWithDreamWeights(counts = deglist_obj,
                                                       formula = form,
                                                       data = results_analysis$alignment_stats,
                                                       plot = T,save.plot = T) 
results_analysis$voom_TD006643$pseudocounts <- deglist_obj$pseudo.counts
results_analysis$voom_TD006643$samples <- deglist_obj$samples
dev.off()
### 
rm(deglist_obj,temp,keep.counts,ix)
#############################################################################################################################################

#############################################################################################################################################
### color pallete
#############################################################################################################################################
#hcl_ramp(4,"paleturquoise","darkblue") rm(hcl_ramp)
ramps_blue <- c(colorRamp2(breaks = c(1,4),colors = c('paleturquoise',"blue4"))(1:4),
                colorRamp2(breaks = c(1,4),colors = c("plum1","red4"))(1:4))

results_analysis$color_pal <- list(
  `TRII DZ50` = ramps_blue[7],     
  `TRII TGFB` = ramps_blue[6],
  `TRII TGFB_DZ50` = ramps_blue[8],
  `TRII Untreated`= ramps_blue[5],
  `VCaP DZ50` = ramps_blue[3],    
  `VCaP TGFB` = ramps_blue[2],
  `VCaP TGFB_DZ50` = ramps_blue[4],
  `VCaP Untreated` = ramps_blue[1])

rm(ramps_blue)
#############################################################################################################################################

### PCA + heatmap and other figures
#############################################################################################################################################
mydata <- prcomp(t(results_analysis$voom_TD006643$E[,]), scale=TRUE,center = TRUE) 
pca_mts<-as.data.frame(mydata$x)
pca_mts$sampleID <- rownames(pca_mts)
pca1 <- format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2)
pca3 <- format(round(get_eigenvalue(mydata)$variance.percent[3], 2), nsmall = 2)
pca_mts <- merge(pca_mts,results_analysis$alignment_stats,by='sampleID')

pdf(file="figures_TD006643/pca_rnaseq_screeplot_samples_TD006643_normalized.pdf",width = 3,height = 2)
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_classic() + ylim(0,55)
dev.off()


pdf(file="figures_TD006643/pca_rnaseq_PCs1to3_samples_TD006643_normalized.pdf",width = 4,height = 3.5)
plot_grid(
  ggplot(pca_mts, aes(PC2, PC1, color = individualID )) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  cowplot::get_legend(ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point() + theme_classic() + 
                        scale_color_manual(values=results_analysis$color_pal) +guides(guide_legend(ncol=2)) + 
                        theme(legend.position = 'right') +
                        labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment")) 
  ,
  ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  ggplot(pca_mts, aes(PC3, PC2, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC2:", pca2,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  labels = c("","","","","","","","","") 
)
dev.off()

#############################################################################################################################################
### only VCaP

mydata <- prcomp(t(results_analysis$voom_TD006643$E[,grep('VCaP',colnames(results_analysis$voom_TD006643$E))]), scale=TRUE,center = TRUE) 
pca_mts<-as.data.frame(mydata$x)
pca_mts$sampleID <- rownames(pca_mts)
pca1 <- format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2)
pca3 <- format(round(get_eigenvalue(mydata)$variance.percent[3], 2), nsmall = 2)
pca_mts <- merge(pca_mts,results_analysis$alignment_stats,by='sampleID')

pdf(file="figures_TD006643/pca_rnaseq_screeplot_samples_TD006643_normalized_VCaP_Only.pdf",width = 3,height = 2)
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_classic() + ylim(0,35)
dev.off()


pdf(file="figures_TD006643/pca_rnaseq_PCs1to3_samples_TD006643_normalized_VCaP_Only.pdf",width = 4,height = 3.5)
plot_grid(
  ggplot(pca_mts, aes(PC2, PC1, color = individualID )) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  cowplot::get_legend(ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point() + theme_classic() + 
                        scale_color_manual(values=results_analysis$color_pal) +guides(guide_legend(ncol=2)) + 
                        theme(legend.position = 'right') +
                        labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment")) 
  ,
  ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  ggplot(pca_mts, aes(PC3, PC2, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC2:", pca2,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  labels = c("","","","","","","","","") 
)
dev.off()

#############################################################################################################################################
### only TRII

mydata <- prcomp(t(results_analysis$voom_TD006643$E[,grep('TRII',colnames(results_analysis$voom_TD006643$E))]), scale=TRUE,center = TRUE) 
pca_mts<-as.data.frame(mydata$x)
pca_mts$sampleID <- rownames(pca_mts)
pca1 <- format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2)
pca3 <- format(round(get_eigenvalue(mydata)$variance.percent[3], 2), nsmall = 2)
pca_mts <- merge(pca_mts,results_analysis$alignment_stats,by='sampleID')

pdf(file="figures_TD006643/pca_rnaseq_screeplot_samples_TD006643_normalized_TRII_Only.pdf",width = 3,height = 2)
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_classic() + ylim(0,50)
dev.off()

pdf(file="figures_TD006643/pca_rnaseq_PCs1to3_samples_TD006643_normalized_TRII_Only.pdf",width = 4,height = 3.5)
plot_grid(
  ggplot(pca_mts, aes(PC2, PC1, color = individualID )) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  cowplot::get_legend(ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point() + theme_classic() + 
                        scale_color_manual(values=results_analysis$color_pal) +guides(guide_legend(ncol=2)) + 
                        theme(legend.position = 'right') +
                        labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment")) 
  ,
  ggplot(pca_mts, aes(PC3, PC1, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  ggplot(pca_mts, aes(PC3, PC2, color = individualID)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=results_analysis$color_pal) + 
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC2:", pca2,"%", sep=""), color="Cell line &\nTreatment") 
  ,
  labels = c("","","","","","","","","") 
)
dev.off()

rm(pca_mts,mydata,curWarnings)
#############################################################################################################################################


#############################################################################################################################################
### differential expression
#############################################################################################################################################
identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)

form <- ~ 0 + individualID # + (1|individualID) 

L =  makeContrastsDream( form, results_analysis$alignment_stats, 
                         contrasts = c('`individualIDVCaP Untreated` - `individualIDTRII Untreated`',
                                       '`individualIDVCaP DZ50` - `individualIDTRII DZ50`',
                                       '`individualIDVCaP TGFB` - `individualIDTRII TGFB`',
                                       '`individualIDVCaP TGFB_DZ50` - `individualIDTRII TGFB_DZ50`',
                                       
                                       '`individualIDVCaP Untreated` - `individualIDVCaP DZ50`',
                                       '`individualIDVCaP Untreated` - `individualIDVCaP TGFB`',
                                       '`individualIDVCaP Untreated` - `individualIDVCaP TGFB_DZ50`',
                                       '`individualIDVCaP DZ50` - `individualIDVCaP TGFB`',
                                       '`individualIDVCaP DZ50` - `individualIDVCaP TGFB_DZ50`',
                                       '`individualIDVCaP TGFB` - `individualIDVCaP TGFB_DZ50`',
                                       
                                       '`individualIDTRII Untreated` - `individualIDTRII DZ50`',
                                       '`individualIDTRII Untreated` - `individualIDTRII TGFB`',
                                       '`individualIDTRII Untreated` - `individualIDTRII TGFB_DZ50`',
                                       '`individualIDTRII DZ50` - `individualIDTRII TGFB`',
                                       '`individualIDTRII DZ50` - `individualIDTRII TGFB_DZ50`',
                                       '`individualIDTRII TGFB` - `individualIDTRII TGFB_DZ50`'))

fitmm = dream( results_analysis$voom_TD006643, form, results_analysis$alignment_stats, L)
fitmm = eBayes(fitmm)

fitmAV = dream( results_analysis$voom_TD006643, form, results_analysis$alignment_stats)
marginal_means_model <- fitmAV$coefficients

lmfreq_results <- list() ; for ( i in 1:ncol(L) ) { lmfreq_results[[i]] <- topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(fitmm$coefficients)[i]
lmfreq_results[[i]]$Marker <- rownames( topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") ) }

differential_expression_pr <- do.call(rbind,lmfreq_results)
rownames(differential_expression_pr) <- NULL
differential_expression_pr$nLogFDR <- -log10(differential_expression_pr$adj.P.Val)
differential_expression_pr$Comparison <- gsub('individualID','',differential_expression_pr$Comparison)

results_analysis$de_v1 <- differential_expression_pr
results_analysis$fitmm_v1 <- fitmm
results_analysis$fitmAV_v1 <- fitmAV
results_analysis$fitmAV_v1 <- marginal_means_model

rm(L,fitmm,lmfreq,form,i,simple_model_test,differential_expression_pr,lmfreq_results,x,variance_olink,ii,ix,iix,my_vars,count,differential_expression_pr,fitmAV,marginal_means_model)
#############################################################################################################################################


#############################################################################################################################################
### summary barplot
#############################################################################################################################################
table(results_analysis$de_v1$Comparison[results_analysis$de_v1$adj.P.Val<0.05])
#############################################################################################################################################
stacked_bar_plot <- results_analysis$de_v1
stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05
stacked_bar_plot <- data.frame(table(stacked_bar_plot$Comparison, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
colnames(stacked_bar_plot) <- c('Comparison', 'logFC_dir','sig', 'count')
regional_dmr_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down'] <- 0 - regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down']
regional_dmr_summary$Comparison <- as.character(regional_dmr_summary$Comparison) 
regional_dmr_summary$Comparison <- factor(regional_dmr_summary$Comparison,levels = unique(results_analysis$de_v1$Comparison))
#############################################################################################################################################
pdf(file="figures_TD006643/stack_deg_signals_v1.pdf",width = 5,height = 3)
ggplot(regional_dmr_summary) + 
  aes(x=count, y=Comparison, fill=logFC_dir, label=abs(count)) + 
  geom_bar(stat='identity', width=0.5) +
  theme_classic() + theme(axis.text.y = element_text(face='bold'),
                          axis.text.x = element_text(face='bold'),
                          axis.title = element_text(face = "bold"),
                          legend.title= element_text(face='bold'),
                          strip.text=element_text(face='bold'),
                          plot.title=element_text(face='bold')) + 
  scale_fill_manual(values=c('firebrick', 'steelblue'), name='LogFC Direction') +
  labs(x='no. DEGs (FDR < 0.05)', y='Contrast') + theme(legend.position="bottom")
dev.off()

rm(regional_dmr_summary,stacked_bar_plot,pca1,pca2,pca3)
#############################################################################################################################################


#############################################################################################################################################
### Hallmark pathways
#############################################################################################################################################
library(fgsea)
library(GSVA)
library(ggbeeswarm)

### MsigDB hallmark genes
my_msigdb_oncogenic <- gmtPathways("/Users/gonzae34/Documents/msig_pathways/h.all.v6.2.symbols.gmt")
my_msigdb_immune <- gmtPathways("/Users/gonzae34/Documents/msig_pathways/c7.all.v7.1.symbols.gmt")

ExprDat <- results_analysis$voom_TD006643$E 
rownames(ExprDat) <- results_analysis$voom_TD006643$genes$symbol
ExprDat <- ExprDat[!is.na(rownames(ExprDat)),]

gsva.de.oncogenic <- gsva(ExprDat, my_msigdb_oncogenic, verbose=FALSE)
gsva.de.immune <- gsva(ExprDat, my_msigdb_immune, verbose=FALSE)

pdf(file='figures_TD006643/hallmarker_signatures_oncogenic_gsva.pdf',width = 8.5, height = 9)
pheatmap::pheatmap(gsva.de.oncogenic,scale = 'none',angle_col = 90, main='OncoGenic')
dev.off()

x <- melt(gsva.de.oncogenic)
colnames(x)[2] <- 'sampleID'
x <- merge(x,results_analysis$alignment_stats,by='sampleID')

x$Var1 <- gsub("HALLMARK_","",x$Var1)

my_comp <- list( c("TRII Untreated","TRII DZ50"),
                 c("TRII Untreated","TRII TGFB"),
                 c("TRII Untreated","TRII TGFB_DZ50"),
                 
                 c("VCaP Untreated","VCaP DZ50"),
                 c("VCaP Untreated","VCaP TGFB"),
                 c("VCaP Untreated","VCaP TGFB_DZ50"),
                 
                 c("TRII Untreated","VCaP Untreated"))

pdf(file='figures_TD006643/hallmarker_signatures_oncogenic_gsva_boxplots_t_test.pdf',width = 4.5, height = 4.5)
for( i in 1:length(unique(x$Var1) )){
  y <- x[x$Var1 %in% unique(x$Var1)[i], ]
  print(
  ggplot(data=y)+aes(x=individualID,y=value) + theme_classic() + 
    geom_boxplot(aes(fill=individualID),show.legend = F)+ scale_fill_manual(values = results_analysis$color_pal)+
    geom_quasirandom()+
    stat_compare_means(comparisons = my_comp,method = 't.test')+
    rotate_x_text(angle = 45) +
    labs(x='Cell line * Treatment',y='GSVA Score',title = unique(x$Var1)[i]))
}
dev.off()

pdf(file='figures_TD006643/hallmarker_signatures_oncogenic_gsva_boxplots_wilcox.pdf',width = 4.5, height = 4.5)
for( i in 1:length(unique(x$Var1) )){
y <- x[x$Var1 %in% unique(x$Var1)[i], ]
print(
ggplot(data=y)+aes(x=individualID,y=value) + theme_classic() + 
  geom_boxplot(aes(fill=individualID),show.legend = F)+ scale_fill_manual(values = results_analysis$color_pal)+
  geom_quasirandom()+
  stat_compare_means(comparisons = my_comp,method = 'wilcox.test')+
  rotate_x_text(angle = 45) +
  labs(x='Cell line * Treatment',y='GSVA Score',title = unique(x$Var1)[i]))
}
dev.off()

pdf(file='figures_TD006643/hallmarker_signatures_oncogenic_gsva_boxplots_paired_t_test.pdf',width = 4.5, height = 4.5)
for( i in 1:length(unique(x$Var1) )){
  y <- x[x$Var1 %in% unique(x$Var1)[i], ]
  print(
    ggplot(data=y)+aes(x=individualID,y=value) + theme_classic() + 
      geom_boxplot(aes(fill=individualID),show.legend = F)+ scale_fill_manual(values = results_analysis$color_pal)+
      geom_quasirandom()+
      stat_compare_means(comparisons = my_comp,method = 't.test',paired = TRUE)+
      rotate_x_text(angle = 45) +
      labs(x='Cell line * Treatment',y='GSVA Score',title = unique(x$Var1)[i]))
}
dev.off()

pdf(file='figures_TD006643/hallmarker_signatures_oncogenic_gsva_boxplots_paired_wilcox.pdf',width = 4.5, height = 4.5)
for( i in 1:length(unique(x$Var1) )){
  y <- x[x$Var1 %in% unique(x$Var1)[i], ]
  print(
    ggplot(data=y)+aes(x=individualID,y=value) + theme_classic() + 
      geom_boxplot(aes(fill=individualID),show.legend = F)+ scale_fill_manual(values = results_analysis$color_pal)+
      geom_quasirandom()+
      stat_compare_means(comparisons = my_comp,method = 'wilcox.test',paired = TRUE)+
      rotate_x_text(angle = 45) +
      labs(x='Cell line * Treatment',y='GSVA Score',title = unique(x$Var1)[i]))
}
dev.off()

results_analysis$hallmark_signature_list <- list(immune=my_msigdb_immune,oncogenic=my_msigdb_oncogenic)
results_analysis$gsva.de.immune <- gsva.de.immune
results_analysis$gsva.de.oncogenic <- gsva.de.oncogenic
#############################################################################################################################################

#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#rm(list=setdiff(ls(), "results_analysis"))
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################

#############################################################################################################################################
### biotype table
#############################################################################################################################################
#library(ensembldb)
#library(AnnotationHub)
gtf <- rtracklayer::import('/Users/gonzae34/Documents/genomic_references/Homo_sapiens.GRCh38.102.gtf.gz')
gtf_df=as.data.frame(gtf) ; rm(gtf)
gtf_df <- gtf_df[,c('gene_id','gene_name','gene_biotype',"seqnames")] ### "start","end","width","strand"
gtf_df <- unique(gtf_df)
results_analysis$gene_gtf <- gtf_df
colnames(results_analysis$gene_gtf)[1] <- 'ensembl'
rm(gtf_df)
###
results_analysis$gene_gtf <- results_analysis$gene_gtf[results_analysis$gene_gtf$ensembl %in% results_analysis$de_v1$ensembl,]
###
table(results_analysis$de_v1$ensembl %in% results_analysis$gene_gtf$ensembl)
results_analysis$de_v1[results_analysis$de_v1$ensembl %in% results_analysis$gene_gtf$ensembl,]
#############################################################################################################################################


#############################################################################################################################################
### Volcano Plots
#############################################################################################################################################
results_analysis$de_v1$symbol_V2 <- results_analysis$de_v1$symbol
results_analysis$de_v1$symbol_V2[is.na(results_analysis$de_v1$symbol)] <- results_analysis$de_v1$ensembl[is.na(results_analysis$de_v1$symbol)]

pdf(file='figures_TD006643/DE_volcano_plots_v1.pdf',width = 4.5, height = 4.5)

for( i in 1:length(unique(results_analysis$de_v1$Comparison))){
  
  x <- results_analysis$de_v1[results_analysis$de_v1$Comparison %in% unique(results_analysis$de_v1$Comparison)[i],]
  x$color <- x$logFC>0 
  x$color[x$adj.P.Val>0.05] <- 'NS'
  print(
  ggplot(data=x) + aes(x=logFC,y=nLogFDR) + theme_classic() + geom_point(aes(color=color),show.legend = F) + xlim(-(max(x$logFC)+3),(max(x$logFC)+3))+
    scale_color_manual(values = c("#F8766D","grey50","#619CFF"))+
    geom_text_repel(data=x[c(which(x$logFC>0)[1:6],which(x$logFC<0)[1:6]),],aes(label=symbol_V2),max.overlaps = 100,show.legend = F,force = 37)+
    geom_vline(xintercept = 0,linetype=2,color='deeppink') + geom_hline(yintercept = 1.3,linetype=2,color='deeppink') +
    labs(x='Log2FC',y='-Log10(FDR)',title = unique(results_analysis$de_v1$Comparison)[i])
  )
  
  write.csv(file=paste('tables/results_DE',gsub("\\.\\.","\\.",gsub("^X","",make.names(unique(results_analysis$de_v1$Comparison)[i]))),'.csv',sep=""),x,row.names=FALSE)
  
}
dev.off()


for( i in 1:length(unique(results_analysis$de_v1$Comparison))) {
  
    x <- results_analysis$de_v1[results_analysis$de_v1$Comparison %in% unique(results_analysis$de_v1$Comparison)[i],]
    x$nLogFDR<-NULL
    x$Marker<-NULL
    
    if( table(x$adj.P.Val<0.05 & abs(x$logFC)>1)[2] < 100 ) {
      print(table(x$adj.P.Val<0.1 & abs(x$logFC)>0.5))  
      x <- x[which(x$adj.P.Val<0.1 & abs(x$logFC)>1),] 
      write.csv(file=paste('Clean_Tables_For_Prerna/results_simplified_DE',gsub("\\.\\.","\\.",gsub("^X","",make.names(unique(results_analysis$de_v1$Comparison)[i]))),'csv',sep=""),x,row.names=FALSE)  
      
    } else {
      x <- x[which(x$adj.P.Val<0.05 & abs(x$logFC)>1),] 
      write.csv(file=paste('Clean_Tables_For_Prerna/results_simplified_DE',gsub("\\.\\.","\\.",gsub("^X","",make.names(unique(results_analysis$de_v1$Comparison)[i]))),'csv',sep=""),x,row.names=FALSE)  
    }

  }


#############################################################################################################################################

#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#rm(list=setdiff(ls(), "results_analysis"))
#save.image(file = 'RData/analysis_pr_TD006643.RData')
load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################

#############################################################################################################################################
### Heatmaps Plots
#############################################################################################################################################
identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)
###
my_msigdb_oncogenic <- gmtPathways("/Users/gonzae34/Documents/msig_pathways/h.all.v6.2.symbols.gmt")
my_msigdb_immune <- gmtPathways("/Users/gonzae34/Documents/msig_pathways/c7.all.v7.1.symbols.gmt")

### what is the treatment effect on cell lines
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `TRII Untreated`" & results_analysis$de_v1$adj.P.Val<0.05)
ib <- which(results_analysis$de_v1$Comparison %in% "`VCaP TGFB_DZ50` - `TRII TGFB_DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
ic <- which(results_analysis$de_v1$Comparison %in% "`VCaP DZ50` - `TRII DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
id <- which(results_analysis$de_v1$Comparison %in% "`VCaP TGFB` - `TRII TGFB`" & results_analysis$de_v1$adj.P.Val<0.05)

x <- list(
  #Untreated=results_analysis$de_v1$ensembl[ia],
  TGFB_DZ50=results_analysis$de_v1$ensembl[ib],
  DZ50=results_analysis$de_v1$ensembl[ic],
  TGFB=results_analysis$de_v1$ensembl[id]
)

y <- list(
  #Untreated=results_analysis$de_v1$ensembl[ia][!results_analysis$de_v1$ensembl[ia] %in% Reduce(intersect,x)],
  TGFB_DZ50=results_analysis$de_v1$ensembl[ib][!results_analysis$de_v1$ensembl[ib] %in% Reduce(intersect,x)],
  DZ50=results_analysis$de_v1$ensembl[ic][!results_analysis$de_v1$ensembl[ic] %in% Reduce(intersect,x)],
  TGFB=results_analysis$de_v1$ensembl[id][!results_analysis$de_v1$ensembl[id] %in% Reduce(intersect,x)]
)

library(Vennerable)

pdf(file='figures_TD006643/venn_signatures_intersect_cap_vs_tri_treatments.pdf',width = 4,height = 4)
plot(Venn(x))
dev.off()
#plot(Venn(y))

#ix <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `TRII Untreated`")
#x <- results_analysis$de_v1[ix,]

#"EPI" "E-cadherin, B-catenin, Cytokeratin" #c("CDH1","CTNNB1","KRTs")
#"MES" "N-cadherin, Vimentin, Fibronectin, Snail, Slug, Twist or Matrix Metalloproteinases-2,3,9" # c("CDH1","KLRG1,"FN1","SNAI1","SNAI2",SNAI3","TWIST1","TWIST2")

results_analysis$voom_TD006643$genes$symbol_u <- make.unique(results_analysis$voom_TD006643$genes$symbol)
rownames(results_analysis$voom_TD006643$E) <- results_analysis$voom_TD006643$genes$symbol_u

ix <- grep('^CDH1$|CTNNB1|^KRT|KLRG1|^FN1|^SNAI|^TWIST|^PCA|^KLK',rownames(results_analysis$voom_TD006643$E),value = T)

identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)

df_ann <- data.frame(#Line=results_analysis$alignment_stats$cell_line,
                     #Treat=results_analysis$alignment_stats$treatment,
                     LiTr=results_analysis$alignment_stats$individualID)
rownames(df_ann) <- results_analysis$alignment_stats$sampleID

color_ann <- list(LiTr=unlist(results_analysis$color_pal))
               #Line = c(TRII='#FFBBFFFF',VCaP='#AFEEEEFF'),
               #Treat = c(TRII='#FFBBFFFF',VCaP='#AFEEEEFF'),)

pdf(file='figures_TD006643/heatmap_known_markers.pdf',width = 8,height = 7)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()

iy <- grep('vcap',df_ann$LiTr,ignore.case = T)

pdf(file='figures_TD006643/heatmap_known_markers_vcap.pdf',width = 6,height = 7)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()

pdf(file='figures_TD006643/heatmap_known_markers_trii.pdf',width = 6,height = 7)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,-iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()

#############################################################################################################################################
### untreated
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `TRII Untreated`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP Untreated","TRII Untreated")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_trii_untreated.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### TGFB_DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP TGFB_DZ50` - `TRII TGFB_DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP TGFB_DZ50","TRII TGFB_DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_trii_TGFB_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP DZ50` - `TRII DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP DZ50","TRII DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_trii_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### TGFB
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP TGFB` - `TRII TGFB`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP TGFB","TRII TGFB")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_trii_TGFB.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### VCaP TGFB
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `VCaP TGFB`" & results_analysis$de_v1$adj.P.Val<0.2)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
#x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP Untreated","VCaP TGFB")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_TGFB.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### VCaP DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `VCaP DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP Untreated","VCaP DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### VCaP TGFB_DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `VCaP TGFB_DZ50`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("VCaP Untreated","VCaP TGFB_DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_vcap_TGFB_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### TRII TGFB
ia <- which(results_analysis$de_v1$Comparison %in% "`TRII Untreated` - `TRII TGFB`" & results_analysis$de_v1$adj.P.Val<0.05)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("TRII Untreated","TRII TGFB")

pdf(file='figures_TD006643/heatmap_top50_markers_TRII_TGFB.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### TRII DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`TRII Untreated` - `TRII DZ50`" & results_analysis$de_v1$adj.P.Val<0.2)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
#x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("TRII Untreated","TRII DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_TRII_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################
### VCaP TGFB_DZ50
ia <- which(results_analysis$de_v1$Comparison %in% "`TRII Untreated` - `TRII TGFB_DZ50`" & results_analysis$de_v1$adj.P.Val<0.1)
x <- results_analysis$de_v1[ia,]
dim(x)
summary(x$logFC)
x <- x[which(abs(x$logFC)>1 & !is.na(x$symbol)),]
#write.csv(file)
ix <- which(results_analysis$voom_TD006643$genes$ensembl %in% x$ensembl[c(which(x$logFC>0)[1:25],which(x$logFC<0)[1:25])])
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% results_analysis$voom_TD006643$genes$symbol_u[ix])

iy <- df_ann$LiTr %in% c("TRII Untreated","TRII TGFB_DZ50")

pdf(file='figures_TD006643/heatmap_top50_markers_TRII_TGFB_DZ50.pdf',width = 5,height = 9)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[ix,iy],
                   scale = 'row',
                   angle_col= 90,
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()
#############################################################################################################################################



#############################################################################################################################################
### Pathways
#############################################################################################################################################

dbs <- listEnrichrDbs()
dbs <- dbs$libraryName

ix <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `VCaP DZ50`" & results_analysis$de_v1$adj.P.Val<0.05 & results_analysis$de_v1$logFC>1)
my_genes <- results_analysis$de_v1$symbol[ix][!is.na(results_analysis$de_v1$symbol[ix])]
down_regulation_VCaP_dz50 <- enrichr( genes=my_genes, databases = dbs )

ix <- which(results_analysis$de_v1$Comparison %in% "`VCaP Untreated` - `VCaP DZ50`" & results_analysis$de_v1$adj.P.Val<0.05 & results_analysis$de_v1$logFC< -1)
my_genes <- results_analysis$de_v1$symbol[ix][!is.na(results_analysis$de_v1$symbol[ix])]
up_regulation_VCaP_dz50 <- enrichr( genes=my_genes, databases = dbs )

head(down_regulation_VCaP_dz50$Reactome_2016)
head(up_regulation_VCaP_dz50$Reactome_2016)

head(down_regulation_VCaP_dz50$GO_Molecular_Function_2023)
head(up_regulation_VCaP_dz50$GO_Molecular_Function_2023)

#############################################################################################################################################
"PARP-1"
#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################
table(results_analysis$de_v1$Comparison)

ix <- grep('^SRC',results_analysis$de_v1$symbol_V2[results_analysis$de_v1$adj.P.Val<0.05],ignore.case = T)

results_analysis$de_v1[results_analysis$de_v1$adj.P.Val<0.05,][ix,]

ix <- grep('^PARP1$',results_analysis$de_v1$symbol_V2[results_analysis$de_v1$adj.P.Val<0.05],ignore.case = T)

results_analysis$de_v1[results_analysis$de_v1$adj.P.Val<0.05,][ix,]

ix <- grep('^CFL1$',results_analysis$de_v1$symbol_V2[results_analysis$de_v1$adj.P.Val<0.05],ignore.case = T)

results_analysis$de_v1[results_analysis$de_v1$adj.P.Val<0.05,][ix,]

#############################################################################################################################################
#############################################################################################################################################
### ARTICLE REVIEW AND FIGURES
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### FIGURE 1
#> TGF-β induces a differential response in prostate cancer cells. 
#Figure 1. (apoptosis) #1. Prerna share the pdfs of her figures, #2. redo figures that apoptosis, heatmaps, volcano plots.
#############################################################################################################################################
### HAETMAP

identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)
#results_analysis$gsva.de.oncogenic
# subset apoptosis
ix <- which(rownames(results_analysis$voom_TD006643$E) %in% my_msigdb_oncogenic$HALLMARK_APOPTOSIS)
tmp_df <- results_analysis$voom_TD006643$E[ix,]

ix <- which(rownames(tmp_df) %in% unique(results_analysis$de_v1$symbol_V2[results_analysis$de_v1$P.Value<0.01]))
tmp_df <- tmp_df[ix,]

ann_col <- data.frame(Treatment=results_analysis$alignment_stats$treatment,
                      Cell=results_analysis$alignment_stats$cell_line)
rownames(ann_col) <- results_analysis$alignment_stats$sampleID

#summary(rowVars(tmp_df))
#ix <- which(rowVars(tmp_df)>1.5)

### FOR TRII
iy <- which(results_analysis$alignment_stats$cell_line %in% 'TRII')

#pheatmap::pheatmap(tmp_df[,iy],annotation_col = ann_col,scale='row',show_colnames = F,
#                   color = colorRampPalette(c("steelblue","white","firebrick"))(20))
#dev.off()
ann_col_v2 <- ann_col[iy,]

z <- t(scale(t(tmp_df[,iy]))) # transpose, scale (z-scoring) / normalization, tranpose
#z <- z[,match(rownames(ann_col_v2),colnames(z))] # reorder z based on ann_col
identical(rownames(ann_col_v2),colnames(z))

#rownames(ann_row) <- ann_row$marker
#ann_row <- ann_row[match(rownames(z),rownames(ann_row)),]
#ann_row$marker<-NULL
#identical(rownames(ann_row),rownames(z))

### ROWS
#row_ha = rowAnnotation(Type = anno_block(gp = gpar(fill = c('grey75','black','grey45')), 
#                                         labels = c('Lipoproteins','Hallmark','Adhesion'),
#                                         labels_gp = gpar(col = "black", fontsize = 10)) )

### results_analysis$color_pal

### COLUMNS
column_ha = HeatmapAnnotation(
  Treatment = anno_block(gp = gpar(fill = c('#BA4A50FF','#E082A4FF','#8B0000FF','#FFBBFFFF')), 
                    labels = c('DZ50','TGFB','TGFB_DZ50','Untreated'),
                    labels_gp = gpar(col = "black", fontsize = 10)) )

### ann_row$Type <- factor(x = ann_row$Type,levels = c('Lipoproteins','Adhesion','Hallmark')) ##

ann_col_v2$Treatment <- factor(ann_col_v2$Treatment,levels = c('DZ50','TGFB','TGFB_DZ50','Untreated'))

pdf(file='figures_TD006643/heatmap_apoptosis_TRII.pdf',width = 6, height = 10)
draw(
  ComplexHeatmap::Heatmap(z, 
                          col= colorRampPalette(c("steelblue", "white","firebrick"))(255),
                          show_row_names = TRUE , show_column_names = FALSE,
                          #column_title_gp = gpar(fontsize = 1, fontface = "plain"), #column_names_gp = gpar(fontsize = 5, fontface = "plain"),
                          column_names_rot = 45,
                          column_split = ann_col_v2$Treatment, #column_km = 7,
                          #row_split =  ann_row, #row_split = 4, #row_km = 5,
                          clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
                          cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                          top_annotation = column_ha, 
                          #left_annotation = row_ha,
                          border = TRUE, use_raster = FALSE, name = "LogFC") ,
  
  heatmap_legend_side = "left", annotation_legend_side = "left" , padding = unit(c(2, 2, 2, 20), "mm")
)
dev.off()


### FOR VCaP
iy <- which(results_analysis$alignment_stats$cell_line %in% 'VCaP')
ann_col_v2 <- ann_col[iy,]

z <- t(scale(t(tmp_df[,iy]))) # transpose, scale (z-scoring) / normalization, tranpose
identical(rownames(ann_col_v2),colnames(z))

### results_analysis$color_pal

### COLUMNS
column_ha = HeatmapAnnotation(
  Treatment = anno_block(gp = gpar(fill = c('#5852ADFF','#879ECEFF','#00008BFF','#AFEEEEFF')), 
                         labels = c('DZ50','TGFB','TGFB_DZ50','Untreated'),
                         labels_gp = gpar(col = "black", fontsize = 10)) )

### ann_row$Type <- factor(x = ann_row$Type,levels = c('Lipoproteins','Adhesion','Hallmark')) ##

ann_col_v2$Treatment <- factor(ann_col_v2$Treatment,levels = c('DZ50','TGFB','TGFB_DZ50','Untreated'))

pdf(file='figures_TD006643/heatmap_apoptosis_VCaP.pdf',width = 6, height = 10)
draw(
  ComplexHeatmap::Heatmap(z, 
                          col= colorRampPalette(c("steelblue", "white","firebrick"))(255),
                          show_row_names = TRUE , show_column_names = FALSE,
                          #column_title_gp = gpar(fontsize = 1, fontface = "plain"), #column_names_gp = gpar(fontsize = 5, fontface = "plain"),
                          column_names_rot = 45,
                          column_split = ann_col_v2$Treatment, #column_km = 7,
                          #row_split =  ann_row, #row_split = 4, #row_km = 5,
                          clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
                          cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                          top_annotation = column_ha, 
                          #left_annotation = row_ha,
                          border = TRUE, use_raster = FALSE, name = "LogFC") ,
  
  heatmap_legend_side = "left", annotation_legend_side = "left" , padding = unit(c(2, 2, 2, 20), "mm")
)
dev.off()
#############################################################################################################################################

#############################################################################################################################################
#FIGURE 2
#> DZ50 kills them regardless of TGFB. (to be added)
#FIGURE 2.
#1. Figures with galactXXX inh of TGFB.
#2. Add figures for:
### What is different between VCaP vs TRII.
### VCaP are more sensitive. Die regardless. Sinergy effect of TGFB+DZ50.
results_analysis$de_v1$Comparison <- make.names(results_analysis$de_v1$Comparison)
###
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....TRII.Untreated.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB.....TRII.TGFB.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.DZ50.....TRII.DZ50.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB_DZ50.....TRII.TGFB_DZ50.')

### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....TRII.Untreated.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_untreated_VCaP_vs_TRII.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#FFBBFFFF','#AFEEEEFF'))
dev.off()
#############################################################################################################################################
library(org.Hs.eg.db)
### get entrezids on the fly
goana_df <- goana(results_analysis$fitmm_v1, coef = 1, 
                  geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                  FDR = 0.05)
###
write.csv(file='RData/goana_untreated_vcap_vs_trii.csv',goana_df)
###
for_figure_pathway <- rbind(goana_df[order(goana_df$P.Up,decreasing = F),][1:10,],
                            goana_df[order(goana_df$P.Down,decreasing = F),][1:10,])

for_figure_pathway$direction <- c(rep('VCaP',10),rep('TRII',10))
for_figure_pathway$score <- c(-log10(for_figure_pathway$P.Up[1:10]) , -(-log10(for_figure_pathway$P.Down[11:20]) ) )

### pathway
pdf(file='figures_TD006643/barplot_GO_untreated_VCaP_vs_TRII.pdf',width = 5,height = 4)
ggplot(data=for_figure_pathway) + aes(x=scale(score),y=reorder(Term,score),fill=direction) + theme_classic()+
  geom_bar(stat='identity',show.legend = F)+ 
  scale_fill_manual(values = c('#FFBBFFFF','#AFEEEEFF'))
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB.....TRII.TGFB.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_TGFB_VCaP_vs_TRII.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#E082A4FF','#879ECEFF'))
dev.off()
#############################################################################################################################################
### get entrezids on the fly
goana_df <- goana(results_analysis$fitmm_v1, coef = 3, 
                  geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                  FDR = 0.05)
###
write.csv(file='RData/goana_TGFB_vcap_vs_trii.csv',goana_df)
###
for_figure_pathway <- rbind(goana_df[order(goana_df$P.Up,decreasing = F),][1:10,],
                            goana_df[order(goana_df$P.Down,decreasing = F),][1:10,])

for_figure_pathway$direction <- c(rep('VCaP',10),rep('TRII',10))
for_figure_pathway$score <- c(-log10(for_figure_pathway$P.Up[1:10]) , -(-log10(for_figure_pathway$P.Down[11:20]) ) )

for_figure_pathway$Term[for_figure_pathway$Term %in% "steroid dehydrogenase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"] <- "dehydrogenase ADP as acceptor"

### pathway
pdf(file='figures_TD006643/barplot_GO_TGFB_VCaP_vs_TRII.pdf',width = 5,height = 4)
ggplot(data=for_figure_pathway) + aes(x=scale(score),y=reorder(Term,score),fill=direction) + theme_classic()+
  geom_bar(stat='identity',show.legend = F)+ 
  scale_fill_manual(values = c('#E082A4FF','#879ECEFF'))
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.VCaP.DZ50.....TRII.DZ50.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_dz50_VCaP_vs_TRII.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#BA4A50FF','#5852ADFF'))
dev.off()
#############################################################################################################################################
### get entrezids on the fly
goana_df <- goana(results_analysis$fitmm_v1, coef = 2, 
                  geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                  FDR = 0.05)
###
write.csv(file='RData/goana_dz50_vcap_vs_trii.csv',goana_df)
###
for_figure_pathway <- rbind(goana_df[order(goana_df$P.Up,decreasing = F),][1:10,],
                            goana_df[order(goana_df$P.Down,decreasing = F),][1:10,])

for_figure_pathway$direction <- c(rep('VCaP',10),rep('TRII',10))
for_figure_pathway$score <- c(-log10(for_figure_pathway$P.Up[1:10]) , -(-log10(for_figure_pathway$P.Down[11:20]) ) )

### pathway
pdf(file='figures_TD006643/barplot_GO_dz50_VCaP_vs_TRII.pdf',width = 5,height = 4)
ggplot(data=for_figure_pathway) + aes(x=scale(score),y=reorder(Term,score),fill=direction) + theme_classic()+
  geom_bar(stat='identity',show.legend = F)+ 
  scale_fill_manual(values = c('#BA4A50FF','#5852ADFF'))
dev.off()


#############################################################################################################################################
### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB_DZ50.....TRII.TGFB_DZ50.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_dz50+tgfb_VCaP_vs_TRII.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#8B0000FF','#00008BFF'))
dev.off()
#############################################################################################################################################
### get entrezids on the fly
goana_df <- goana(results_analysis$fitmm_v1, coef = 4, 
                  geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                  FDR = 0.05)
###
write.csv(file='RData/goana_dz50+tgfb_vcap_vs_trii.csv',goana_df)
###
for_figure_pathway <- rbind(goana_df[order(goana_df$P.Up,decreasing = F),][1:10,],
                            goana_df[order(goana_df$P.Down,decreasing = F),][1:10,])
###
for_figure_pathway$direction <- c(rep('VCaP',10),rep('TRII',10))
for_figure_pathway$score <- c(-log10(for_figure_pathway$P.Up[1:10]) , -(-log10(for_figure_pathway$P.Down[11:20]) ) )
### pathway
pdf(file='figures_TD006643/barplot_GO_dz50+tgfb_VCaP_vs_TRII.pdf',width = 5,height = 4)
ggplot(data=for_figure_pathway) + aes(x=scale(score),y=reorder(Term,score),fill=direction) + theme_classic()+
  geom_bar(stat='identity',show.legend = F)+ 
  scale_fill_manual(values = c('#8B0000FF','#00008BFF'))
dev.off()
#############################################################################################################################################


#############################################################################################################################################
###
#############################################################################################################################################
'CFL1'
COFILIN

table(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....TRII.Untreated.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB.....TRII.TGFB.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.DZ50.....TRII.DZ50.')
table(results_analysis$de_v1$Comparison %in% 'X.VCaP.TGFB_DZ50.....TRII.TGFB_DZ50.')

results_analysis$de_v1[which(results_analysis$de_v1$symbol %in% 'CFL1' & results_analysis$de_v1$adj.P.Val<0.05),]

identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)

results_analysis$alignment_stats$condition <-paste(results_analysis$alignment_stats$treatment,results_analysis$alignment_stats$cell_line,sep='___')

my_condition <- unique(results_analysis$alignment_stats$condition)

my_cor_results <- list()

for(i in 1:length(my_condition)){
  
  tmp_df <- results_analysis$voom_TD006643$E[,which(results_analysis$alignment_stats$condition %in% my_condition[i] )]
  tmp_df <- tmp_df[ -grep("^NA..",rownames(tmp_df)) , ]
  tmp_df <- tmp_df[ -grep("NA",rownames(tmp_df)) , ]
  tmp_df <- tmp_df[ !is.na(rownames(tmp_df)) , ]
  
  my_cor <- cor(t(tmp_df))
  my_cor <- melt(my_cor)
  my_cor <- my_cor[abs(my_cor$value)>0.8,]
  my_cor <- my_cor[which(my_cor$Var1 %in% 'CFL1'),]
  my_cor <- my_cor[order(abs(my_cor$value),decreasing = T),]
  my_cor <- my_cor[c(which(my_cor$value>0)[1:11],which(my_cor$value<0)[1:10]),]
  
  my_cor_results[[i]] <- my_cor
}

for(i in 1:length(my_condition)){my_cor_results[[i]]$comparison <- my_condition[i]}

my_cor <- do.call(rbind,my_cor_results)
my_cor <- as.character(unique(c(my_cor$Var1,my_cor$Var2)))
my_cor <- my_cor[my_cor %in% unique(results_analysis$de_v1$symbol[results_analysis$de_v1$adj.P.Val<0.05])]

summary(rowMeans(results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% my_cor,]))
summary(rowMedians(results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% my_cor,]) )

x <- results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% my_cor,][which(rowMedians(results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% my_cor,]) > -2),]

colnames(x)[6]
colnames(x)[9]
x <- x[,c(1,2,3,4,5,9,7,8,6,10,11,12:ncol(x))]



ann_col <- data.frame(TreatCell=paste(results_analysis$alignment_stats$treatment,results_analysis$alignment_stats$cell_line))
rownames(ann_col) <- colnames(results_analysis$voom_TD006643$E)
                                                                                            
ann_colors<-list(TreatCell=c(
  "DZ50 TRII"='#BA4A50FF',
  "TGFB TRII"='#E082A4FF',
  "TGFB_DZ50 TRII"='#8B0000FF',
  "Untreated TRII"='#FFBBFFFF',
  "Untreated VCaP"='#AFEEEEFF',
  "DZ50 VCaP"='#5852ADFF',
  "TGFB VCaP"='#879ECEFF',
  "TGFB_DZ50 VCaP"='#00008BFF'),
  Cor=c(`1`='steelblue',`2`='firebrick'))
  

out <- pheatmap::pheatmap(t(x),angle_col = 90,cluster_rows = F,
                    col= colorRampPalette(c("steelblue", "white","firebrick"))(25),clustering_method = 'complete',
                    annotation_row = ann_col,annotation_colors = ann_colors,
                    scale = 'column')

ann_row <- data.frame(Cor=as.character(cutree(out$tree_col,2)))
rownames(ann_row) <- out$tree_col$labels

pdf(file='figures_TD006643/heatmap_genes_correlated_to_CFL1.pdf',width = 17,height = 6)
pheatmap::pheatmap(t(x),angle_col = 90,cluster_rows = F,
                   col= colorRampPalette(c("steelblue", "white","firebrick"))(25),clustering_method = 'complete',
                     annotation_row = ann_col,annotation_colors = ann_colors,annotation_col = ann_row,
                   scale = 'column')
dev.off()

y <- as.data.frame(read_tsv(file = 'string_interactions_short.tsv'))

write.csv(file='genes_correlated_to_coff.csv',ann_row)
### rownames(x) %in% unique(c(y$`#node1`,y$node2))


z <- results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% unique(c(y$`#node1`,y$node2)),]

pdf(file='figures_TD006643/heatmap_PPI_CFL1.pdf',width = 6,height = 2.9)
pheatmap::pheatmap((z),angle_col = 90,cluster_rows = T,show_colnames = F,
                   col= colorRampPalette(c("steelblue", "white","firebrick"))(25),clustering_method = 'complete',
                   annotation_col = ann_col,annotation_colors = ann_colors,#annotation_col = ann_row,
                   scale = 'row')
dev.off()
#############################################################################################################################################

#############################################################################################################################################
#Intersection of EMT and anoikis. (to be added)
#Figure 5. EMT. blots
#TRII response to TGFb in EMT markers.
#############################################################################################################################################
EMT_markers <- list(
  epithelial=c('CDH1', 'DSP', 'OCLN', 'CRB3'),
  mesenchymal=c('VIM', 'CDH2', 'FOXC2', 'SNAI1', 'SNAI2','SNAI3', 'TWIST1','FN1', 'ITGB6', 'MMP2', 'MMP3', 'MMP9', 'SOX10', 'GSC', 'ZEB1', 'ZEB2','TWIST2'))


z <- results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% unique(unlist(EMT_markers)),]
#pdf(file='figures_TD006643/heatmap_EMT_genes.pdf',width = 6,height = 2.9)
pheatmap::pheatmap((z),angle_col = 90,cluster_rows = T,show_colnames = F,
                   col= colorRampPalette(c("steelblue", "white","firebrick"))(25),clustering_method = 'complete',
                   annotation_col = ann_col,annotation_colors = ann_colors,#annotation_col = ann_row,
                   scale = 'row')
dev.off()
###
z<-melt(z)
###
ann_col$Var2 <- rownames(ann_col)
###
z <- merge(z,ann_col,by='Var2')

z$EMT <- 'Mesenchymal'
z$EMT[which(z$Var1 %in% c('CDH1', 'DSP', 'OCLN', 'CRB3'))] <- 'Epithelial'
head(z)

z$cell <- 'TRII'
z$cell[grep('VCaP',z$TreatCell)] <- 'VCaP'

z$treat <- z$TreatCell
z$treat <- gsub(' TRII| VCaP','',z$treat)

my_comparisons <- list(c('Untreated','TGFB'),c('Untreated','DZ50'),c('Untreated','TGFB_DZ50'))
my_comparisons2 <- my_comparisons

pdf(file='figures_TD006643/boxplots_EMT_genes_for_figure4_a.pdf',width = 10,height = 3.25)
ggplot(data=z[z$cell %in% 'TRII',])+aes(x=treat,y=value,fill=TreatCell) + theme_classic() + 
  facet_wrap(~cell+Var1,ncol=9,scales = 'free_x') + geom_boxplot(show.legend = F) + rotate_x_text() +
  scale_fill_manual(values = list("DZ50 TRII"='#BA4A50FF',"TGFB TRII"='#E082A4FF',"TGFB_DZ50 TRII"='#8B0000FF',"Untreated TRII"='#FFBBFFFF',
                                  "Untreated VCaP"='#AFEEEEFF',"DZ50 VCaP"='#5852ADFF',"TGFB VCaP"='#879ECEFF',"TGFB_DZ50 VCaP"='#00008BFF')) +
  stat_compare_means(comparisons = my_comparisons,method = 't.test',paired = T) +
    stat_compare_means(comparisons = my_comparisons2,method = 't.test',paired = T)
dev.off()

pdf(file='figures_TD006643/boxplots_EMT_genes_for_figure4_b.pdf',width = 10,height = 3.25)
ggplot(data=z[z$cell %in% 'VCaP',])+aes(x=treat,y=value,fill=TreatCell) + theme_classic() + 
  facet_wrap(~cell+Var1,ncol=9,scales = 'free_x') + geom_boxplot(show.legend = F) + rotate_x_text() +
  scale_fill_manual(values = list("DZ50 TRII"='#BA4A50FF',"TGFB TRII"='#E082A4FF',"TGFB_DZ50 TRII"='#8B0000FF',"Untreated TRII"='#FFBBFFFF',
                                  "Untreated VCaP"='#AFEEEEFF',"DZ50 VCaP"='#5852ADFF',"TGFB VCaP"='#879ECEFF',"TGFB_DZ50 VCaP"='#00008BFF')) +
  stat_compare_means(comparisons = my_comparisons,method = 't.test',paired = T) +
  stat_compare_means(comparisons = my_comparisons2,method = 't.test',paired = T)
dev.off()

#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################


#############################################################################################################################################
### Volcanos per treatment and signature intersection
#############################################################################################################################################
### At the transcriptome level, there was little effect of DZ50 on TRII.
unique(results_analysis$de_v1$Comparison)

### Venns
x <- list(
DZ50=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.TRII.Untreated.....TRII.DZ50.' & results_analysis$de_v1$adj.P.Val<0.05)],
TFGB=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.TRII.Untreated.....TRII.TGFB.' & results_analysis$de_v1$adj.P.Val<0.05)],
DZ50_TGFB=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.TRII.Untreated.....TRII.TGFB_DZ50.' & results_analysis$de_v1$adj.P.Val<0.05)])

pdf(file='figures_TD006643/venn_intersection_TRII_DEGs.pdf',width = 7,height = 7)
plot(Venn(x))
dev.off()

x <- list(
  DZ50=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....VCaP.DZ50.' & results_analysis$de_v1$adj.P.Val<0.05)],
  TFGB=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....VCaP.TGFB.' & results_analysis$de_v1$adj.P.Val<0.05)],
  DZ50_TGFB=results_analysis$de_v1$ensembl[which(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....VCaP.TGFB_DZ50.' & results_analysis$de_v1$adj.P.Val<0.05)])

pdf(file='figures_TD006643/venn_intersection_VCaP_DEGs.pdf',width = 7,height = 7)
plot(Venn(x))
dev.off()

### Volcanos
###  results_analysis$de_v1

unique(results_analysis$de_v1$Comparison)

results_analysis$color_pal

### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.TRII.Untreated.....TRII.TGFB_DZ50.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_untreated_TRII_vs_combination.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#8B0000FF','#FFBBFFFF'))
dev.off()


### subset
ix <- which(results_analysis$de_v1$Comparison %in% 'X.VCaP.Untreated.....VCaP.TGFB_DZ50.')
###
tmp_df <- results_analysis$de_v1[ix,]
### only known genes
tmp_df <- tmp_df[!is.na(tmp_df$symbol),]
###
tmp_df$color <- ''
tmp_df$color[tmp_df$logFC>1 & tmp_df$adj.P.Val<0.05] <- 'VCaP'
tmp_df$color[tmp_df$logFC < -1 & tmp_df$adj.P.Val<0.05] <- 'TRII'
###
ix <- c(which(tmp_df$logFC>0)[1:11],which(tmp_df$logFC<0)[1:11])
### volcano
pdf(file='figures_TD006643/volcano_plot_top15_untreated_VCaP_vs_combination.pdf',width = 6,height = 4)
ggplot(data=tmp_df)+aes(x=logFC,y=-log10(adj.P.Val),color=color)+theme_classic()+
  geom_point(show.legend = F) + 
  xlim(-12.2,12.2)+
  geom_vline(xintercept = -1,linetype=2,color='red')+
  geom_vline(xintercept = 1,linetype=2,color='red')+
  geom_hline(yintercept = 1.3,linetype=2,color='red') +
  geom_text_repel(data=tmp_df[ix,],aes(label=symbol),color='black',max.overlaps = 100,force = 50,force_pull = 20) +
  scale_color_manual(values = c('grey50','#00008BFF','#AFEEEEFF'))
dev.off()

### pathway venn
library(org.Hs.eg.db)

### TRII

### get entrezids on the fly
dz50 <- goana(results_analysis$fitmm_v1, coef = 11, 
                  geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                  FDR = 0.05)

### get entrezids on the fly
tgfb <- goana(results_analysis$fitmm_v1, coef = 12, 
              geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
              FDR = 0.05)

### get entrezids on the fly
dz50_tgfb <- goana(results_analysis$fitmm_v1, coef = 13, 
              geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
              FDR = 0.05)

### store
pathway_results_trii<-list(dz50=dz50,tgfb=tgfb,dz50_tgfb=dz50_tgfb)

x <- list(
  DZ50=dz50$Term[dz50$N > 4 & dz50$P.Up<0.05 | dz50$N > 4 & dz50$P.Down<0.05],
  TGFB=tgfb$Term[tgfb$N > 4 & tgfb$P.Up<0.05 | tgfb$N > 4 & tgfb$P.Down<0.05],
  DZ50_TGFB=dz50_tgfb$Term[dz50_tgfb$N > 4 & dz50_tgfb$P.Up<0.05 | dz50_tgfb$N > 4 & dz50_tgfb$P.Down<0.05])

pdf(file='figures_TD006643/venn_intersection_TRII_pathways.pdf',width = 7,height = 7)
plot(Venn(x))
dev.off()

### Reduce(intersect,x) nothing
pdf(file='figures_TD006643/venn_intersection_TRII_pathways_upset.pdf',width = 7,height = 7)
UpSetR::upset(fromList(x))
dev.off()

m_trii = make_comb_mat(x)
extract_comb(m_trii, "011")
extract_comb(m_trii, "101")
extract_comb(m_trii, "110")

### VCap

### get entrezids on the fly
dz50 <- goana(results_analysis$fitmm_v1, coef = 5, 
              geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
              FDR = 0.05)

### get entrezids on the fly
tgfb <- goana(results_analysis$fitmm_v1, coef = 6, 
              geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
              FDR = 0.05)

### get entrezids on the fly
dz50_tgfb <- goana(results_analysis$fitmm_v1, coef = 7, 
                   geneid = as.character(mapIds(org.Hs.eg.db, keys = as.character(rownames(results_analysis$fitmm_v1)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) ,
                   FDR = 0.05)

### store
pathway_results_vcap<-list(dz50=dz50,tgfb=tgfb,dz50_tgfb=dz50_tgfb)

y <- list(
  DZ50=dz50$Term[dz50$N > 4 & dz50$P.Up<0.05 | dz50$N > 4 & dz50$P.Down<0.05],
  TGFB=tgfb$Term[tgfb$N > 4 & tgfb$P.Up<0.05 | tgfb$N > 4 & tgfb$P.Down<0.05],
  DZ50_TGFB=dz50_tgfb$Term[dz50_tgfb$N > 4 & dz50_tgfb$P.Up<0.05 | dz50_tgfb$N > 4 & dz50_tgfb$P.Down<0.05])

pdf(file='figures_TD006643/venn_intersection_VCaP_pathways.pdf',width = 7,height = 7)
plot(Venn(y))
dev.off()

### Reduce(intersect,x) nothing
pdf(file='figures_TD006643/venn_intersection_VCaP_pathways_upset.pdf',width = 7,height = 7)
UpSetR::upset(fromList(y))
dev.off()

m_vcap = make_comb_mat(y)
m_vcap
extract_comb(m_vcap, "011")
extract_comb(m_vcap, "110")



Reduce(intersect,list(a=extract_comb(m_vcap, "101"),b=c(extract_comb(m_trii, "011"),extract_comb(m_trii, "101"))))


c(extract_comb(m_trii, "011"),extract_comb(m_trii, "101"))

pathways_combined <- c(c(pathway_results_trii$dz50_tgfb$Term[pathway_results_trii$dz50_tgfb$N > 4 & pathway_results_trii$dz50_tgfb$P.Down<0.005 | pathway_results_trii$dz50_tgfb$N > 4 & pathway_results_trii$dz50_tgfb$P.Up<0.005]),
c(pathway_results_vcap$dz50_tgfb$Term[pathway_results_vcap$dz50_tgfb$N > 4 & pathway_results_vcap$dz50_tgfb$P.Down<0.001 | pathway_results_vcap$dz50_tgfb$N > 4 & pathway_results_vcap$dz50_tgfb$P.Up<0.001]))
pathways_combined <- unique(pathways_combined)

x <- rbind(do.call(rbind,pathway_results_trii),do.call(rbind,pathway_results_trii))

x <- x[which(x$Term %in% pathways_combined),]

x <- x %>% group_by(Term,N,P.Up,P.Down) %>% summarise()



###
for_figure_pathway <- rbind(goana_df[order(goana_df$P.Up,decreasing = F),][1:10,],
                            goana_df[order(goana_df$P.Down,decreasing = F),][1:10,])
###
for_figure_pathway$direction <- c(rep('VCaP',10),rep('TRII',10))
for_figure_pathway$score <- c(-log10(for_figure_pathway$P.Up[1:10]) , -(-log10(for_figure_pathway$P.Down[11:20]) ) )
### pathway
pdf(file='figures_TD006643/barplot_GO_dz50+tgfb_VCaP_vs_TRII.pdf',width = 5,height = 4)
ggplot(data=for_figure_pathway) + aes(x=scale(score),y=reorder(Term,score),fill=direction) + theme_classic()+
  geom_bar(stat='identity',show.legend = F)+ 
  scale_fill_manual(values = c('#8B0000FF','#00008BFF'))
dev.off()


#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################

### Questions Aug 19 2024
#############################################################################################################################################

dz50 <- pathway_results_trii$dz50
dz50_tgfb <- pathway_results_trii$dz50_tgfb

path_isec <- extract_comb(m_trii, "101")

trii_DZ50 <- dz50[dz50$N > 4 & dz50$P.Up<0.05 | dz50$N > 4 & dz50$P.Down<0.05,]
trii_DZ50_TGFB <- dz50_tgfb[dz50_tgfb$N > 4 & dz50_tgfb$P.Up<0.05 | dz50_tgfb$N > 4 & dz50_tgfb$P.Down<0.05,]

trii_DZ50 <- trii_DZ50[trii_DZ50$Term %in% path_isec,]
trii_DZ50_TGFB <- trii_DZ50_TGFB[trii_DZ50_TGFB$Term %in% path_isec,]  


trii_DZ50$zp <- -log10(trii_DZ50$P.Up)
trii_DZ50$zp[trii_DZ50$Down > 0] <- -log10(trii_DZ50$P.Down)[trii_DZ50$Down > 0]

trii_DZ50_TGFB$zp <- -log10(trii_DZ50_TGFB$P.Up)
trii_DZ50_TGFB$zp[trii_DZ50_TGFB$Down > 0] <- -log10(trii_DZ50_TGFB$P.Down)[trii_DZ50_TGFB$Down > 0]
trii_DZ50$Treat <- 'DZ50'
trii_DZ50_TGFB$Treat <- 'DZ50+TGFB'

trii_DZ50_TGFB$zp <- -trii_DZ50_TGFB$zp

big_table <- rbind(trii_DZ50_TGFB,trii_DZ50)

head(big_table)

pdf(file='figures_TD006643/overlap_btw_pathways_trii_prerna.pdf',width = 10,height = 5)
ggplot(data=big_table)+aes(y=reorder(Term,zp,median),x=zp,fill=Treat)+theme_classic()+geom_bar(stat='identity') +
scale_fill_manual(values = as.character(color_ann$LiTr[which(names(color_ann$LiTr) %in% c("TRII DZ50","TRII TGFB_DZ50"))]) )
dev.off()

pdf(file='figures_TD006643/overlap_btw_pathways_trii_prerna_n_size.pdf',width = 5,height = 5)
ggplot(data=big_table)+aes(y=reorder(Term,zp,median),x=N/2)+theme_classic()+geom_bar(stat='identity')
dev.off()


library(biomaRt)
library(GO.db)
library(org.Hs.eg.db)

retrieved <- list()
for(i in 1:length(rownames(trii_DZ50_TGFB))){
  retrieved[[i]] <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=rownames(trii_DZ50_TGFB)[1], columns="SYMBOL")  
  retrieved[[i]] <- unique(retrieved[[i]]$SYMBOL)
}


#table((table(retrieved$SYMBOL)>2))
#table(retrieved$SYMBOL)[table(retrieved$SYMBOL)>1]
#unique(retrieved$SYMBOL[])
#colSums(table(retrieved$SYMBOL,retrieved$EVIDENCEALL))

# LOGFC LOGFC inersetction


x <- list(
  DZ50=results_analysis$de_v1[results_analysis$de_v1$Comparison %in% "X.TRII.Untreated.....TRII.DZ50.",],
  DZ50_TGFB=results_analysis$de_v1[results_analysis$de_v1$Comparison %in% "X.TRII.Untreated.....TRII.TGFB_DZ50.",])

x$DZ50 <- x$DZ50$symbol[x$DZ50$P.Value<0.005]
x$DZ50_TGFB <- x$DZ50_TGFB$symbol[x$DZ50_TGFB$P.Value<0.005]

x$DZ50_TGFB <- x$DZ50_TGFB[!is.na(x$DZ50_TGFB)]
x$DZ50 <- x$DZ50[!is.na(x$DZ50)]

table(my_genes %in% unlist(x) )

pdf(file='figures_TD006643/venn_intersection_TRII_genes.pdf',width = 7,height = 7)
plot(Venn(x))
dev.off()


dz50_tgfb <- results_analysis$de_v1[results_analysis$de_v1$Comparison %in% "X.TRII.Untreated.....TRII.TGFB_DZ50.",]
dz50 <- results_analysis$de_v1[results_analysis$de_v1$Comparison %in% "X.TRII.Untreated.....TRII.DZ50.",]

dz50_tgfb<-dz50_tgfb[!is.na(dz50_tgfb$symbol),]
dz50<-dz50[!is.na(dz50$symbol),]

identical(dz50$symbol,dz50_tgfb$symbol)
table(dz50$symbol %in% dz50_tgfb$symbol)
dz50_tgfb <- dz50_tgfb[match(dz50$symbol,dz50_tgfb$symbol),]
identical(dz50$symbol,dz50_tgfb$symbol)

x <- data.frame(logFC_dz50_tgfb = dz50_tgfb$logFC,
           fdr_dz50_tgfb = dz50_tgfb$P.Value,
           symbol_dz50_tgfb = dz50_tgfb$symbol,
           logFC_dz50 = dz50$logFC,
           fdr_dz50 = dz50$P.Value,
           symbol_dz50 = dz50$symbol)
x <- x[!is.na(x$symbol_dz50),]

x$fdr <- 'ns'
x$fdr[x$fdr_dz50_tgfb<0.05 ] <- 'DZ50+TGFB'
x$fdr[x$fdr_dz50<0.005] <- 'DZ50'
x$fdr[x$fdr_dz50<0.005 & x$fdr_dz50_tgfb<0.005] <- 'Intersection'
x$fdr[x$logFC_dz50 > 0 | x$logFC_dz50_tgfb > 0] <- 'Untreated'

x$fdr <- factor(x$fdr,levels = c('ns','Untreated','DZ50+TGFB','DZ50','Intersection'))
table(x$fdr)

pdf(file='figures_TD006643/overlap_btw_logfc_trii_prerna_n_size.pdf',width = 7,height = 5)
ggplot()+#data=x
  aes(x=-logFC_dz50,y=-logFC_dz50_tgfb,color=fdr)+
  theme_classic()+
  geom_point(data=x[!x$fdr %in% 'Intersection',],alpha=0.5,shape=16)+
  geom_point(data=x[x$fdr %in% 'Intersection',],alpha=1,shape=16)+
  geom_hline(yintercept = 0,linetype=2,color='red')+
  geom_vline(xintercept = 0,linetype=2,color='red')+
  scale_color_manual(values = c('grey30','grey70','grey80','grey90','red')) +
  ylim(c(-2,6)) + xlim(c(-2,6)) + 
  geom_text_repel(data=x[x$fdr %in% 'Intersection',],aes(label=symbol_dz50),show.legend = F,max.overlaps = 25,force = 7,force_pull = 1)
dev.off()
  
my_genes <- x$symbol_dz50[x$fdr %in% 'Intersection']



identical(colnames(results_analysis$voom_TD006643$E),results_analysis$alignment_stats$sampleID)

df_ann <- data.frame(#Line=results_analysis$alignment_stats$cell_line,
  #Treat=results_analysis$alignment_stats$treatment,
  LiTr=results_analysis$alignment_stats$individualID)
rownames(df_ann) <- results_analysis$alignment_stats$sampleID

color_ann <- list(LiTr=unlist(results_analysis$color_pal))
#Line = c(TRII='#FFBBFFFF',VCaP='#AFEEEEFF'),
#Treat = c(TRII='#FFBBFFFF',VCaP='#AFEEEEFF'),)

iy <- grep('TRII',df_ann$LiTr,ignore.case = T)

pdf(file='figures_TD006643/heatmap_intersection_markers_trii.pdf',width = 5,height = 6.5)
pheatmap::pheatmap(mat = results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% my_genes,iy],
                   scale = 'row',
                   angle_col= 90,
                   clustering_method = 'mcquitty',
                   annotation_col = df_ann,
                   annotation_colors = color_ann)
dev.off()


x <- as.data.frame(results_analysis$voom_TD006643$E[results_analysis$voom_TD006643$genes$symbol %in% c('PARP1'),])
colnames(x) <- 'value'
x$samples <- rownames(x)

df_ann$samples <- rownames(df_ann)

x <- merge(df_ann,x)
library(ggplot2)

pdf(file='fig_parp1.pdf',width = 4,height = 3)
ggplot(data=x) + aes(x=LiTr,y=value) + geom_boxplot() + rotate_x_text(angle = 45) + labs(x='Condition',y='Expression',title = 'PARP1')
dev.off()
#############################################################################################################################################



#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################



#Edgar, after my discussions with Prerna, we hope you can re-do one of the Figures as the heatmaps are not in vertical orientation. 
#This is for Figure 1, panels e and f, that be more  “read-ready” for the readers is presented vertical with names of genes listed 
#on the left (instead of below as they are now)and treatment listed on top.


#############################################################################################################################################
#Figure 1 e.
#############################################################################################################################################



#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image(file = 'RData/analysis_pr_TD006643.RData')
#load(file = 'RData/analysis_pr_TD006643.RData')
#############################################################################################################################################