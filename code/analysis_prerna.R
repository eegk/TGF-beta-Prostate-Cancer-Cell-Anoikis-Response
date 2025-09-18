### libraries and functions
#############################################################################################################################################
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
library(enrichR)
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
#############################################################################################################################################

### Work space
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_dogra/projects_kyprianou/prerna_natasha/")
#############################################################################################################################################

#############################################################################################################################################
### get log (in minerva)
#############################################################################################################################################
temp = list.files(pattern="^Log.final.out$",recursive = TRUE,full.name = TRUE)
temp <- temp[-grep('_STARpass1',temp)]
my_log_files = lapply(temp, read.table, sep="\t",header=TRUE,fill=TRUE) 
names(my_log_files) <- temp
save.image("logs.RData")
#############################################################################################################################################

#############################################################################################################################################
### stats
#############################################################################################################################################
load('RData/logs.RData')
###
names(my_log_files) <- gsub('./align_star/','',names(my_log_files))
names(my_log_files) <- gsub('.fastq.trimmed_fastq/Log.final.out','',names(my_log_files))
names(my_log_files) <- gsub('\\.\\/','',names(my_log_files))
names(my_log_files) <- make.names(names(my_log_files))
names(my_log_files) <- gsub('\\.','_',names(my_log_files))
###

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
temp_params <- melt(as.matrix(alignment_stats))
### alignment stats figure
pdf(file="figures/alignment_statistics_TD06527.pdf",width = 3,height = 3)
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
counts_TD06527 <- read.table(file='count_matrices_TD06527/counts.all.txt',header = T)
#table(counts_TD06527$gene_biotype)
results_analysis$raw_counts_TD06527 <- counts_TD06527 ; rm(counts_TD06527)
colnames(results_analysis$alignment_stats)
my_names <- colnames(results_analysis$raw_counts_TD06527)[13:22] 
my_names <- gsub('\\.','_',my_names)
my_names <- gsub('_fastq_trimmed_fastq_Aligned_out_bam','',my_names)
my_names <- gsub('__','',my_names)
table(colnames(results_analysis$alignment_stats) %in% my_names)
colnames(results_analysis$raw_counts_TD06527)[13:22]  <- my_names
###
results_analysis$counts_TD06527 <- results_analysis$raw_counts_TD06527[,13:22]
rownames(results_analysis$counts_TD06527) <- results_analysis$raw_counts_TD06527$Geneid
results_analysis$annotation_TD06527 <- results_analysis$raw_counts_TD06527[,1:12]
rm(my_names)
#############################################################################################################################################


### normalization
#############################################################################################################################################
### investigate
table(rowSums(cpm(results_analysis$counts_TD06527)>1)>=1)
table(rowSums(cpm(results_analysis$counts_TD06527)>1)>=3)
table(rowSums(cpm(results_analysis$counts_TD06527)>1)>=5)
table(rowSums(cpm(results_analysis$counts_TD06527)>5)>=3)
table(rowSums(cpm(results_analysis$counts_TD06527)>5)>=5)
table(rowSums(cpm(results_analysis$counts_TD06527)>9)>=5)
###
keep.counts <- which(rowSums(cpm(results_analysis$counts_TD06527)>=5)>=3)
###
deglist_obj <- DGEList(counts = results_analysis$counts_TD06527[keep.counts,],
                       genes = results_analysis$annotation_TD06527[keep.counts,])
###
deglist_obj <- calcNormFactors(deglist_obj, method = "TMM") 
deglist_obj = estimateCommonDisp(deglist_obj, verbose=TRUE)
### voom
results_analysis$voom_TD06527 <- voom(deglist_obj, plot = TRUE, save.plot = TRUE) 
results_analysis$voom_TD06527$pseudocounts <- deglist_obj$pseudo.counts
results_analysis$voom_TD06527$samples <- deglist_obj$samples
dev.off()
rm(deglist_obj)
#############################################################################################################################################


#############################################################################################################################################
### naming
#############################################################################################################################################
colnames(results_analysis$voom_TD06527)
results_analysis$voom_TD06527$samples$labels <- rownames(results_analysis$voom_TD06527$samples)
results_analysis$voom_TD06527$samples$condition <- tidyr::separate(data.frame(results_analysis$voom_TD06527$samples$labels),1, sep="_", c("a","b",'c','d'))$a
results_analysis$voom_TD06527$samples$condition[grep('TGFB',results_analysis$voom_TD06527$samples$condition)] <- c("TGFB_DZ50")
#############################################################################################################################################


#############################################################################################################################################
#ix <- which(results_analysis$voom_TD06527$samples$condition %in% 'Undetermined')
#pheatmap::pheatmap(results_analysis$voom_TD06527$E[,-ix],scale = 'row')
#dev.off()
###
identical(results_analysis$voom_TD06527$samples$labels,colnames(results_analysis$voom_TD06527$E))
###
library(org.Hs.eg.db)
results_analysis$voom_TD06527$genes$entrez <- as.character(mapIds(org.Hs.eg.db, keys = as.character(results_analysis$voom_TD06527$genes$Geneid), column="ENTREZID", keytype="ENSEMBL", multiVals="first")) 
###
design = model.matrix( ~ 0 + condition, data = results_analysis$voom_TD06527$samples )
### names
colnames(design) <- gsub("condition","",colnames(design))
### contr.matrix
contr.matrix <- makeContrasts( 'DZ50_vs_Un' = (DZ50 - Untreated),
                               'TGFB_DZ50_vs_Un' = (TGFB_DZ50 - Untreated),
                               'DZ50_vs_TGFB_DZ50' = (DZ50 - TGFB_DZ50),
                                levels = colnames(design))
### make model
vfit <- lmFit(results_analysis$voom_TD06527, design)
### fit contrasts
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
### bayesian statistics
efit <- eBayes(vfit) ; summary(decideTests(efit))

### loop to catch the results
lmfreq_results <- list()
for ( xzy in 1:ncol(summary(decideTests(efit))) ) { 
  lmfreq_results[[xzy]] <- topTable(efit, coef=xzy, n=Inf, adjust.method="BH") 
  lmfreq_results[[xzy]]$Comparison <- colnames(summary(decideTests(efit)))[xzy]
  lmfreq_results[[xzy]]$Marker <- rownames( topTable(efit, coef=xzy,n=Inf,adjust.method="BH") ) 
  rownames(lmfreq_results[[xzy]]) <- NULL }      
###
de_analysis_results <- do.call(rbind,lmfreq_results)

### KEGG & GO enrichment analysis
### loop to catch the results
lmfreq_go_results <- list()
lmfreq_kegg_results <- list()

for ( xzy in 1:ncol(summary(decideTests(efit))) ) { 
  
  lmfreq_go_results[[xzy]] <- goana(efit, coef=xzy, geneid = "entrez", species="Hs")
  lmfreq_go_results[[xzy]]$Comparison <- colnames(summary(decideTests(efit)))[xzy]
  lmfreq_go_results[[xzy]]$Marker <- rownames( lmfreq_go_results[[xzy]] )
  rownames(lmfreq_go_results[[xzy]]) <- NULL
  
  lmfreq_kegg_results[[xzy]] <- kegga(efit, coef=xzy, geneid = "entrez")
  lmfreq_kegg_results[[xzy]]$Comparison <- colnames(summary(decideTests(efit)))[xzy]
  lmfreq_kegg_results[[xzy]]$Marker <- rownames( lmfreq_kegg_results[[xzy]] ) 
  rownames(lmfreq_kegg_results[[xzy]]) <- NULL
}      

lmfreq_go_results <- do.call(rbind,lmfreq_go_results)
lmfreq_kegg_results <- do.call(rbind,lmfreq_kegg_results)

results_analysis$de_analysis_results <- de_analysis_results
results_analysis$lmfreq_go_results <- lmfreq_go_results
results_analysis$lmfreq_kegg_results <- lmfreq_kegg_results
### clean up
rm(vfit,efit,lmfreq_results,design,pca1,pca2,pca3,xzy,keep.counts,mydata,contr.matrix,lmfreq_kegg_results,lmfreq_go_results,de_analysis_results)
###
#############################################################################################################################################


#############################################################################################################################################
### Volcanos
#############################################################################################################################################

my_comp <- unique(results_analysis$de_analysis_results$Comparison)
#colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_colors <- list(DZ50='#D55E00',Un='#56B4E9',TGFB_DZ50='#009E73')

### comp1
ix <- which(results_analysis$de_analysis_results$Comparison %in% my_comp[1])

a<-results_analysis$de_analysis_results[ix,]
ix <- c(which(a$logFC < -0.5 & a$adj.P.Val < 0.05)[1:10],
        which(a$logFC > 0.5 & a$adj.P.Val < 0.05)[1:10])

a$colors <- 'NA'
a$colors[a$logFC>0 & a$adj.P.Val<0.05] <- 'DZ50'
a$colors[a$logFC<0 & a$adj.P.Val<0.05] <- 'Untreated'
a$logFC[a$logFC>2] <- 2
a$logFC[a$logFC < -2] <- -2

pdf(file="figures/volcano_DZ50_vs_Un_TD06527.pdf",width = 3,height = 3)
ggplot(a)+aes(y=-log10(adj.P.Val),x=logFC,color=colors)+geom_point(show.legend = F,size=0.5)+theme_classic()+
  geom_hline(yintercept = 1.3,color='deeppink',linetype=2)+
  geom_vline(xintercept = 0,color='deeppink',linetype=2)+
  scale_color_manual(values = c('#D55E00','black','#56B4E9'))+
  geom_text_repel(data=a[ix,],aes(label = gene_name), max.overlaps = 10000, size = 3, force = 3,label.size = 3,color='black')+ 
  labs(x='Log2FC',y='-Log10(FDR)')
dev.off()

### comp2
ix <- which(results_analysis$de_analysis_results$Comparison %in% my_comp[2])

my_comp[2]

a<-results_analysis$de_analysis_results[ix,]
ix <- c(which(a$logFC < -0.5 & a$adj.P.Val < 0.05)[1:10],
        which(a$logFC > 0.5 & a$adj.P.Val < 0.05)[1:10])

a$colors <- 'NA'
a$colors[a$logFC>0 & a$adj.P.Val<0.05] <- 'TGFB_DZ50'
a$colors[a$logFC<0 & a$adj.P.Val<0.05] <- 'Untreated'
a$logFC[a$logFC>2] <- 2
a$logFC[a$logFC < -2] <- -2

pdf(file="figures/volcano_TGFB_DZ50_vs_Un_TD06527.pdf",width = 3,height = 3)
ggplot(a)+aes(y=-log10(adj.P.Val),x=logFC,color=colors)+geom_point(show.legend = F,size=0.5)+theme_classic()+
  geom_hline(yintercept = 1.3,color='deeppink',linetype=2)+
  geom_vline(xintercept = 0,color='deeppink',linetype=2)+
  scale_color_manual(values = c('black','#009E73','#56B4E9'))+
  geom_text_repel(data=a[ix,],aes(label = gene_name), max.overlaps = 10000, size = 3, force = 3,label.size = 3,color='black')+ 
  labs(x='Log2FC',y='-Log10(FDR)')
dev.off()

### comp3
ix <- which(results_analysis$de_analysis_results$Comparison %in% my_comp[3])

my_comp[3]

a<-results_analysis$de_analysis_results[ix,]
ix <- c(which(a$logFC < -0.5 & a$adj.P.Val < 0.05)[1:10],
        which(a$logFC > 0.5 & a$adj.P.Val < 0.05)[1:10])

a$colors <- 'NA'
a$colors[a$logFC>0 & a$adj.P.Val<0.05] <- 'DZ50'
a$colors[a$logFC<0 & a$adj.P.Val<0.05] <- 'TGFB_DZ50'
a$logFC[a$logFC>2] <- 2
a$logFC[a$logFC < -2] <- -2

pdf(file="figures/volcano_DZ50_vs_TGFB_DZ50_TD06527.pdf",width = 3,height = 3)
ggplot(a)+aes(y=-log10(adj.P.Val),x=logFC,color=colors)+geom_point(show.legend = F,size=0.5)+theme_classic()+
  geom_hline(yintercept = 1.3,color='deeppink',linetype=2)+
  geom_vline(xintercept = 0,color='deeppink',linetype=2)+
  scale_color_manual(values = c('#009E73','black','#D55E00'))+
  geom_text_repel(data=a[ix,],aes(label = gene_name), max.overlaps = 10000, size = 3, force = 3,label.size = 3,color='black')+ 
  labs(x='Log2FC',y='-Log10(FDR)')
dev.off()
#############################################################################################################################################


#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image('RData/analysis_prerna.RData')
load('RData/analysis_prerna.RData')
#############################################################################################################################################



### PCA + heatmap and other figures
#############################################################################################################################################
ix <- which(!results_analysis$voom_TD06527$samples$condition %in% "Undetermined")

mydata <- prcomp(t(results_analysis$voom_TD06527$E[,ix]), scale=TRUE,center = TRUE) 
pca_mts<-as.data.frame(mydata$x)
pca_mts$sample_id <- rownames(pca_mts)
pca1 <- format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2)
pca3 <- format(round(get_eigenvalue(mydata)$variance.percent[3], 2), nsmall = 2)
pca_mts$treatment <- results_analysis$voom_TD06527$samples$condition[ix]

pdf(file="figures/pca_rnaseq_screeplot_samples_td06527.pdf",width = 3,height = 2)
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_classic() + ylim(0,50)
dev.off()

pdf(file="figures/pca_rnaseq_profiles_samples_td06527.pdf",width = 4,height = 3.5)
plot_grid(
  ggplot(pca_mts, aes(PC2, PC1, color = treatment)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=c("orangered","olivedrab3","#0072B2")) +
    labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Type / Origin") 
  ,
  cowplot::get_legend(ggplot(pca_mts, aes(PC3, PC1, color = treatment)) + geom_point() + theme_classic() + 
                        scale_color_manual(values=c("orangered","olivedrab3","#0072B2")) + guides(guide_legend(ncol=2)) + 
                        theme(legend.position = 'right') +
                        labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Type / Origin") )
  ,
  ggplot(pca_mts, aes(PC3, PC1, color = treatment)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=c("orangered","olivedrab3","#0072B2")) +
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Type / Origin") 
  ,
  ggplot(pca_mts, aes(PC3, PC2, color = treatment)) + geom_point(show.legend = FALSE) + theme_classic() + 
    scale_color_manual(values=c("orangered","olivedrab3","#0072B2")) +
    labs(x=paste("PC3:",pca3,"%",sep=""), y=paste("PC2:", pca2,"%", sep=""), color="Type / Origin") 
  ,
  labels = c("","","","","","","","","") 
)
dev.off()

pdf(file="figures/pca_rnaseq_profiles_samples_td06527_simple.pdf",width = 2.25,height = 2.50)
ggplot(pca_mts, aes(PC2, PC1, color = treatment)) + geom_point(show.legend = TRUE) + theme_classic() + 
  scale_color_manual(values=c("orangered","olivedrab3","#0072B2")) +
  labs(x=paste("PC2:",pca2,"%",sep=""), y=paste("PC1:", pca1,"%", sep=""), color="Type / Origin") +
  theme(legend.position="bottom")
dev.off()
#############################################################################################################################################


#############################################################################################################################################
### summary barplot
#############################################################################################################################################
table(results_analysis$de_analysis_results$Comparison[results_analysis$de_analysis_results$adj.P.Val<0.05])
#############################################################################################################################################
stacked_bar_plot <- results_analysis$de_analysis_results
stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05
stacked_bar_plot <- data.frame(table(stacked_bar_plot$Comparison, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
colnames(stacked_bar_plot) <- c('Comparison', 'logFC_dir','sig', 'count')
regional_dmr_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down'] <- 0 - regional_dmr_summary$count[regional_dmr_summary$logFC_dir=='Down']
regional_dmr_summary$Comparison <- as.character(regional_dmr_summary$Comparison) 
#############################################################################################################################################
pdf(file="figures/stack_deg_signals.pdf",width = 4,height = 2)
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
#############################################################################################################################################
library(Vennerable)
a <- results_analysis$de_analysis_results
y <- list( `TGFB+DZ50 vs Un` = a$Marker[which(a$Comparison %in% "TGFB_DZ50_vs_Un" & a$adj.P.Val<0.05)],
           `DZ50 vs TGFB+DZ50` = a$Marker[which(a$Comparison %in% "DZ50_vs_TGFB_DZ50" & a$adj.P.Val<0.05)],
           `DZ50 vs Un` = a$Marker[which(a$Comparison %in% "DZ50_vs_Un" & a$adj.P.Val<0.05)])
x_RNAS <- unique(unlist(y))
pdf(file="figures/venn_plot.pdf",width = 4,height = 4)
plot(Venn(y))
dev.off()
#############################################################################################################################################
### make comb
#list_to_matrix(lt)
m1 = make_comb_mat(y)
comb_name(m1)
set_size(m1)
extract_comb(m1, "100")
#############################################################################################################################################


#############################################################################################################################################
### summary table of DEGs
#############################################################################################################################################
results_analysis$voom_TD06527$results_summary_table <- as.data.frame(results_analysis$voom_TD06527$E)
results_analysis$voom_TD06527$results_summary_table$Geneid <- rownames(results_analysis$voom_TD06527$results_summary_table)
rownames(results_analysis$voom_TD06527$results_summary_table) <- NULL
ix <- which(colnames(results_analysis$voom_TD06527$genes) %in% c('Geneid','entrez','gene_name','gene_biotype'))

results_analysis$voom_TD06527$results_summary_table <- merge(results_analysis$voom_TD06527$results_summary_table,results_analysis$voom_TD06527$genes[,ix],by='Geneid')

# 'DZ50_vs_TGFB_DZ50
ix <- which(results_analysis$de_analysis_results$Comparison %in% 'DZ50_vs_TGFB_DZ50')
iy <- which(colnames(results_analysis$de_analysis_results) %in% c('Geneid','logFC','AveExpr','P.Value','adj.P.Val'))
a <- results_analysis$de_analysis_results[ix,iy]
colnames(a) <- c("Geneid","logFC_DZ50_vs_TGFB_DZ50","AveExpr_DZ50_vs_TGFB_DZ50","P.Value_DZ50_vs_TGFB_DZ50","adj.P.Val_DZ50_vs_TGFB_DZ50")
dim(results_analysis$voom_TD06527$results_summary_table)
results_analysis$voom_TD06527$results_summary_table <- merge(results_analysis$voom_TD06527$results_summary_table,a,by='Geneid')
dim(results_analysis$voom_TD06527$results_summary_table)

# 'DZ50_vs_Un
ix <- which(results_analysis$de_analysis_results$Comparison %in% 'DZ50_vs_Un')
iy <- which(colnames(results_analysis$de_analysis_results) %in% c('Geneid','logFC','AveExpr','P.Value','adj.P.Val'))
a <- results_analysis$de_analysis_results[ix,iy]
colnames(a) <- c("Geneid","logFC_DZ50_vs_Un","AveExpr_DZ50_vs_Un","P.Value_DZ50_vs_Un","adj.P.Val_DZ50_vs_Un")
dim(results_analysis$voom_TD06527$results_summary_table)
results_analysis$voom_TD06527$results_summary_table <- merge(results_analysis$voom_TD06527$results_summary_table,a,by='Geneid')
dim(results_analysis$voom_TD06527$results_summary_table)

# 'TGFB_DZ50_vs_Un
ix <- which(results_analysis$de_analysis_results$Comparison %in% 'TGFB_DZ50_vs_Un')
iy <- which(colnames(results_analysis$de_analysis_results) %in% c('Geneid','logFC','AveExpr','P.Value','adj.P.Val'))
a <- results_analysis$de_analysis_results[ix,iy]
colnames(a) <- c("Geneid","logFC_TGFB_DZ50_vs_Un","AveExpr_TGFB_DZ50_vs_Un","P.Value_TGFB_DZ50_vs_Un","adj.P.Val_TGFB_DZ50_vs_Un")
dim(results_analysis$voom_TD06527$results_summary_table)
results_analysis$voom_TD06527$results_summary_table <- merge(results_analysis$voom_TD06527$results_summary_table,a,by='Geneid')
dim(results_analysis$voom_TD06527$results_summary_table)

# signatures
results_analysis$voom_TD06527$results_summary_table$unique_to_TGFB_DZ50_vs_Un <- results_analysis$voom_TD06527$results_summary_table$Geneid %in% extract_comb(m1, "100")
results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 <- results_analysis$voom_TD06527$results_summary_table$Geneid %in% extract_comb(m1, "010")
results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un <- results_analysis$voom_TD06527$results_summary_table$Geneid %in% extract_comb(m1, "001")

### review
results_analysis$voom_TD06527$results_summary_table$adj.P.Val_DZ50_vs_Un[which(results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un %in% 'TRUE')]

### 
results_analysis$voom_TD06527$results_summary_table$signatures <- 'NA'
results_analysis$voom_TD06527$results_summary_table$signatures[results_analysis$voom_TD06527$results_summary_table$unique_to_TGFB_DZ50_vs_Un %in% 'TRUE'] <- "TGFB_DZ50_vs_Un"
results_analysis$voom_TD06527$results_summary_table$signatures[results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 %in% 'TRUE'] <- "DZ50_vs_TGFB_DZ50"
results_analysis$voom_TD06527$results_summary_table$signatures[results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un %in% 'TRUE'] <- "DZ50_vs_Un"
###
write.csv(file='tables/results_summary_table_full.csv',results_analysis$voom_TD06527$results_summary_table,row.names=FALSE)
###
ix <- which(results_analysis$voom_TD06527$results_summary_table$unique_to_TGFB_DZ50_vs_Un %in% 'TRUE' | 
        results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 %in% 'TRUE' | 
        results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un %in% 'TRUE')
###
write.csv(file='tables/results_summary_table_simple.csv',results_analysis$voom_TD06527$results_summary_table[ix,],row.names=FALSE)
### clear
rm(a,mydata,pca_mts,regional_dmr_summary,stacked_bar_plot)
#############################################################################################################################################


### Expression of uniquely significant genes
#############################################################################################################################################
## adjust
xx <- which(results_analysis$voom_TD06527$results_summary_table$unique_to_TGFB_DZ50_vs_Un %in% 'TRUE' | # & results_analysis$voom_TD06527$results_summary_table$logFC_TGFB_DZ50_vs_Un > 0  
                    results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 %in% 'TRUE' | 
                    results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un %in% 'TRUE' ) # & results_analysis$voom_TD06527$results_summary_table$logFC_DZ50_vs_Un > 0

generic_expression_matrix <- results_analysis$voom_TD06527$results_summary_table[xx,2:11]
generic_expression_matrix <- generic_expression_matrix[,-3]

results_analysis$voom_TD06527$results_summary_table$gene_nameV2 <- paste(1:nrow(results_analysis$voom_TD06527$results_summary_table),results_analysis$voom_TD06527$results_summary_table$gene_name,sep='_')

rownames(generic_expression_matrix) <- results_analysis$voom_TD06527$results_summary_table$gene_nameV2[xx]

### ann
annotation_row <- rowAnnotation( Signature = results_analysis$voom_TD06527$results_summary_table$signatures[xx],
                                 col = list(Signature=c(TGFB_DZ50_vs_Un='deeppink1',
                                                        DZ50_vs_Un='skyblue1',
                                                        DZ50_vs_TGFB_DZ50='lightgreen')))
                                 #Biotype = results_analysis$voom_TD06527$results_summary_table$gene_biotype[xx]

genericVar <-  t(scale(t(generic_expression_matrix)))

### colum
annotation_column = HeatmapAnnotation( 
  Treatment = factor(results_analysis$voom_TD06527$samples$condition[-3],levels = c('Untreated','TGFB_DZ50','DZ50')),
  col = list(Treatment=c(TGFB_DZ50='olivedrab3',
                         DZ50='orangered',
                         Untreated='#0072B2')
  ) )

#genericVar[genericVar > 3.5] <- 3.5
#genericVar[genericVar < -3.5] <- -3.5

pdf(file='figures/expression_profiles_log2cpm_unique_signatures.pdf', width = 7, height = 5)
draw(
  Heatmap( genericVar, 
           col = colorRampPalette(c("steelblue", "white","firebrick"))(255),
           show_row_names = FALSE , 
           show_column_names = TRUE, 
           cluster_columns = TRUE,
           cluster_rows =  TRUE , 
           #column_split = df_col,
           row_names_gp = gpar(fontsize = 5),
           #column_km = 3,
           cluster_column_slices = FALSE,
           column_names_rot = 45,
           clustering_method_rows = 'ward.D2',
           clustering_method_columns = 'ward.D2',
           top_annotation = annotation_column,
           right_annotation = annotation_row,
           border = TRUE,
           use_raster = FALSE,
           name = "Z-Score")
  ,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()
#############################################################################################################################################

### Expression of top 15 
#############################################################################################################################################
## adjust
xx <- c(which(results_analysis$voom_TD06527$results_summary_table$unique_to_TGFB_DZ50_vs_Un %in% 'TRUE' & 
              results_analysis$voom_TD06527$results_summary_table$logFC_TGFB_DZ50_vs_Un > 0.5 &
              results_analysis$voom_TD06527$results_summary_table$adj.P.Val_TGFB_DZ50_vs_Un < 0.005 &
                results_analysis$voom_TD06527$results_summary_table$gene_biotype %in% 'protein_coding')
,
which(results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_Un %in% 'TRUE' & 
        results_analysis$voom_TD06527$results_summary_table$logFC_DZ50_vs_Un > 0.30 &
        results_analysis$voom_TD06527$results_summary_table$adj.P.Val_DZ50_vs_Un < 0.05 &
        results_analysis$voom_TD06527$results_summary_table$gene_biotype %in% 'protein_coding')
,
which(results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 %in% 'TRUE' & 
        results_analysis$voom_TD06527$results_summary_table$logFC_DZ50_vs_TGFB_DZ50 > 0.30 &
        results_analysis$voom_TD06527$results_summary_table$adj.P.Val_DZ50_vs_TGFB_DZ50 < 0.05 &
        results_analysis$voom_TD06527$results_summary_table$gene_biotype %in% 'protein_coding')
,
which(results_analysis$voom_TD06527$results_summary_table$unique_to_DZ50_vs_TGFB_DZ50 %in% 'TRUE' & 
        results_analysis$voom_TD06527$results_summary_table$logFC_DZ50_vs_TGFB_DZ50 < -0.35 &
        results_analysis$voom_TD06527$results_summary_table$adj.P.Val_DZ50_vs_TGFB_DZ50 < 0.05  &
        results_analysis$voom_TD06527$results_summary_table$gene_biotype %in% 'protein_coding') )

xx <- unique(xx)

length(unique(xx))

generic_expression_matrix <- results_analysis$voom_TD06527$results_summary_table[xx,2:11]
generic_expression_matrix <- generic_expression_matrix[,-3]

results_analysis$voom_TD06527$results_summary_table$gene_nameV2 <- paste(1:nrow(results_analysis$voom_TD06527$results_summary_table),results_analysis$voom_TD06527$results_summary_table$gene_name,sep='_')

rownames(generic_expression_matrix) <- results_analysis$voom_TD06527$results_summary_table$gene_name[xx]

### ann
annotation_row <- rowAnnotation( Signature = results_analysis$voom_TD06527$results_summary_table$signatures[xx],
                                 col = list(Signature=c(TGFB_DZ50_vs_Un='deeppink1',
                                                        DZ50_vs_Un='skyblue1',
                                                        DZ50_vs_TGFB_DZ50='lightgreen')))
#Biotype = results_analysis$voom_TD06527$results_summary_table$gene_biotype[xx]

genericVar <-  t(scale(t(generic_expression_matrix)))

### colum
annotation_column = HeatmapAnnotation( 
  Treatment = factor(results_analysis$voom_TD06527$samples$condition[-3],levels = c('Untreated','TGFB_DZ50','DZ50')),
  col = list(Treatment=c(TGFB_DZ50='olivedrab3',
                         DZ50='orangered',
                         Untreated='#0072B2')
  ) )

#genericVar[genericVar > 3.5] <- 3.5
#genericVar[genericVar < -3.5] <- -3.5

pdf(file='figures/expression_profiles_log2cpm_top_unique_signatures.pdf', width = 6, height = 6)
draw(
  Heatmap( genericVar, 
           col = colorRampPalette(c("steelblue", "white","firebrick"))(255),
           show_row_names = TRUE , 
           show_column_names = TRUE, 
           cluster_columns = TRUE,
           cluster_rows =  TRUE , 
           #column_split = df_col,
           row_names_gp = gpar(fontsize = 5),
           #column_km = 3,
           cluster_column_slices = FALSE,
           column_names_rot = 45,
           clustering_method_rows = 'ward.D2',
           clustering_method_columns = 'ward.D2',
           top_annotation = annotation_column,
           #right_annotation = annotation_row,
           border = TRUE,
           use_raster = FALSE,
           name = "Z-Score")
  ,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()
#############################################################################################################################################


#############################################################################################################################################
### NEW pathway figure
#############################################################################################################################################
library(rrvgo)

### select markers
ix <- which(results_analysis$lmfreq_go_results$Comparison %in% 'DZ50_vs_Un' & 
              results_analysis$lmfreq_go_results$P.Down > 0.05 & 
              results_analysis$lmfreq_go_results$P.Up < 0.05)

simMatrix <- calculateSimMatrix(results_analysis$lmfreq_go_results$Marker[ix],
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10((results_analysis$lmfreq_go_results$P.Up[ix])), results_analysis$lmfreq_go_results$Marker[ix])

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

pdf(file='figures/rrvgo_DZ50_vs_Un.pdf', width = 7, height = 5)
treemapPlot(reducedTerms)
dev.off()

### select markers
ix <- which(results_analysis$lmfreq_go_results$Comparison %in% 'TGFB_DZ50_vs_Un' & 
              results_analysis$lmfreq_go_results$P.Down > 0.05 & 
              results_analysis$lmfreq_go_results$P.Up < 0.05)

simMatrix <- calculateSimMatrix(results_analysis$lmfreq_go_results$Marker[ix],
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10((results_analysis$lmfreq_go_results$P.Up[ix])), results_analysis$lmfreq_go_results$Marker[ix])

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

pdf(file='figures/rrvgo_TGFB_DZ50_vs_Un.pdf', width = 7, height = 5)
treemapPlot(reducedTerms)
dev.off()

### select markers
ix <- which(results_analysis$lmfreq_go_results$Comparison %in% 'DZ50_vs_TGFB_DZ50' & 
              results_analysis$lmfreq_go_results$P.Down > 0.05 & 
              results_analysis$lmfreq_go_results$P.Up < 0.05)

simMatrix <- calculateSimMatrix(results_analysis$lmfreq_go_results$Marker[ix],
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10((results_analysis$lmfreq_go_results$P.Up[ix])), results_analysis$lmfreq_go_results$Marker[ix])

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

pdf(file='figures/rrvgo_DZ50_vs_TGFB_DZ50_higher_in_DZ50.pdf', width = 7, height = 5)
treemapPlot(reducedTerms)
dev.off()

### select markers
ix <- which(results_analysis$lmfreq_go_results$Comparison %in% 'DZ50_vs_TGFB_DZ50' & 
              results_analysis$lmfreq_go_results$P.Down < 0.05 & 
              results_analysis$lmfreq_go_results$P.Up > 0.05)

simMatrix <- calculateSimMatrix(results_analysis$lmfreq_go_results$Marker[ix],
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10((results_analysis$lmfreq_go_results$P.Up[ix])), results_analysis$lmfreq_go_results$Marker[ix])

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

pdf(file='figures/rrvgo_DZ50_vs_TGFB_DZ50_higher_in_TGFB.pdf', width = 7, height = 5)
treemapPlot(reducedTerms)
dev.off()
#############################################################################################################################################



#############################################################################################################################################
### checkpoint
#############################################################################################################################################
#save.image('RData/analysis_prerna.RData')
#load('RData/analysis_prerna.RData')
#############################################################################################################################################