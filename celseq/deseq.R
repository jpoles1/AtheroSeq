library(DESeq2)
library(ggplot2)
library(biomaRt)
library(goseq)

species.db <- "mmusculus_gene_ensembl"
db.version <- "grcm38"
contrasts <- list( c( "Experimental.Group", "PA", "RA" ),
                   c( "Experimental.Group", "PA", "HA" ),
                   c( "Experimental.Group", "RA", "HA" ),
                   c( "Experimental.Group", "PA", "BA" ),
                   c( "Experimental.Group", "RA", "BA" ))
fdr <- 0.01
lfc <- 0

# read in metadata
metadata <- read.csv("data/rmetadata.csv")
# generate object for DESeq2
study_data <- read.csv("data/rcount_matrix.csv")
rownames(study_data) = study_data$Sample.ID
rownames(metadata) = rownames(study_data)
study_data$Sample.ID = NULL
cds <- DESeqDataSetFromMatrix(countData = t(study_data),
                              colData = metadata,
                              design = ~ Experimental.Group)
cds <- cds[ which( rowSums( counts( cds ) ) > 1 ), ]
counts <- counts( cds )

# NBM fit
dds <- DESeq( cds )
saveRDS( dds, file = "output/deseq2_fit.RDS" )
plotDispEsts( dds, main = "Dispersion plot" )

# create database for annotation before running loops
mart <- useMart( "ENSEMBL_MART_ENSEMBL", dataset = species.db ) 
# host = "www.biomart.org" )
# host = paste0( db.version, ".ensembl.org" ) )
filters <- c("external_gene_name", "mgi_symbol")
attributes <- c("external_gene_name", "mgi_symbol",  "name_1006", "description")

# loop through all the contrasts 
for (fac in 1:length(contrasts)) {
  currcon <- contrasts[[fac]]
  res <- results( dds, contrast = currcon )
  
  pdf( paste0("output/", currcon[2], "_vs_", currcon[3], ".pdf" ) )
  DESeq2::plotMA( res, main = paste0( currcon[2], " vs ", currcon[3] ),
                  ylim = c( -5, 5 ), alpha = fdr )
  dev.off()
  
  # All differential genes
  all <- subset( res, padj < fdr )
  all <- all[ order( all$log2FoldChange, decreasing = T ), ]
  # all annotation
  all.genes <- rownames( all ) # rownames = ensembl gene ids
  all.genes <- data.frame( all.genes )
  colnames( all.genes )[1] <- "Ensembl"
  rows = all.genes[,1]
  all.anno <- getBM( attributes = attributes, filters = filters, 
                     values = list("external_gene_name" = rows, "mgi_symbol"=rows),
                     mart = mart) # get all annotation based on ensembl ids
  all.anno = all.anno[all.anno$name_1006 != "",]
  all.anno = aggregate(name_1006~external_gene_name, all.anno, toString)
  all.anno$name_1006 = gsub(",", ";", all.anno$name_1006)
  # maintain un-annotated genes in BioMart, otherwise it will just get excluded in annotation data frame 
  # resulting in un-identified rows in the stats data frame
  all.anno2 <- merge( all.genes, all.anno, by.x = "Ensembl", by.y = "external_gene_name",
                      all.x = T )
  all.anno2 <- all.anno2[ match( all.genes[,1], all.anno2$Ensembl ), ] # reorder the annotation by LFC order
  all.anno.wStats <- data.frame( all.anno2, all )
  write.csv( all.anno.wStats, paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_ALL.csv" ),
               sep="\t", quote=F, row.names=F, col.names=T )
  write.csv( all.anno.wStats[,2], paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_ALL_GeneList.csv" ),
               sep="\t", quote=F, row.names=F, col.names=F )
  
  # Upregulated genes only
  up.anno.wStats <- subset(all.anno.wStats, padj < fdr & log2FoldChange > lfc )
  up.anno.wStats <- up.anno.wStats[ order( up.anno.wStats$log2FoldChange, decreasing = T ), ]
  
  write.csv( up.anno.wStats, paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_up.csv" ),
               sep="\t", quote=F, row.names=F, col.names=T )
  write.csv( up.anno.wStats[,2], paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_up_GeneList.csv" ),
               sep="\t", quote=F, row.names=F, col.names=F )
  
  # Downregulated genes only
  down.anno.wStats <- subset(all.anno.wStats, padj < fdr & log2FoldChange < lfc )
  down.anno.wStats <- down.anno.wStats[ order( down.anno.wStats$log2FoldChange, decreasing = F ), ]
  
  write.csv( down.anno.wStats, paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_down.csv" ),
               sep="\t", quote=F, row.names=F, col.names=T )
  write.csv( down.anno.wStats[,2], paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_down_GeneList.csv" ),
               sep="\t", quote=F, row.names=F, col.names=F )
  #Differential Gene Ontologies
  DE_genes <- as.integer(res$padj < fdr)
  names(DE_genes) <- rownames(res)
  diff_GO = goseq(nullp(DE_genes[complete.cases(DE_genes)], "mm9", "geneSymbol"), "mm9", "geneSymbol")
  diff_GO = diff_GO[diff_GO$over_represented_pvalue < fdr | diff_GO$under_represented_pvalue < fdr,]
  #all.go <- getgo(rows, "mm9", "geneSymbol",fetch.cats=c("GO:CC","GO:BP","GO:MF"))
  write.csv( diff_GO, paste0("output/", currcon[2], "_vs_", currcon[3], "_FDR", fdr, "_GO.csv" ),
             sep="\t", quote=F, row.names=F, col.names=T )
}
  