#Title: Differential expression analysis for Mn-induced neurotoxicity
# Description: Takes in a counts matrix or raw outputputs from STAR aligner
#               . sample metadata, runs DESeq2 and outputs results + volcanoPlot
# Author: George W Kagugube
# Date: 2023 - 2024
# Usage:
#       currently run within rstudio or R
###############################################################################
## Clear the environmental space here to start afresh 
## Clear environmental space
rm(list = ls())

## Load the helper functions for this analysis
suppressPackageStartupMessages({
  library(optparse)
  source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")
})

options(
  stringAsFactors = FALSE,
  scipen = 999
)

## Set seed for reproducibility here
set.seed(101)

###############################################################################
## Setup command line arguments (not needed to run this analysis but can be integrated)
options_list <- list(
  make_option(c('-c', '--counts', type = 'character',help = 'Counts file (tsv)')),
  make_option(c("-m", "--meta"),    type = "character", help = "Sample metadata (TSV)"),
  make_option(c("-o", "--outdir"),  type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$meta) || is.null(opt$outdir)) {
  stop("Missing required arguments. Run with --help for usage.", call. = FALSE)
}

###############################################################################
## Set the directory with the expression datasets here 
setwd('/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/data')

# This is the output folder for the final analysis
#output_dir = '/Users/gwk/Desktop/Thesis Figures/DifferentialGeneExpression/'
output_dir = '/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/data'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

## Load the datasets to be analysed here 
countMatrix <- read.csv('raw_count_matrix.csv', row.names = 1)
samples <- read.csv("./samplieinformation.csv", row.names = 1)
###############################################################################
## Split the datasets into WT and Home
wt_sample_info <- samples |>
  filter(Genotype == 'wt')
wtcountmtx <- countMatrix |>
  select(row.names(wt_sample_info))

mut_sample_info <- samples |>
  filter(Genotype == 'mutant')
mutcountmtx <- countMatrix |>
  select(row.names(mut_sample_info))

## Filter before creating the DESeq object
##Although filtering with DESeq is not reccomended, here we remove all genes whose
## row sum is less than 10
keep <- rowSums(countMatrix > 10) >= 6
countMatrix <- countMatrix[keep,]

keep_wt <- rowSums(wtcountmtx > 10) >= 6
wtcountmtx <- wtcountmtx[keep_wt,]

keep_mut <- rowSums(mutcountmtx > 10) >= 6
mutcountmtx <- mutcountmtx[keep_mut,]
###############################################################################
# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Group)

## Wild type only datasets, see PCA sample clustering for more insight into this choice
ddswt <- DESeqDataSetFromMatrix(countData = as.matrix(wtcountmtx),
                                colData = wt_sample_info,
                                design = ~ Group)

## Mutant only datasets, see PCA sample clustering for more insight into this choice
ddsmut <- DESeqDataSetFromMatrix(countData = as.matrix(mutcountmtx),
                                colData = mut_sample_info,
                                design = ~ Group)
###############################################################################
## Relevel the conditions here
dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
#dds$Group <- relevel(dds$Group, ref = 'mut_unexposed')

# Relevel the Wild type data here so the unexposed is the reference group
ddswt$Group <- relevel(ddswt$Group, ref = 'wt_unexposed')
ddsmut$Group <- relevel(ddsmut$Group, ref = 'mut_unexposed')
###############################################################################
## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
# ## check out the individual functions below in the environmnetal space
# ## calculate the pca values here
vsd <- vst(ddswt)

## Plot visualise the PCA
plotPCA(vsd, intgroup="Group")
# 
# ## Normalise the count data
# all_norm_data <- normalisation_func(dds)
# wt_norm_data <- normalisation_func(ddswt)
# mut_norm_data <- normalisation_func(ddsmut)
# 
# ## calculate the principle component analysis
# pca_all <- principle_component(dds)
#   
# # Create a data frame that can used going forward from here on 
# count_all <- cbind(samples, pca_all$x)
# 
# ## Visualise the data here 
# ggplot(data = count_all) +
#   geom_point(aes(x=PC1, y=PC2, colour = Group), size=5) +
#   theme_minimal() +
#   labs(x = 'PC1: 29% variance',
#        y = 'PC2: 23% varience') +
#   theme(
#     axis.text = element_text(size = 20),
#     axis.text.x = element_text(angle = 90, vjust = 0.5),
#     axis.title.x = element_text(size = 15, vjust = 0.5),
#     axis.text.y = element_text(angle = 90, vjust = 0.5),
#     axis.title.y = element_text(size = 15, vjust = 0.5)
#   )
###############################################################################  
# ## Perform differenPC10## Perform differential gene expression analysis here
dds <- DESeq(dds)
ddsmut <- DESeq(ddsmut)
ddswt <- DESeq(ddswt)
###############################################################################

mutantUnexposed |>
  as.data.frame() |>
  filter(padj < 0.05) |>
  arrange(desc(log2FoldChange)) |>
  nrow()
## Extract Group comparisons here
## ========================= Mutant,Exposed vs Mutant unexposed ===============
resMut <- results(ddsmut)
summary(resMut)
resultsNames(ddsmut)
mutantExposed <- results(ddsmut, name = "Group_mut_exposed_vs_mut_unexposed")
mutantExposed <- lfcShrink(ddsmut, coef = 2, type = 'apeglm')
mutantExposed <- addDirectionlabel(mutantExposed)
mutantExposed <- annot_data(mutantExposed)
head(mutantExposed)
# Remove all genes that are not annotated here
mutantExposed <- mutantExposed[!grepl('LOC', rownames(mutantExposed)), ]
mutantExposed <- mutantExposed[!grepl('si:', rownames(mutantExposed)), ]
mutantExposed <- mutantExposed[!grepl('zgc:', rownames(mutantExposed)), ]
mutantExposed <- mutantExposed[!grepl('wu:', rownames(mutantExposed)), ]
# visualise the data using a volcano plot here
volcanoPlot(mutantExposed, xlimlimit = c(-2.5,3.5),
            ylimlimit = c(0,15))

## Export the the DGE here
write.csv(mutantExposed,
          file = paste0(output_dir,
                        'CleanAnalysis_without_noise/mutantExposed/mutExposed_dge.csv'))
###############################################################################
## =======================Wild type Exposed vs Unexposed =======================
resWT <- results(ddswt)
resWT
summary(resWT)
resultsNames(ddswt)
WTExposed <- results(ddswt, name = "Group_wt_exposed_vs_wt_unexposed")
WTExposed <- lfcShrink(ddswt, coef = 2, type = 'apeglm')
WTExposed <- addDirectionlabel(WTExposed)
WTExposed <- annot_data(WTExposed)
summary(WTExposed)
# Remove all genes that are not annotated here
WTExposed <- WTExposed[!grepl('LOC', rownames(WTExposed)), ]
WTExposed <- WTExposed[!grepl('si:', rownames(WTExposed)), ]
WTExposed <- WTExposed[!grepl('zgc:', rownames(WTExposed)), ]
WTExposed <- WTExposed[!grepl('wu:', rownames(WTExposed)), ]
# visualise the data using a volcano plot here
volcanoPlot(WTExposed,
            xlimlimit = c(-1, 5))

## Export the the DGE here
write.csv(WTExposed,
          file = paste0(output_dir,
                        'CleanAnalysis_without_noise/wtExposed/wtExposed_dge.csv'))
###############################################################################
## ============= Combined data analysis for interactions and mutant effect =====
res <- results(dds)
res
summary(res)
resultsNames(dds)

## ================================= Mutant only effects =======================
mutantUnexposed <- results(dds, name = "Group_mut_unexposed_vs_wt_unexposed")
mutantUnexposed <- lfcShrink(dds = dds, coef = 3, type = 'apeglm')
mutantUnexposed <- addDirectionlabel(mutantUnexposed)
mutantUnexposed <- annot_data(mutantUnexposed)
head(mutantUnexposed)
# Remove all genes that are not annotated here
mutantUnexposed <- mutantUnexposed[!grepl('LOC', rownames(mutantUnexposed)), ]
mutantUnexposed <- mutantUnexposed[!grepl('si:', rownames(mutantUnexposed)), ]
mutantUnexposed <- mutantUnexposed[!grepl('zgc:', rownames(mutantUnexposed)), ]
mutantUnexposed <- mutantUnexposed[!grepl('wu:', rownames(mutantUnexposed)), ]

## Remove all genes that are not annotated in the reference genome
volcanoPlot(mutantUnexposed,
            xlimlimit = c(-5,8),
            ylimlimit = c(0, 20))

## Export the datset here
write.csv(mutantUnexposed,
          file = paste0(output_dir,
                        'CleanAnalysis_without_noise/MutantUnexposed/mutUnexposed_dge.csv'))

###############################################################################

