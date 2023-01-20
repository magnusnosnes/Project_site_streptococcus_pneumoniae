
```{r include=F}
ggplot(both_genes %>% filter(n_amino_acids == 328), aes(x = penicillin_G)) +
  geom_density() +
  ggtitle("Gene length versus penicillin levels") +
  xlab("Gene length (amino acids)") +theme_light()


#Info about the genes
table(na.omit(both_genes$gene_of_interest))
c1 = na.omit(both_genes)
length(unique(c1$Identity))
```

```{r, include=F, eval=F}

#   ____________________________________________________________________________
#   Exploratory section:                                                    ####

BiocManager::install("GenomicFeatures")
browseVignettes(("GenomicFeatures"))

test = GenomicFeatures::makeTxDbFromGFF(gff_files[1])

# We have some different options
# GenomicFeatures
# rtracklayer
# 


#load library
library(rtracklayer)
test = rtracklayer::readGFFAsGRanges(gff_files[1])
test2 = rtracklayer::readGFF(gff_files[1])
test3 = rtracklayer::readGFFPragmas(gff_files[1])
db = GenomicFeatures::makeTxDbFromGFF(gff_files[1])
# Create a SQLite database
db <- makeDbPackageFromGFF(gff_files[1], format = "gff3", type = "SQLite")
gene_info <- select(db, keys = "KhpB", columns = c("start", "end"), keytype = "gene")

# Extract the information for the gene of interest
gene_info <- select(db, keys = "gene_name", columns = c("start", "end"))

#Calculate the length of the gene
gene_length <- nchar(gene_info$end - gene_info$start)


# install.packages("devtools")
devtools::install_github("thackl/thacklr")
devtools::install_github("thackl/gggenomes")

library(gggenomes)
# Read features with it
test5 = gggenomes::read_gff3(gff_files[1])
results$`Alignment length:`








```

```{r , include=FALSE}

#   ____________________________________________________________________________
#   Old code                                                                ####


library(Biostrings)
sequences_of_interest <- readDNAStringSet("/Users/magnusnygardosnes/Desktop/Postdoc/Pneumo/Sequences_of_interest_NMBU/SOI.txt")

#Pseudocode of what we want to do:
# List all the isoaltes we have MIC data. 
# Blast search in these fasta files for the Elo_R deletion.
# Blast search in the fasta files for full Elo_R sequence.
# Compare the MIC of full Elo_R's and reduced Elo_R
# Look at the combination with variants of pbp2b. 

eloR_full = read.table("/Users/magnusnygardosnes/Desktop/Postdoc/Pneumo/Sequences_of_interest_NMBU/blast_output_eloR_full.txt")
colnames(eloR_full) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
eloR_reduced = read.table("/Users/magnusnygardosnes/Desktop/Postdoc/Pneumo/Sequences_of_interest_NMBU/blast_output_eloR_reduced.txt")
colnames(eloR_reduced)= c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(eloR_full)

#Look at full matches for the full eloR
full_matches_eloR_full = eloR_full[eloR_full$pident == 100,]
#Extract their sequence identifiers
eloR_full_sequence_identifiers = unique(substring(full_matches_eloR_full$sseqid, 1, as.integer(regexpr("_", full_matches_eloR_full$sseqid)[1])-1))

#eloR_reduced
full_matches_eloR_reduced = eloR_reduced[eloR_reduced$pident == 100,]
#Extract their sequence identifiers
eloR_reduced_sequence_identifiers = unique(substring(full_matches_eloR_reduced$sseqid, 1, as.integer(regexpr("_", full_matches_eloR_reduced$sseqid)[1])-1))


#Read MIC data. 
penicillin = read.table("/Users/magnusnygardosnes/Desktop/Postdoc/Pneumo/06_Pyseer/phenotype.tsv", header=T)
eloR_full_indexes = unlist(lapply(eloR_full_sequence_identifiers, FUN=function(x) which(x==penicillin$samples)))
eloR_reduced_indexes = unlist(lapply(eloR_reduced_sequence_identifiers, FUN=function(x) which(x==penicillin$samples)))
eloR_full_MIC = penicillin$penicillin_G[eloR_full_indexes]
eloR_reduced_MIC = penicillin$penicillin_G[eloR_reduced_indexes]

MIC = c(eloR_full_MIC,eloR_reduced_MIC)
gene_variants = c(rep("eloR_full", length(eloR_full_MIC)), rep("eloR_reduced", length(eloR_reduced_MIC)))
dat1 = cbind(as.numeric(MIC),gene_variants)
colnames(dat1)=c("MIC", "gene_variants")
dat1 = as.data.frame(dat1)

library(ggplot2)

ggplot(dat1, aes(x = gene_variants, y = as.numeric(MIC))) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_jitter(size = 2, alpha = 0.5, color = "red", ) +
  ggtitle("MIC by gene variants") +
  xlab("Gene variants") +
  ylab("MIC")+theme_light()+coord_flip()


ggplot(dat1, aes(x = gene_variants, y = as.numeric(MIC))) +
  geom_boxplot(fill = "#E6E6E6", color = "black") +
  geom_jitter(size = 2, alpha = 0.5, color = "#FF5722") +
  ggtitle("MIC by gene variants") +
  xlab("Gene variants") +
  ylab("MIC")+
  theme_light()+
  coord_flip()

colors <- c("#000000", "#1A1A1A", "#333333", "#4D4D4D", "#666666", "#808080", "#999999", "#B3B3B3", "#CCCCCC", "#E6E6E6", "#FFFFFF", 
                     "#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", 
                     "#CDDC39", "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E", "#607D8B", "#E91E63", "#9C27B0", "#3F51B5", 
                     "#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A",
                     "#CDDC39", "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#9E9E9E", "#607D8B", "#795548", "#9E9E9E", "#607D8B", "#E91E63",
                     "#9C27B0", "#3F51B5", "#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4", "#009688",
                     "#4CAF50", "#8BC34A", "#CDDC39", "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#9E9E9E", "#607D8B", "#795548", "#9E9E9E",
                     "#607D8B")
                     

annotation_first_eloR_reduced = read.delim("/Users/magnusnygardosnes/Desktop/Postdoc/Pneumo/pneumo_project_USA_ERR/Spades_assembly/Filtered/ERR876551.gff3", sep="\t", header=F, comment.char="#")
colnames(annotation_first_eloR_reduced) = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
eloR_reduced

# In other words: khpA or KhpA is eloR.

# Next steps: 

```
