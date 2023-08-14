library(decontam)
library(phyloseq)
packageVersion("decontam") # 1.1.2 when this was put together
#colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case
asv_tab<-read.table("../../Microbiome_DADA2_Results/ASVs_counts.tsv", header=T, row.names=1,
                    check.names=F, sep="\t")
asv_tax <- as.matrix(read.table("../../Microbiome_DADA2_Results/ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
asv_fasta<-Biostrings::readDNAStringSet("../../Microbiome_DADA2_Results/ASVs.fa")
asv_seqs <- (as.vector(asv_fasta))
asv_headers <- names(asv_seqs)

for (i in 1:length(asv_headers)) {
  asv_headers[i] <- paste(">", asv_headers[i], sep="")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_tree<-read_tree_greengenes("../../Microbiome_DADA2_Results/ASVs-gg_placement.tog.tre")
tree_remove<-asv_tree$tip.label[!grepl("ASV",asv_tree$tip.label)]
asv_tree<-drop.tip(asv_tree,tree_remove)


map_file <- "../../Microbiome_DADA2_Results/metadata.txt" #Define pathway for mapping file
metdat<-as.data.frame(read.delim(map_file))
metdat$AgeClass<-ifelse(metdat$Species=="Antarctica",ifelse(grepl("J",metdat$Sample),"Juvuenile","Adult"),NA)
asv_tab<-asv_tab[names(asv_tab) %in% metdat$Sample]

library(phyloseq)

meta_reformed<-(metdat[match(colnames(asv_tab), metdat$Sample),])
rownames(meta_reformed)<-1:nrow(meta_reformed)
meta_reformed$Sample_or_Control<-ifelse(meta_reformed$Sample_or_Control=="Sample",FALSE,TRUE)
contam_df <- isContaminant(t(asv_tab), neg=meta_reformed$Sample_or_Control,threshold = 0.1)

table(contam_df$contaminant) # identified 6 as contaminants
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

