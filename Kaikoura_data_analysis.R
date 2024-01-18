
setwd("~/Dropbox/PhD Project/Kaikoura_MS/Kaikoura_Analysis/")
library(tidyverse)
library(microbiome)
library(phyloseq)
count_tab<-read.table("../../Microbiome_DADA2_Results/ASVs_counts-no-contam.tsv", header=T, row.names=1,
                      check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("../../Microbiome_DADA2_Results/ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
tax_tab_lotus<-as.matrix(read.delim("../../Microbiome_DADA2_Results/ASVs-no-contam.fa.hier"))
tax_tab_lotus<-gsub("[?]",NA,tax_tab_lotus)
rownames(tax_tab_lotus)<-tax_tab_lotus[,1]
tax_tab_lotus<-tax_tab_lotus[,2:8]
colnames(tax_tab_lotus)<-colnames(tax_tab)
tax_tab<-tax_tab_lotus
tax_tab<-as.data.frame(tax_tab_lotus)
tax_tab$domain<-paste0("k__",tax_tab$domain)
tax_tab$phylum<-paste0("p__",tax_tab$phylum)
tax_tab$class<-paste0("c__",tax_tab$class)
tax_tab$order<-paste0("o__",tax_tab$order)
tax_tab$family<-paste0("f__",tax_tab$family)
tax_tab$genus<-paste0("g__",tax_tab$genus)
tax_tab$species<-paste0("s__",tax_tab$species)
tax_tab_taxa<-paste(tax_tab$domain,tax_tab$phylum,tax_tab$class,tax_tab$order,
                    tax_tab$family,tax_tab$genus,tax_tab$species,sep="; ")
tax_tab_taxa<-data.frame(tax_tab_taxa,rownames(tax_tab))
tax_tab2keep<-tax_tab_taxa$rownames.tax_tab.[!(grepl("Chloroplast|Mitochondria|Eukaryota",(tax_tab_taxa$tax_tab_taxa)))]
tax_tab_lotus_filt<-tax_tab_lotus[rownames(tax_tab_lotus) %in% tax_tab2keep,]

fastas<-Biostrings::readDNAStringSet("../../Microbiome_DADA2_Results/ASVs-no-contam.fa")
library(ape)

otu_rel<-read.delim("../../Microbiome_DADA2_Results/99_otu_taxonomy.txt",head=F)
library(castor)
cast_tree<-phyloseq::read_tree("../../Microbiome_DADA2_Results/ASV_retree2_maxit_0fasttree_gtr.tree")

tree_remove<-cast_tree$tip.label[!grepl("ASV",cast_tree$tip.label)]
cast_tree<-drop.tip(cast_tree,tree_remove)
library(phyloseq)


ASV_list<-cast_tree$tip.label[!grepl("ASV",cast_tree$tip.label)]
library(phyloseq)
metdat<-as.data.frame(read.delim("../../Microbiome_DADA2_Results/metadata.txt"))
metdat$AgeClass<-ifelse(metdat$Species=="Antarctica",ifelse(grepl("J",metdat$Sample),"Juvenile","Adult"),NA)
meta_reformed<-metdat[metdat$Sample %in% intersect(metdat$Sample,colnames(count_tab)),]
meta_reformed$AgeClass[meta_reformed$Species=="Willana"]<-"Adult"
count_tab<-count_tab[,colnames(count_tab) %in% meta_reformed$Sample]
meta_reformed$SampleType[meta_reformed$Sample=="CCWater4"]<-"Seawater"

rownames(meta_reformed)<-meta_reformed$Sample
count_tab_phy<-otu_table(count_tab,taxa_are_rows = T)

tax_tab_phy_lotus<-tax_table(tax_tab_lotus)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy_lotus, sample_data(meta_reformed),fastas,cast_tree)
ASV_physeq<-prune_taxa(tax_tab2keep,ASV_physeq)

library(tidyverse)
keep<-meta_reformed %>% 
  filter(!Location %in% c("Control")) %>% 
  filter(!SampleType %in% c("Pooled","Unknown")) %>% 
  filter(!grepl("RODI",Sample)) %>% 
  filter(!Species =="ISCA") %>%
  filter(SamplingYear == "2020")


ASV_physeq<-prune_samples(ASV_physeq@sam_data$Sample[!ASV_physeq@sam_data$Location =="Munida"],ASV_physeq)
ASV_physeq<-prune_samples(ASV_physeq@sam_data$Sample[!ASV_physeq@sam_data$Species =="Willana"],ASV_physeq)

ASV_physeq<-prune_samples(keep$Sample, ASV_physeq)

ASV_physeq<-rarefy_even_depth(ASV_physeq,7000)


metadata_hier<-meta_reformed[meta_reformed$Sample %in% ASV_physeq@sam_data$Sample,]

library(devtools)
library(tidyverse)            # cohesive and intuitive collection of packages for data science
library(sp)                   #  creating spatial objects from tables, CRS(), merge()
library(rgdal)                # reading and writing shapefiles, projection
library(raster)               # working with rasterdata, also for calculating area (loads sp)
library(rgeos)                # for function gBuffer to create buffer
library(sf)                   # simple feature support, succeeds sp
library(tmap)
library(marmap)
mean_vel <- raster("../../Microbiome_biogeography/EnvironmentalLayers/Present.Surface.Current.Velocity.Mean.tif.BOv2_1.tif")
sal <- raster("../../Microbiome_biogeography/EnvironmentalLayers/Present.Surface.Salinity.Mean.asc")
temp <- raster("../../Microbiome_biogeography/EnvironmentalLayers/Present.Surface.Temperature.Mean.asc")
tidal_range<-raster("../../Microbiome_biogeography/EnvironmentalLayers/ECDS0243-001-V1.0_Annual average cycle amplitude.tif")
sst_extent_NZ<-raster("../../Kelp_Rafting_MS/fill_stats_bsst_good.tif")
sst_extent_NZ<-terra::crop(sst_extent_NZ,c(163,175.4,-52,-39))

temp <- crop(temp, extent(sst_extent_NZ))
temp <- resample(temp,sst_extent_NZ)
sal <- crop(sal, extent(sst_extent_NZ))
sal <- resample(sal,sst_extent_NZ)
mean_vel <- crop(mean_vel, extent(sst_extent_NZ))
mean_vel <- resample(mean_vel,sst_extent_NZ)
tidal_range<-crop(tidal_range,extent(sst_extent_NZ))
tidal_range <- resample(tidal_range,sst_extent_NZ)
library(terra)
preds<-stack(mean_vel,sal,temp,tidal_range)

ASV_kelp_richness<-estimate_richness(ASV_physeq)
ASV_kelp_richness$Sample<-ASV_physeq@sam_data$Sample
ASV_kelp_richness<-inner_join(ASV_kelp_richness,sample_data(ASV_physeq),by="Sample")
t.test(ASV_kelp_richness$Observed[ASV_kelp_richness$SampleType=="Blade"],
       ASV_kelp_richness$Observed[ASV_kelp_richness$SampleType=="Palmate"],alternative="greater")

library(ggsignif)


richplot<-ggplot(ASV_kelp_richness) +
  aes(x = Location, y = Observed,fill=SampleType)+
  geom_boxplot() +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  theme_minimal() +  theme(
    panel.background = element_rect(fill='white'), #transparent panel bg
    plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    legend.background = element_rect(fill='white'), #transparent legend bg
    legend.box.background = element_rect(fill='white'), #transparent legend panel
    text = element_text(size=18,face="bold"))# +


aggregate(ASV_kelp_richness$Observed~ASV_kelp_richness$SampleType,FUN=mean)
kruskal.test(Observed~SampleType,data=ASV_kelp_richness[ASV_kelp_richness$SampleType %in% c("Palmate","Blade"),])
kruskal.test(Observed~AgeClass,data=ASV_kelp_richness[ASV_kelp_richness$SampleType %in% c("Palmate","Blade"),])

kruskal.test(Observed~Species,data=ASV_kelp_richness[ASV_kelp_richness$Species %in% c("Seawater","Antarctica"),])
kruskal.test(Observed~Species,data=ASV_kelp_richness[ASV_kelp_richness$Species %in% c("Substrate","Antarctica"),])
richplot


ASV_physeq@sam_data$Uplift<-ifelse(ASV_physeq@sam_data$Location %in% c("Kaikoura","CapeCampbell"),
                                   "NoUplift","Uplift")

ASV_Water<-prune_samples(ASV_physeq@sam_data$Sample[ASV_physeq@sam_data$SampleType=="Seawater"],
                         ASV_physeq)
ASV_Water<-prune_taxa(names(taxa_sums(ASV_Water))[taxa_sums(ASV_Water)>0],ASV_Water)
ASV_Sub<-prune_samples(ASV_physeq@sam_data$Sample[ASV_physeq@sam_data$SampleType=="Substrate"],
                       ASV_physeq)
ASV_Sub<-prune_taxa(names(taxa_sums(ASV_Sub))[taxa_sums(ASV_Sub)>0],ASV_Sub)

ASV_Kelp<-prune_samples(ASV_physeq@sam_data$Sample[ASV_physeq@sam_data$SampleType %in% c("Palmate","Blade")],
                        ASV_physeq)

ASV_Kelp<-prune_taxa(names(taxa_sums(ASV_Kelp))[taxa_sums(ASV_Kelp)>0],ASV_Kelp)

set.seed(234731325)
GP.ord_EQ_All <- ordinate(ASV_physeq, "NMDS", distance="bray",k=2,trymax=30)
x<-as.data.frame(GP.ord_EQ_All$points[,1:2])
x$Sample<-rownames(x)
x<-inner_join(x,metadata_hier,by="Sample")
x$Uplift<-ifelse(x$Location %in% c("Kaikoura","CapeCampbell"),
                 "NoUplift","Uplift")
x$SampleGroup<-paste(x$Location,x$SampleType)
library(ggpubr)
library(ggplot2)
nmds_plot<-ggplot(x) +
  aes(x = MDS1, y = MDS2, fill = SampleType, colour = SampleType) +
  geom_point(shape = "circle", size = 3.45) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  theme_minimal() + stat_ellipse() + facet_wrap(~Uplift) + theme_classic2()

nmds_plot


library(vegan)
bdiv_SW<-betadiver(t(ASV_Water@otu_table@.Data),triangular=F,method=1)
bdiv_res_SW<-(betadisper(bdiv_SW,
                         ASV_Water@sam_data$Uplift,type="centroid",bias.adjust=T))
boxplot(bdiv_res_SW)
anova(bdiv_res_SW)

permanova_water <- adonis2(phyloseq::distance(ASV_Water,method="bray",binary=T) ~ ASV_Water@sam_data$Uplift,
                           permutations=999)

bdiv_sub<-betadiver(t(ASV_Sub@otu_table@.Data),triangular=F,method=1)
bdiv_res_sub<-(betadisper(bdiv_sub,
                          ASV_Sub@sam_data$Uplift,type="centroid",bias.adjust=T))

boxplot(bdiv_res_sub)
anova(bdiv_res_sub)
permanova_sub <- adonis2(phyloseq::distance(ASV_Sub,method="bray",binary=T) ~ ASV_Sub@sam_data$Uplift,
                         permutations=999)
permanova_nested <- adonis2(phyloseq::distance(ASV_Kelp,method="bray",binary=T)~ ASV_Kelp@sam_data$Uplift/ASV_Kelp@sam_data$Location/ASV_Kelp@sam_data$AgeClass/ASV_Kelp@sam_data$SampleType,
                            permutations=999)

bdiv_noup<-betadiver(t(ASV_Kelp@otu_table@.Data),triangular=F,method=1)
bdiv_res_tissue<-(betadisper(bdiv_noup,
                             ASV_Kelp@sam_data$SampleType,type="centroid",bias.adjust=T))
boxplot(bdiv_res_tissue)
anova(bdiv_res_tissue)

bdiv_loc<-betadiver(t(ASV_Kelp@otu_table@.Data),triangular=F,method=1)
bdiv_res_loc<-(betadisper(bdiv_loc,
                          ASV_Kelp@sam_data$Location,type="centroid",bias.adjust=T))

boxplot(bdiv_res_loc)
anova(bdiv_res_loc)
bdiv_kelp<-betadiver(t(ASV_Kelp@otu_table@.Data),triangular=F,method=1)
bdiv_res_density<-(betadisper(bdiv_kelp,
                              ASV_Kelp@sam_data$Uplift,type="centroid",bias.adjust=T))
anova(bdiv_res_density)

bdiv_res_age<-(betadisper(bdiv_kelp,
                          ASV_Kelp@sam_data$AgeClass,type="centroid",bias.adjust=T))
boxplot(bdiv_res_age)
anova(bdiv_res_age)
bdiv_res_tiss<-(betadisper(bdiv_kelp,
                           ASV_Kelp@sam_data$SampleType,type="centroid",bias.adjust=T))

boxplot(bdiv_res_tiss)
anova(bdiv_res_tiss)
bdiv_dens<-betadiver(t(ASV_Kelp@otu_table@.Data),triangular=F,method=1)

dist.to.centroid <- betadisper(bdiv_dens, ASV_Kelp@sam_data$Location,type="centroid",bias.adjust=T)


eq_dfs<-data.frame(x=as.numeric(ASV_Kelp@sam_data$Long),y=as.numeric(ASV_Kelp@sam_data$Lat))
eq_dfs$x[eq_dfs$x==174.294088]<-174.304343
eq_dfs$y[eq_dfs$y==-41.711808]<--41.754868
eq_dfs$y[eq_dfs$y==-42.223364]<--42.243364
presvals_eq <- raster::extract(preds,eq_dfs)
presvals_eq<-as.data.frame(presvals_eq)
presvals_eq$Uplift<-ASV_Kelp@sam_data$Uplift
uplift <- ASV_Kelp@sam_data$Uplift
pc_eq_env<-prcomp(presvals_eq[,1:4])
summary(pc_eq_env)
pc_eq_env<-as.data.frame(pc_eq_env$x)
pc_eq_env$Uplift<-ASV_Kelp@sam_data$Uplift
bat <- getNOAA.bathy(165, 179, -49, -32, res = 1, keep=TRUE) ## Bring in a bathymetric map to allow us to do
dist_ant<-as.matrix(sample_data(ASV_Kelp))
dist_ant<-dist_ant[,c(14,13)]
dist_ant<-apply(X = dist_ant,FUN = as.numeric,MARGIN = 2)
colnames(dist_ant)<-c("long","lat")
dist_ant<-as.data.frame(dist_ant)
trans2 <- trans.mat(bat, min.depth=0, max.depth =-2000) # Convert bathymetric object to transmition matrix with max depth of 2000m
ant_dist2 <- lc.dist(trans2, dist_ant, res = "dist",meters=F) # calculated least cost distance matrix
eq_all_wat<-metadata_hier[!metadata_hier$Location == "Munida",]
eq_all_wat<-eq_all_wat[!eq_all_wat$Species == "Willana",]
eq_all_kelp<-eq_all_wat[eq_all_wat$Species == "Antarctica",]
ASV_EQ_Kelp_Wat_Sub<-prune_samples(eq_all_wat$Sample,ASV_physeq)
ASV_EQ_Kelp<-prune_samples(eq_all_kelp$Sample,ASV_EQ_Kelp_Wat_Sub)
ASV_EQ_Kelp@sam_data$Uplift<-ifelse(ASV_EQ_Kelp@sam_data$Location %in% c("Kaikoura","CapeCampbell"),
                                    "NoUplift","Uplift")
ASV_Up<-prune_samples(ASV_EQ_Kelp@sam_data$Sample[ASV_EQ_Kelp@sam_data$Uplift == "Uplift"],
                      ASV_physeq)

ASV_Up<-prune_taxa(taxa_sums(ASV_Up)>0,ASV_Up)

ASV_NoUp<-prune_samples(ASV_EQ_Kelp@sam_data$Sample[ASV_EQ_Kelp@sam_data$Uplift == "NoUplift"],
                        ASV_physeq)
ASV_NoUp<-prune_taxa(taxa_sums(ASV_NoUp)>0,ASV_NoUp)
eqdbRDA<-rda(bdiv_res_age$distances~.,pc_eq_env[,c(1,2,5)] ,scale=F)
anova(eqdbRDA,by="terms")

varpart_eq<-(varpart(bdiv_res_age$distances,~PC1,~PC2,~Uplift,data=pc_eq_env))
varpart_eq$part$indfract$Adj.R.square<-varpart_eq$part$indfract$Adj.R.square*100
plot(varpart_eq,
     Xnames = c("PC1", "PC2","Uplift"), # name the partitions
     bg = c("seagreen3", "mediumpurple","red"), alpha = 80, # colour the circles
     digits = 3, # only show 2 digits
     cex = 1.5)

permutest(eqdbRDA, pairwise = TRUE, permutations = 999)

paste2 <- function(...,sep="") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  gsub(paste0("(^",sep,"|",sep,"$)"),"",
       gsub(paste0(sep,sep),sep,
            do.call(paste,c(L,list(sep=sep)))))
}
library(microbiome)
tax_tab_lotus_temp<-as.data.frame(tax_tab_lotus)
tax_tab_lotus_temp$domain<-paste2("k__",tax_tab_lotus_temp$domain)
tax_tab_lotus_temp$phylum<-paste2("p__",tax_tab_lotus_temp$phylum)
tax_tab_lotus_temp$class<-paste2("c__",tax_tab_lotus_temp$class)
tax_tab_lotus_temp$order<-paste2("o__",tax_tab_lotus_temp$order)
tax_tab_lotus_temp$family<-paste2("f__",tax_tab_lotus_temp$family)
tax_tab_lotus_temp$genus<-paste2("g__",tax_tab_lotus_temp$genus)
tax_tab_lotus_temp$species<-paste2("s__",tax_tab_lotus_temp$species)
tax_tab_lotus_temp$Rels<-paste(tax_tab_lotus_temp$domain,tax_tab_lotus_temp$phylum,
                               tax_tab_lotus_temp$class,tax_tab_lotus_temp$order,
                               tax_tab_lotus_temp$family,tax_tab_lotus_temp$genus,
                               tax_tab_lotus_temp$species,sep="; ")
tax_tab_lotus_temp$Name<-rownames(tax_tab_lotus_temp)

otutable_kelp<-as.data.frame(ASV_Kelp@otu_table@.Data)
otutable_kelp$taxonomy<-tax_tab_lotus_temp$Rels[match(rownames(otutable_kelp),tax_tab_lotus_temp$Name)]
otutable_kelp$taxonomy<-gsub("\t"," ",otutable_kelp$taxonomy)

#write.table(otutable_kelp,"for_faprotax.tsv",quote = F,sep="\t")


#./collapse_table.py -i for_faprotax.txt -o functional_table_nonnorm_26Apr3.tsv -g FAPROTAX.txt -c "#" -d "taxonomy" --omit_columns 0 --column_names_are_in last_comment_line -r report_nonnorm_apr26.txt -n none -v
functional_table<-read.table("./Submission/Revised_Apr_23/functional_table_nonnorm_26Apr3(1).tsv")
colnames(functional_table)<-c("functionalgroup",colnames(otutable_kelp)[1:195])

forfaprotax<-as.data.frame(read.delim("./Kaikoura_Analysis/for_faprotax.tsv"))
## its because some taxa can have multiple functional assignments
rownames(functional_table)<-functional_table$functionalgroup;functional_table$functionalgroup<-NULL

faprotax_report<-readLines("./Submission/Revised_Apr_23/report_nonnorm_apr26.txt")
faprotax_report<-faprotax_report[122:length(faprotax_report)]
faprotax_report<-faprotax_report[!grepl("#",faprotax_report)]
faprotax_report<-faprotax_report[nzchar(faprotax_report)]
taxa_with_funcs<-unique(faprotax_report)

for (i in seq_along(faprotax_report)) {
  if (grepl("^#", faprotax_report[i])) {
    faprotax_report <- faprotax_report[-seq(max(i-2,1),i)]
  }
}


faprotax_report<-unique(faprotax_report)
forfaprotax<-aggregate(. ~ forfaprotax$taxonomy, data = forfaprotax[,1:ncol(forfaprotax)-1], FUN = sum)
func_counts<-forfaprotax[trimws(forfaprotax$`forfaprotax$taxonomy`) %in% trimws(taxa_with_funcs),]
func_counts$`forfaprotax$taxonomy`<-NULL
mean(colSums(func_counts)/7000) #Rarefaction depth = 7000

ASV_function <- phyloseq(otu_table(functional_table,taxa_are_rows = T),
                              sample_data(ASV_Kelp))

create_core_func_comm<-function(physeq,nReads){
  library(tidyverse)
  library(reshape2)
  library(vegan)
  library(ggsci)
  theme_set(theme_light())
  nReads=nReads
  otu<-physeq@otu_table@.Data
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('otu')
  map<-sample_data(physeq)
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Sample, abun, -otu) %>%
    left_join(map, by = 'Sample') %>%
    group_by(otu, Location) %>%
    summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
              coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
              detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
    group_by(otu) %>%
    summarise(sumF=sum(plot_freq),
              sumG=sum(coreSite),
              nS=length(Location)*2,
              Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 
  
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by='otu') %>%
    transmute(otu=otu,
              rank=Index) %>%
    arrange(desc(rank))
  
  BCaddition <- NULL
  
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  for(i in 2:30){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))  
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')
  
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC=mean(BC)) %>%
    arrange(-desc(MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
  elbow <- which.max(BC_ranked$fo_diffs)
  
  lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))
  
  graph<-ggplot(BC_ranked[1:30,], aes(x=factor(BC_ranked$rank[1:30], levels=BC_ranked$rank[1:30]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
    geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='ranked OTUs',y='Bray-Curtis similarity') +
    annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
    annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")
  
  otu_ranked<-otu_ranked
  last_call<-lastCall
  output<-list(graph,otu_ranked,lastCall)
  output
}

library(dplyr)
asv_func_sampdat<-sample_data(ASV_function)
asv_func_sampdat$Group<-paste(asv_func_sampdat$SampleType,asv_func_sampdat$Specific.Sample)
new_df <- asv_func_sampdat[asv_func_sampdat$SampleType %in% c("Blade","Palmate"),]
new_df <- new_df %>% group_by(Group) %>% slice_sample(n=1)

ASV_function_kelp_noup<-prune_samples(new_df$Sample[new_df$Uplift=="NoUplift"],ASV_function)
ASV_function_kelp_up<-prune_samples(new_df$Sample[new_df$Uplift=="Uplift"],ASV_function)

noup_core_funcs<-create_core_func_comm(ASV_function_kelp_noup,3150)
up_core_funcs<-create_core_func_comm(ASV_function_kelp_up,3150)
noup_core_funcs[[2]][1:2,]
up_core_funcs[[2]][1:2,]


GP.ord_EQ_All <- ordinate(ASV_function, "NMDS", distance="bray",k=4,trymax=100)
x_func<-as.data.frame(GP.ord_EQ_All$points[,1:2])
x_func$Sample<-rownames(x_func)
x_func<-inner_join(x_func,metadata_hier,by="Sample")
x_func$Uplift<-ifelse(x_func$Location %in% c("Kaikoura","CapeCampbell"),
                 "NoUplift","Uplift")

x_func$SampleGroup<-paste(x_func$Location,x_func$SampleType)
library(ggplot2)

library(ggplot2)
ggplot(x_func) +
  aes(x = MDS1, y = MDS2, fill = Uplift, colour = Uplift,group=Uplift) +
  geom_point(shape = "circle", size = 3.45) +
  scale_fill_brewer(palette = "Dark2", direction = 1) +
  scale_color_brewer(palette = "Dark2", direction = 1) +
  theme_minimal() + stat_ellipse()


library(microbiome)
ASV_function_kelp<-ASV_function#prune_samples(ASV_Kelp@sam_data$Sample, ASV_function)
dist_dat_fun<-phyloseq::distance(ASV_function,method="bray",binary=T)
bdiv_res_fun<-(betadisper(dist_dat_fun,ASV_function_kelp@sam_data$Uplift))


bdiv_dist_fun<-data.frame(Distance=bdiv_res_fun$distances)

bdiv_dist_fun$pop <-ASV_function_kelp@sam_data$Location

xy.df<-as.data.frame(t(combn(c("Kaikoura","CapeCampbell","Waipapa","Ward"), 2)))
xy.list <- split(xy.df, seq(nrow(xy.df)))
xy.list <- lapply(xy.list,unlist)
bdiv_dist<-data.frame(Distance=bdiv_res_loc$distances)

kelp_unifrac<-phyloseq::distance(ASV_Kelp,method="unifrac")
bdiv_upl<-betadiver(t(ASV_Kelp@otu_table@.Data),triangular=F,method=1)
bdiv_upl<-(betadisper(bdiv_upl,
                      ASV_Kelp@sam_data$Location,type="centroid",bias.adjust=T))

bdiv_unifrac<-(betadisper(kelp_unifrac,
                          ASV_Kelp@sam_data$Location,type="centroid",bias.adjust=T))
bdiv_dist$pop <-ASV_Kelp@sam_data$Location
bdiv_dist$UniFrac<-bdiv_unifrac$distances

BDivPlot_Unifrac<-ggplot(bdiv_dist) +
  aes(x = pop, y = UniFrac,fill = pop) +
  geom_boxplot() +
  theme_classic2()+ theme(text = element_text(size=18,face="bold"),legend.position = "none") + geom_signif(comparisons = xy.list, step_increase = 0.1,
                                                                                                           map_signif_level=TRUE,
                                                                                                           test="wilcox.test") +ylab("UniFrac Distance to Centroid") +
  scale_fill_manual(
    values = c(CapeCampbell = "#117733",
               Kaikoura = "#02ba48",
               Waipapa = "#6699CC",
               Ward = "#0052a3")
  ) + ylim(0.3,0.9)
BDivPlot_Unifrac


BDivPlot_func<-ggplot(bdiv_dist_fun) +
  aes(x = pop, y = Distance,fill = pop) +
  geom_boxplot() +
  theme_classic2()+ theme(text = element_text(size=18,face="bold"),legend.position = "none") + 
  geom_signif(comparisons = xy.list, step_increase = 0.02,
map_signif_level=TRUE,
test="wilcox.test") +ylab("Functional Distance to Centroid") +
  scale_fill_manual(
    values = c(CapeCampbell = "#117733",
               Kaikoura = "#02ba48",
               Waipapa = "#6699CC",
               Ward = "#0052a3")
  ) + ylim(0,1)
BDivPlot_Unifrac
BDivPlot_func
library(cowplot)
uni_func<-align_plots(BDivPlot_Unifrac,BDivPlot_func)

uni_func<-plot_grid(uni_func[[1]],uni_func[[2]],nrow=2)

create_core_comm<-function(physeq,nReads){
  library(tidyverse)
  library(reshape2)
  library(vegan)
  library(ggsci)
  theme_set(theme_light())
  nReads=nReads
  otu<-physeq@otu_table@.Data
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('otu')
  map<-sample_data(physeq)
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Sample, abun, -otu) %>%
    dplyr::left_join(map, by = 'Sample') %>%
    group_by(otu, Location) %>%
    summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
              coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
              detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
    group_by(otu) %>%
    dplyr::summarise(sumF=sum(plot_freq),
                     sumG=sum(coreSite),
                     nS=length(Location)*2,
                     Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 
  
  otu_ranked <- occ_abun %>%
    dplyr::left_join(PresenceSum, by='otu') %>%
    transmute(otu=otu,
              rank=Index) %>%
    arrange(desc(rank))
  
  BCaddition <- NULL
  
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  for(i in 2:500){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- dplyr::left_join(BCaddition, df_a, by=c('x_names'))  
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- dplyr::left_join(BCaddition,df_full, by='x_names')
  
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    dplyr::summarise(MeanBC=mean(BC)) %>%
    arrange(-desc(MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- dplyr::left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
  elbow <- which.max(BC_ranked$fo_diffs)
  
  lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))
  
  graph<-ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
    geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='ranked OTUs',y='Bray-Curtis similarity') +
    annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
    annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")
  
  otu_ranked<-otu_ranked
  last_call<-lastCall
  output<-list(graph,otu_ranked,lastCall,BC_ranked)
  output
}
ASV_Up_core<-create_core_comm(ASV_Up,7000)
ASV_NoUp_core<-create_core_comm(ASV_NoUp,7000)

ASV_Up_core[[2]][1:4,1]
tax_tab[rownames(tax_tab) %in%ASV_Up_core[[2]][1:4,1],]
ASV_NoUp_core[[2]][1:3,]
tax_tab[rownames(tax_tab) %in%ASV_NoUp_core[[2]][1:2,1],]

library(iCAMP)#;library(dada2)
library(phyloseq)
Uplift_BNTI<-iCAMP::bNTI.cm(t(ASV_Up@otu_table@.Data),
                            cophenetic.phylo(ASV_Up@phy_tree),
                            nworker = 3,sig.index="bNTI",rand = 10)

Uplift_RC<-iCAMP::RC.cm(t(ASV_Up@otu_table@.Data),nworker = 3,rand = 10)

NoUplift_BNTI<-iCAMP::bNTI.cm(t(ASV_NoUp@otu_table@.Data),
                            cophenetic.phylo(ASV_NoUp@phy_tree),
                            nworker = 3,sig.index="bNTI",rand = 10)
NoUplift_RC<-iCAMP::RC.cm(t(ASV_NoUp@otu_table@.Data),nworker = 3,rand = 10)


NoUP_RC_Dat<-reshape2::melt(NoUplift_RC$index)
UP_RC_Dat<-reshape2::melt(Uplift_RC$index)

RC_Dat<-rbind(NoUP_RC_Dat,UP_RC_Dat)

Up_BNTI_dat<-reshape2::melt(Uplift_BNTI$index)
NoUp_BNTI_dat<-reshape2::melt(NoUplift_BNTI$index)

BNTI_Dat <- rbind(NoUp_BNTI_dat,Up_BNTI_dat)


BNTI_Dat<-rbind(NoUp_BNTI_dat,Up_BNTI_dat)

Comm_Proc_Dat<-cbind(BNTI_Dat,RC_Dat)
colnames(Comm_Proc_Dat)<-c("Var1","Var2","BNTI",
                           #"Var1x","Var2x","BNRI",
                           "Var1y","Var2y","RC")
Comm_Proc_Dat<-Comm_Proc_Dat[,c(1,2,3,6)]
Comm_Proc_Dat$Uplift<-c(rep("LowUplift",dim(NoUp_BNTI_dat)[1]),rep("HighUplift",dim(Up_BNTI_dat)[1]))

Comm_Proc_Dat$Pop1<-ASV_Kelp@sam_data$Location[match(Comm_Proc_Dat$Var1, 
                                                     ASV_Kelp@sam_data$Sample)]
Comm_Proc_Dat$Pop2<-ASV_Kelp@sam_data$Location[match(Comm_Proc_Dat$Var2, 
                                                     ASV_Kelp@sam_data$Sample)]
Comm_Proc_Dat$Tiss1<-ASV_Kelp@sam_data$SampleType[match(Comm_Proc_Dat$Var1, 
                                                     ASV_Kelp@sam_data$Sample)]
Comm_Proc_Dat$Tiss2<-ASV_Kelp@sam_data$SampleType[match(Comm_Proc_Dat$Var2, 
                                                     ASV_Kelp@sam_data$Sample)]
Comm_Proc_Dat$GroupSample<-ifelse(Comm_Proc_Dat$Tiss1==Comm_Proc_Dat$Tiss2,"Same","Different")
Comm_Proc_Dat<-Comm_Proc_Dat[Comm_Proc_Dat$GroupSample=="Same",]

library(ggpubr)
BNTI_Plot<-ggplot(Comm_Proc_Dat) +
  aes(x = BNTI, fill = Uplift, colour = Uplift) +
  geom_density(adjust = 1L,alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = 1) +
  scale_color_brewer(palette = "Set1", direction = 1) +
  theme_classic2() + ylab("Density")+
  geom_vline(xintercept = -2, linetype="dotted", 
                                  color = "black", size=1.5) +
  theme( #transparent legend panel
    text = element_text(size=18,face="bold")) + xlim(-7,2)



Comm_Proc_Dat_NonSel<-Comm_Proc_Dat[abs(Comm_Proc_Dat$BNTI) < 2,]
RC_Plot<-ggplot(Comm_Proc_Dat_NonSel) +
  aes(x = RC, fill = Uplift, colour = Uplift) +
  geom_density(adjust = 1L,alpha=0.7) +
  scale_fill_brewer(palette = "Set1", direction = 1) +
  scale_color_brewer(palette = "Set1", direction = 1) +
  theme_classic2() + ylab("Density")+ xlim(-1,1) +
  geom_vline(xintercept = 0.95, linetype="dotted", 
             color = "black", size=1.5)+ 
  geom_vline(xintercept = -.95, linetype="dotted", 
             color = "black", size=1.5)  +
  theme(text = element_text(size=18,face="bold")) + xlab("Raup-Crick") 

bnti_plots<-ggarrange(BNTI_Plot,RC_Plot, nrow=2)


RC_Plot
library(tidyverse)
PropHomSel_LowUp<-nrow(Comm_Proc_Dat %>% filter(BNTI <= -2) %>% filter(Uplift=="LowUplift"))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="LowUplift"))
PropHetSel_LowUp<-nrow(Comm_Proc_Dat %>% filter(BNTI >= 2) %>% filter(Uplift=="LowUplift"))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="LowUplift"))
round(PropHomSel_LowUp*100,1)

PropHomSel_HighUp<-nrow(Comm_Proc_Dat %>% filter(BNTI <= -2) %>% filter(Uplift=="HighUplift"))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="HighUplift"))
PropHetSel_HighUp<-nrow(Comm_Proc_Dat %>% filter(BNTI >= 2) %>% filter(Uplift=="HighUplift"))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="HighUplift"))
round(PropHomSel_HighUp*100,1)





PropDrift_HighUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                         filter(Uplift=="HighUplift") %>% 
                         filter(abs(RC)<0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="HighUplift"))


PropDrift_LowUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                         filter(Uplift=="LowUplift") %>% 
                         filter(abs(RC)<0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="LowUplift"))


PropDispLim_LowUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                        filter(Uplift=="LowUplift") %>% 
                        filter(RC>=0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="LowUplift"))

PropDispLim_HighUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                          filter(Uplift=="HighUplift") %>% 
                          filter(RC>=0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="HighUplift"))
round(PropDispLim_HighUp*100,1)

PropDispHom_HighUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                           filter(Uplift=="HighUplift") %>% 
                           filter(RC<=-0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="HighUplift"))

PropDispHom_LowUp<-nrow(Comm_Proc_Dat %>% filter(abs(BNTI) < 2) %>%
                           filter(Uplift=="LowUplift") %>% 
                           filter(RC<=-0.95))/
  nrow(Comm_Proc_Dat %>% filter(Uplift=="LowUplift"))
PropDrift_LowUp*100

