####################### CAFE OUTPUTS ##################################
## ********************************************************************
# cd data/
# for i in {1..5}; do
# cafe5 -i orthGroups3.tab -t datedSpeciesTree.txt -k $i -p -o ../output/withK$i;
# done

# ## Error Model --------
# cafexp -i orthGroups.tab -t datedSpeciesTree2.txt -e

### Test model fit -----------
library(ggplot2)
# Log-likelihood values including the base model
log_likelihoods <- c(139649, 58104.7, 129994, 57288.1, 58010.4, 57773.8)
model_labels <- c("Base Model", paste0("Gamma k=", 1:5))
num_params <- c(1, 2, 3, 4, 5, 6)  # Base model + increasing complexity
AIC_values <- -2 * log_likelihoods + 2 * num_params
BIC_values <- -2 * log_likelihoods + log(length(log_likelihoods)) * num_params  # BIC accounts for sample size
model_data <- data.frame(Model = model_labels, LogL = log_likelihoods, AIC = AIC_values, BIC = BIC_values, NumParams = num_params)
# Plot Log-Likelihood
ggplot(model_data, aes(x = NumParams, y = LogL, label = Model)) +
  geom_line(color = "dodgerblue") +
  geom_point(size = 3, color = "dodgerblue") +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  ggtitle("Log-Likelihood vs. Model Complexity") +
  xlab("Number of Parameters") +
  ylab("Log-Likelihood")

# Plot AIC and BIC together
ggplot(model_data, aes(x = NumParams)) +
  geom_line(aes(y = AIC, color = "AIC")) +
  geom_point(aes(y = AIC, color = "AIC"), size = 3) +
  geom_line(aes(y = BIC, color = "BIC")) +
  geom_point(aes(y = BIC, color = "BIC"), size = 3) +
  geom_text(aes(y = AIC, label = Model), vjust = -1, size = 3, color = "red") +  
  geom_text(aes(y = BIC, label = Model), vjust = 1, size = 3, color = "purple") +
  scale_color_manual(values = c("AIC" = "red", "BIC" = "purple")) +
  theme_minimal() +
  ggtitle("Model Selection: AIC & BIC vs. Model Complexity") +
  xlab("Number of Parameters") +
  ylab("AIC / BIC") +
  theme(legend.position = "right")

# 1. Log-Likelihood Plot
# •	The higher the log-likelihood, the better the model fits the data.
# •	If log-likelihood stops increasing significantly, adding more complexity (higher  k ) may not be useful.
# 2.	AIC/BIC Plot
# •	The lowest AIC/BIC value indicates the best model.
# •	AIC is more forgiving (allows more parameters).
# •	BIC penalizes complexity (prefers simpler models).
# •	If AIC keeps decreasing, the model may be overfitting.
# •	If AIC and BIC agree on the best model, it’s likely the optimal choice.


# ** Base model shows the best fit **
# Lambda: 0.011789923099377
# Maximum possible lambda for this topology: 0.0238095
# ** with error model **
# Lambda: 0.010903032180589
# Epsilon: 0.023207
# Maximum possible lambda for this topology: 0.0238095

# Input species tree
library(ape)
sp_tree<-read.tree("data/datedSpeciesTree.txt")

### Base Model family trees (all) -------------
# ** this file contains all gene family trees; the output from CAFE tree reconstruction **
bs_model<-read.nexus("data/cafe/run3/base/Base_asr.tre")


### Model Clades ---------------
# ** This file contains the number of families expanded/contracted at each node and clade in the species tree **
clad<-data.table::fread("data/cafe/run3/base/Base_clade_results.txt")

#the probabilities calculated for each clade and significant family
branch_prob<-data.table::fread("data/cafe/run3/base/Base_branch_probabilities.tab")

#likelihood of families
fam_lkh<-data.table::fread("data/cafe/run3/base/Base_family_likelihoods.txt")

### GAMMA models [Finding the best lambda for each family] --------
lamb2<-read.nexus("data/cafe/run3/withK2/Gamma_asr.tre")
bintree<-multi2di(lamb2)
#cons_tree<-consensus(bintree)

### 2.4 **All significant families** -------------

### Significant gene families --------------
# ** this file shows if the rate of evolution of each gene family is significant or not under the BD model **
# ** the file only shows the significance at 0.05 but 0.01 can be added; and was used downstream here **
fam<-data.table::fread("data/cafe/run3/base/Base_family_results.txt")
# adding significant families for p=0.01
fam$sig_at_0.01<-ifelse(fam$pvalue<=0.01,"y","n")
colnames(fam)<-c("FamilyID","pvalue","sig_at_0.05","sig_at_0.01")

sig_fam2.1<-data.table::fread("data/cafe/sig_fams_0.01.txt",h=T)
sig_fam2.1<-data.frame(sig_fam2.1)

### Distribution of significanly evolving gene families -------
dpar<-par(no.readonly = T)
boxplot(sig_fam2.1[,3:7],ylim=c(0,50))
par(mfrow=c(1,5))
for(i in 1:5){
  mn<-round(mean(sig_fam2.1[,i+2]),1)
  plot(0,xlim=c(-2,50),ylim=c(0,.3),type="n",main=paste0(colnames(sig_fam2.1)[i+2]," (mean ",mn,")"),ylab="density",xlab="number of genes")
  polygon(density(sig_fam2.1[,i+2]),col=2,border=F)
  abline(v=mn,col=4)
}
par(dpar)


# filter significant gene families
og_gene_count<-read.table("data/Orthogroups.GeneCount.tsv",h=T)
sig_og<-og_gene_count[match(fam$FamilyID,og_gene_count$Orthogroup),]
sig_og$pval<-fam$pvalue

sig_og_0.05<-sig_og[sig_og$pval<=0.05,]
sig_og_0.01<-sig_og[sig_og$pval<0.01,]


#drop families found only in one species
zero_count<-apply(sig_og_0.01[,2:6],1,function(x)sum(x==0))
sig_og_0.01<-sig_og_0.01[zero_count<3,]

# drop AT=0
sig_og_0.01_AT_noZero<-sig_og_0.01[sig_og_0.01$Athaliana_Araport11.priTranscripts>0,]



## *************************************
## CHARACTERIZE SIGNIFICANT GENE FAMILIES (ORTHO GROUPS) FROM ORTHOFINDER AND CAFE5 ANALYSES -------------------
## Use UniProt (with python) to retrieve functions --------
# ** See STATS section to find the final file **
# 
fun_added<-data.table::fread("data/fast_evolving_arabidopsis_families_classification_mac.tsv",h=T)
ids<-fun_added$at_genes
ids<-unlist(stringr::str_split(ids,"\\,"))
write.table(ids,"data/gene_ids_to_search_functions.txt",sep="\t",row.names=F, quote=F,col.names = F)
# ** infor retrieved with uni_prot.py script **
uni_prot<-read.csv("data/gene_function_annotations.csv")
categories<-uni_prot$Protein.names
categories[grep("defensin-like protein",categories,ignore.case = T)]<-"Defensin-like protein"
categories[grep("(thale cress)",categories,ignore.case = T)]<-"hypothetical protein"
categories[grep("(rape)",categories,ignore.case = T)]<-"hypothetical protein"
categories[grep("Disease resistance protein",categories,ignore.case = T)]<-"Disease resistance protein"
categories[grep("(wild Malaysian banana)",categories,ignore.case = T)]<-"hypothetical protein"
categories[grep("Amino acid transporter",categories,ignore.case = T)]<-"Amino acid transporter"
categories[grep("AGAMOUS-like",categories,ignore.case = T)]<-"AGAMOUS-like"
categories[grep("gamma-bisabolene synthase",categories,ignore.case = T)]<-"(Z)-gamma-bisabolene synthase"
categories[grep("Ac-like transposase",categories,ignore.case = T)]<-"Ac-like transposase"
categories[grep("Aminotransferase-like",categories,ignore.case = T)]<-"Aminotransferase-like protein"
categories[grep("Zinc knuckle",categories,ignore.case = T)]<-"Zinc knuckle family protein"
categories[grep("Zinc finger",categories,ignore.case = T)]<-"Zinc finger-(like) protein"
categories[grep("Wall-associated",categories,ignore.case = T)]<-"Wall-associated receptor kinase-like"
categories[grep("Uncharacterized protein",categories,ignore.case = T)]<-"Uncharacterized protein"
categories[grep("Uncharacterized mitochondrial protein",categories,ignore.case = T)]<-"Uncharacterized mitochondrial protein"
categories[grep("UDP-Glycosyltransferase",categories,ignore.case = T)]<-"UDP-Glycosyltransferase"
categories[grep("Tropinone reductase homolog",categories,ignore.case = T)]<-"Tropinone reductase homolog"
categories[grep("Ubiquitin-like",categories,ignore.case = T)]<-"Ubiquitin-like"
categories[grep("Tropinone",categories,ignore.case = T)]<-"Tropinone"
categories[grep("Transmembrane protein",categories,ignore.case = T)]<-"Transmembrane protein"
categories[grep("Thionin",categories,ignore.case = T)]<-"Plant thionin family protein"
categories[grep("Transcription factor",categories,ignore.case = T)]<-"Transcription factor"
categories[grep("Terpenoid synthase",categories,ignore.case = T)]<-"Terpenoid synthase"
categories[grep("Subtilisin-like protease",categories,ignore.case = T)]<-"Subtilisin-like protease"
categories[grep("SKP1-like protein",categories,ignore.case = T)]<-"SKP1-like protein"
categories[grep("Serine carboxypeptidase",categories,ignore.case = T)]<-"Serine carboxypeptidase(-like) family"
categories[grep("S-protein homolog",categories,ignore.case = T)]<-"S-protein homolog"
categories[grep("Receptor-like protein kinase",categories,ignore.case = T)]<-"Receptor-like protein kinase"
categories[grep("Receptor-like protein(?! kinase)", categories, perl = TRUE)]<-"Receptor-like protein"
categories[grep("Receptor like protein(?! kinase)", categories, perl = TRUE)]<-"Receptor-like protein"
categories[grep("F-box protein", categories, perl = TRUE)]<-"F-box protein"
categories[grep("F-box", categories, perl = TRUE)]<-"F-box domain"
categories[grep("putative", categories, ignore.case = T)]<-gsub("Putative ","",ignore.case = T,x=categories[grep("putative", categories, ignore.case = T)])
categories[grep("probable", categories, ignore.case = T)]<-gsub("probable ","",ignore.case = T,x=categories[grep("probable", categories, ignore.case = T)])
categories[grep("Polyubiquitin",categories,ignore.case = T)]<-"Polyubiquitin"
categories[grep("Pentatricopeptide repeat",categories,ignore.case = T)]<-"Pentatricopeptide repeat-containing protein"
categories[grep("MATH domain and coiled-coil domain-containing protein",categories,ignore.case = T)]<-"MATH domain and coiled-coil domain-containing protein"
categories[grep("LRR receptor-like serine",categories,ignore.case = T)]<-"LRR receptor-like serine/threonine-protein kinase"
categories[grep("Long-chain-alcohol O-fatty-acyltransferase",categories,ignore.case = T)]<-"Long-chain-alcohol O-fatty-acyltransferase"
categories[grep("Leucine-rich repeat",categories,ignore.case = T)]<-"Leucine-rich repeat(-containing) protein"
categories[grep("Jacalin-related lectin",categories,ignore.case = T)]<-"Jacalin-related lectin"
categories[grep("L-type lectin-domain containing receptor kinase",categories,ignore.case = T)]<-"L-type lectin-domain containing receptor kinase"
categories[grep("inactive receptor kinase",categories,ignore.case = T)]<-"inactive receptor kinase"
categories[grep("Heterodimeric geranylgeranyl pyrophosphate synthase large subunit",categories,ignore.case = T)]<-"Heterodimeric geranylgeranyl pyrophosphate synthase large subunit"
categories[grep("Glycosyl hydrolase family",categories,ignore.case = T)]<-"Glycosyl hydrolase family"
categories[grep("Geranylgeranyl pyrophosphate synthase",categories,ignore.case = T)]<-"Geranylgeranyl pyrophosphate synthase"
categories[grep("E3 ubiquitin-protein ligase SINA-like",categories,ignore.case = T)]<-"E3 ubiquitin-protein ligase SINA-like"
categories[grep("Cytochrome P450",categories,ignore.case = T)]<-"Cytochrome P450"
categories[grep("Cysteine/histidine",categories,ignore.case = T)]<-"Cysteine/histidine-rich C1 domain protein"
categories[grep("Cysteine-rich repeat secretory protein",categories,ignore.case = T)]<-"Cysteine-rich repeat secretory protein"
categories[grep("B3 domain-containing protein",categories,ignore.case = T)]<-"B3 domain-containing protein"
categories[grep("Auxin-responsive",categories,ignore.case = T)]<-"Auxin-responsive protein"
categories[grep("Mitochondrial transcription termination",categories,ignore.case = T)]<-"Mitochondrial transcription termination protein"
categories[grep("Ubiquitin transferase",categories,ignore.case = T)]<-"Ubiquitin transferase"
categories[grep("Ubiquitin carboxyl-terminal hydrolase-related protein",categories,ignore.case = T)]<-"Ubiquitin carboxyl-terminal hydrolase-related protein"
write.table(data.frame(sort(unique(categories),decreasing = T)),"/Users/piyalkaru/Desktop/DDORF/Ann/general/family_discovery/unique_functions.txt",row.names = F,quote = F)
## Categorized with python script.
uniPcat<-data.table::fread("data/categorized_gene_functions_extended.csv",h=T)
## add them into orthogroups/families
ogs_genes<-fun_added[,c(2,12)]
og_tm<-NULL
for(i in 1:nrow(ogs_genes)){
  tm<-cbind(unlist(ogs_genes[i,1]),unlist(stringr::str_split(ogs_genes[i,2],"\\,")))
  og_tm<-rbind(og_tm,tm)
}
og_tm<-data.frame(OG=og_tm[,1],at_gene=og_tm[,2])
og_tm$at_gene<-stringr::str_split_fixed(og_tm$at_gene,"\\.",n=2)[,1]
uniPcat$gene_id<-stringr::str_split_fixed(uniPcat$gene_id,"\\.",n=2)[,1]

uniPcat$gene_id <- trimws(as.character(uniPcat$gene_id))
og_tm$at_gene <- trimws(as.character(og_tm$at_gene))
uniPcat$gene_id <- tolower(uniPcat$gene_id)
og_tm$at_gene <- tolower(og_tm$at_gene)
uniPcat$gene_id <- as.character(uniPcat$gene_id)
og_tm$at_gene <- as.character(og_tm$at_gene)

uniPcat$family_id<-og_tm[match(uniPcat$gene_id,og_tm$at_gene),1]
# Regroup into orthogroups
ogs<-unique(uniPcat$family_id)
og_cat<-lapply(ogs,function(x){
  c(x,sort(unique(unlist(uniPcat[uniPcat$family_id==x,"Category"]))))
})
library(plyr)
reOG <- do.call(rbind.fill, lapply(og_cat, function(x) as.data.frame(t(x))))
## add that to the orthogroup table
ogList<-data.table::fread("data/fast_evolving_arabidopsis_families_classification_mac.tsv",h=T)
ogList<-cbind(ogList[,c(2:8,11,12)],reOG[match(ogList$Family.ID,reOG$V1),-1])
ogList[ogList=="Biosynthesis" | ogList=="Metabolic & enzymatic proteins"]<-"Biosynthesis & Metabolism"
ogList[ogList=="Cell Wall Modification"]<-"Structural Dynamics"
ogList[ogList=="Defense & Resistance" | ogList=="Stress Response"]<-"Defense & Stress response"
ogList[ogList=="Genetic Information Processing" | ogList=="Regulation of Gene Expression"]<-"Genetic Information & Regulation"
ogList[ogList=="Signaling & Transport" | ogList=="P-P Interaction"]<-"Signaling and Transport"
ogList[ogList=="Uncharacterized" | ogList=="Uncategorized"]<-"Uncategorized"
write.csv(ogList,"data/ogList_with_functions.csv",row.names = F)


### Plot Pie of FAMILY categories ------------
library(ggplot2)
library(dplyr)
library(ggforce)

fun_added <- read.csv("data/ogList_with_functions.csv")
df_summary <- fun_added[, c(9, 10)] %>%
  group_by(FunctionalCategory) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(
    percentage = count / sum(count),
    percent_label = paste0(round(percentage * 100,1), "%"),
    ymax = cumsum(percentage),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2  # middle of each slice
  )

# Plot
gp<-ggplot(df_summary, aes(ymax = ymax, ymin = ymin, xmax = 2.5, xmin = 1.5, fill = FunctionalCategory)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(x = 2, y = ymid, label = percent_label),
    color = "black",
    size = 4
  ) +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Functional Categories")

print(gp)

### barplot of gene categories per species ------------
library(ggplot2)
library(dplyr)
library(tidyr)
bar_df<-fun_added[,c(2:5,10)]
bar_long <- bar_df %>%
  pivot_longer(
    cols = c(Aalpina, Ahalleri, Alyrata, Athaliana),
    names_to = "Species",
    values_to = "GeneCount"
  )
bar_summary <- bar_long %>%
  group_by(FunctionalCategory, Species) %>%
  summarise(TotalGenes = sum(GeneCount), .groups = "drop")
ggplot(bar_summary, aes(x = FunctionalCategory, y = TotalGenes, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x=NULL,y = "Number of Genes") +
  theme_minimal() +
  theme(
    plot.margin = margin(1, 1, 2, 1, "cm"),  # top, right, bottom, left
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )


## CAFE BASE MODEL lambda -----------
library(ape)
base_tree<-read.tree("data/cafe/run3/base/Base_asr.tre")
def_gene_ids<-data.table::fread("data/defense_and_stress_fast_evol_families_gene_ids_4_species.txt",h=T)
all_231fams<-fun_added
sig_trees<-base_tree[unlist(lapply(all_231fams$Family.ID,function(x)grep(x,names(base_tree))))]

ttree<-sig_trees[[1]]
plot(ttree)
nodelabels(ttree$node.label)

nd_count<-lapply(sig_trees,function(x)x$node.label)
nd_count<-lapply(nd_count,function(x){as.numeric(gsub("\"","",stringr::str_split_fixed(x,n=2,pattern="_")[,2]))})
nd_count<-data.frame(do.call(rbind,nd_count))
colnames(nd_count)<-c("papaya_ara","arabis_arabidopsis","thaliana_others","lyrata_halleri")
colMeans(nd_count)

## plot the distribution of lambda? values for each node
pdf("plots/node_lambda_density_plots2.pdf",w=8,h=5)
par(mfrow=c(2,2))
for(i in 1:ncol(nd_count)){
  x<-nd_count[,i]
  n1<-density(x,adjust=2)
  peak_x <- n1$x[which.max(n1$y)]
  max_y<-max(n1$y)
  plot(n1,axes=F,lwd=3,ylab=NA,xlab=NA,type="n",main=colnames(nd_count)[i],ylim=c(-(max_y/5),max_y))
  polygon(x=n1$x,y=n1$y,col=colorspace::adjust_transparency(2,alpha=peak_x/4),lwd=2)
  lines(x=c(peak_x,peak_x),y=c(0,max_y),lty=2,lwd=3)
  text(x = peak_x , y = -(max_y/10), labels = round(peak_x, 2), font = 2)
}
dev.off()







### GO enrichment assessment ----------------------
library(clusterProfiler)
library(org.At.tair.db)

fun_added <- read.csv("data/ogList_with_functions.csv")
my_genes <- stringr::str_trim(stringr::str_split_fixed(unlist(strsplit(fun_added$at_genes,split=",")),pattern = "\\.",n=2)[,1],"both")

# Run GO enrichment
ego <- enrichGO(gene          = my_genes,
                universe      = keys(org.At.tair.db),
                OrgDb         = org.At.tair.db,
                keyType       = "TAIR",
                ont           = "ALL",         # Biological Process
                pAdjustMethod = "BH",         # Benjamini-Hochberg correction
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.10)

# Extract results and filter for defense/stress
ego_df <- as.data.frame(ego)

# View results
dotplot(ego, showCategory = 50)

# search and add functional categories to GO table
go_ids<-stringr::str_split(ego_df$geneID,pattern = "\\/")

fun_cats<-lapply(go_ids,function(x){
  tm<-NULL
  for(i in seq_along(x)){
    tm<-c(tm,fun_added[grep(x[i],fun_added$at_genes),"FunctionalCategory"])
    vv<-table(tm)
    tm<-names(which.max(vv))
  }
  return(tm)
})


library(dplyr)
library(ggplot2)

ego_df <- as.data.frame(ego)
ego_df$Functional_category<-unlist(fun_cats)


tv<-data.frame(stringr::str_split_fixed(ego_df$GeneRatio,"\\/",n=2))
tv$X1<-as.numeric(tv$X1);tv$X2<-as.numeric(tv$X2)
ego_df$GeneRatio<-tv$X1/tv$X2
write.csv("data/GO_raw_output.csv",row.names=F)

# Create a numbered ID for each GO term (e.g., "GO_1", "GO_2")
ego_df <- ego_df %>%
  group_by(Functional_category) %>%
  arrange(p.adjust) %>%  # Sort by significance within each category
  mutate(
    GO_short = paste0("GO_", row_number()),  # Short label
    Description_full = Description            # Store full description
  ) %>%
  ungroup()

# Maintain order for plotting
ego_df$GO_short <- factor(ego_df$GO_short, levels = unique(ego_df$GO_short))


ggplot(ego_df, aes(x = GeneRatio, y = GO_short, 
                   color = p.adjust, size = Count)) +
  geom_point() +
  facet_grid(Functional_category ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(
    low = "red", high = "blue", 
    name = "Adjusted p-value",
    trans = "log10"
  ) +
  scale_size_continuous(name = "Gene count") +
  labs(x = "Gene Ratio", y = "GO Term") +  # Update axis label
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0)
  )

go_lookup <- ego_df %>%
  select(Functional_category, GO_short, Description_full, p.adjust, GeneRatio, Count) %>%
  arrange(Functional_category, p.adjust)

# save lookup table
write.csv(go_lookup, "data/GO_term_lookup_table.csv", row.names = FALSE)
