#################### Defense and Stress related gene family analysis ###################
#=======================================================================================

########## ANALYSES OF DEFENSE AND STRESS RESPONSE RELATED GENE FAMILIES ################
#=======================================================================================

# gene ids of def&stress gene families
def_gene_ids<-data.table::fread("data/defense_and_stress_fast_evol_families_gene_ids_4_species.txt",h=T)

# THALIANA ----------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Extract gene sequences and save in FASTA  ------------------------------
library(Biostrings)
at_ids<-def_gene_ids[,c(1,5)]
# genome and the gff files
at_genome<-readDNAStringSet("data/TAIR10_allChr.fas")
athal_gff <-data.frame(biomartr::read_gff("data/Araport11_GFF3_genes_transposons.Feb2022.gff"))

for(j in 1:nrow(at_ids)){
  newIds<-stringr::str_split_fixed(unlist(stringr::str_split(at_ids[j,2],pattern="\\,")),pattern="\\.",n=2)[,1]
  gene_list<-list()
  for(i in seq_along(newIds)){
    g_match<-athal_gff[grep(gsub("\\ ","",newIds[i]),athal_gff$attribute),] #AT1G70020
    rr<-which(g_match$type=="gene")[1]
    coords<-range(g_match[rr,4:5])
    chr<-g_match[1,1]
    gene<-at_genome[[grep(chr,names(at_genome))]][min(coords):max(coords)]
    gene_list[[i]]<-gene
  }
  seqinr::write.fasta(gene_list,names = newIds,file.out=paste0("data/DefenseNstress/AT/refseq/",unlist(at_ids[j,1]),".fasta"))
}

### find gene lengths --------
fls<-list.files("data/DefenseNstress/AT/refseq",full.names=T)
at_gen_length<-do.call(rbind,lapply(fls,function(x){dd<-readDNAStringSet(x);return(cbind(names(dd),Biostrings::width(dd)))}))
at_gen_length<-data.frame(id=at_gen_length[,1],length=at_gen_length[,2])
at_gen_length$length<-as.numeric(at_gen_length$length)
at_gen_length$id<-as.character(stringr::str_trim(at_gen_length$id))
write.table(at_gen_length,"data/DefenseNstress/AT/AT_DnS_gene_lengths.txt",row.names=F, quote=F,sep="\t")

## BLAST FOR GENES --------------------------------
reference_dir="data/defensNstress/AT/refSeq"
seq_dir="data/Athal/LR_assemblies/136_LR_fnas"
out_dir="data/defensNstress/AT/blast2/"

mkdir -p $out_dir
for ref_fa in ${reference_dir}/*.fasta; do
ref=`basename $ref_fa | sed 's/.fasta//'`
echo $ref
for file in ${seq_dir}/*.fna; do
samp=`basename $file | sed 's/.fna//'`
echo $samp
blastn -query $ref_fa -subject $file -outfmt '7 std stitle' -out ${out_dir}${ref}_${samp}.blastn ;
done;
done


## CHEKC BLAST HITS  --------------------------
b_files<-list.files("data/DefenseNstress/AT/blast2",full.names=T)
og_ids<-at_ids$Orthogroup

for(j in seq_along(og_ids)){
  og_files<-b_files[grep(og_ids[j],b_files)]
  og_added<-NULL
  for(i in seq_along(og_files)){
    b_hits<-read.table(og_files[i],comment.char = "#",col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","stitle"))
    b_hits$stitle<-stringr::str_split_fixed(gsub(".blastn","",basename(og_files[i])),pattern="_",n=2)[,2]
    og_added<-rbind(og_added,b_hits)
  }
  write.table(og_added,paste0("data/DefenseNstress/AT/blast2/og_sorted/",og_ids[j],"_bhits_v1.txt"),quote=F,row.names=F,sep="\t")
}

## GROUP THE BLAST HITS  -----------------------------
at_gen_length<-read.table("data/DefenseNstress/AT/AT_DnS_gene_lengths.txt",h=T)
og_sorted<-list.files("data/DefenseNstress/AT/blast2/og_sorted/",full.names = T)
pb<-txtProgressBar(min=0,max=length(og_sorted),style = 3,width = 50)
for(m in seq_along(og_sorted)){
  setTxtProgressBar(pb,m)
  all_blast<-read.table(og_sorted[m],header = T)
  all_blast<-all_blast[all_blast$pident>=90,]
  og<-gsub("bhits_v1.txt","",basename(og_sorted[m]))
  all_blast$assembly<- gsub("PacbioHiFiAssembly_","",all_blast$stitle)
  assembly<-unique(all_blast$assembly)
  out<-NULL
  for(i in seq_along(assembly)){
    tmp<-all_blast[all_blast$assembly==assembly[i],]
    gene<-unique(tmp$qseqid)
    for(j in seq_along(gene)){
      gen_len<-at_gen_length[at_gen_length$id==gene[j],2]
      max_len<-gen_len*1.25
      tmp1<-tmp[tmp$qseqid==gene[j],]
      chr<-unique(tmp1$sseqid)
      copy_counter <- 0
      for(k in seq_along(chr)){
        tmp2<-tmp1[tmp1$sseqid==chr[k],]
        cls<-gSoup::find_blast_clusters(tmp2[,c("sstart","send")],max_range = max_len)
        rng<-do.call(rbind,lapply(cls,range))
        tmp3<-data.frame(assembly=assembly[i],gene=gene[j],scaffold=chr[k],start=rng[,1],end=rng[,2],length=rng[,2]-rng[,1],max_length=max_len)
        tmp3<-tmp3[tmp3$length>(gen_len*.75),] # keep hits bigger than 75% of gene length
        if (nrow(tmp3) > 0) {
          tmp3$copy <- seq_len(nrow(tmp3)) + copy_counter
          copy_counter <- copy_counter + nrow(tmp3)  # Update the copy counter
          out <- rbind(out, tmp3)
        }
      }
    }
  }
  write.table(out,paste0("data/DefenseNstress/AT/blast2/length_sorted/",og,"length_sorted2.txt"),quote=F,sep="\t",row.names=F)
}

#### blast plot to visualize hits ------
bp_fl<-list.files("data/DefenseNstress/AT/blast/length_sorted/",full.names = T)
for(i in seq_along(bp_fl)){
  bp_tab<-read.table(bp_fl[i],h=T)
  pdf(paste0("data/DefenseNstress/AT/plots/",substr(basename(bp_fl[i]),1,10),"_blast_plot.pdf"),w=5,h=8)
  gSoup::blast_plot(bp_tab,cols=c("assembly","scaffold","start","end","gene"))
  dev.off()
}

#### find and remove overlapping hits -------
pb<-txtProgressBar(min=1,max=length(bp_fl),style=3,width = 50)
bp_fl<-list.files("data/DefenseNstress/AT/blast2/length_sorted/",full.names = T,pattern = "length_sorted2.txt")
for(k in seq_along(bp_fl)){
  setTxtProgressBar(pb,k)
  bp_tab<-read.table(bp_fl[k],h=T)
  assembly<-unique(bp_tab$assembly)
  
  tm0<-NULL
  for(g in seq_along(assembly)){
    df<-bp_tab[bp_tab$assembly==assembly[g],]
    if(nrow(df)>1){
      for (i in 1:(nrow(df))) {
        overlapping_rows<-c(i)
        if(i<nrow(df)){
          for (j in (i + 1):nrow(df)) {
            if (overlap_percentage(df$start[i], df$end[i], df$start[j], df$end[j]) > 0.7) {
              overlapping_rows<-c(unique(overlapping_rows),j) # hits that overlap more than 70%
            }
          }
        }
        for(rr in seq_along(overlapping_rows)){
          ll<-df[overlapping_rows[rr],2]
          tt<-df[which(df$gene==ll & df$copy==rr),]
          tm0<-rbind(tm0,tt)
        }
        tm0<-tm0[!duplicated(tm0[,4:5]),]
      }
    } else {
      tm0<-rbind(tm0,df)
    }
  }
  write.table(tm0,paste0("data/DefenseNstress/AT/cnv_within/ranges_dup_rem2/",gsub("length_sorted","pseu_dup_rem",basename(bp_fl[k]))),row.names = F,quote=F,sep="\t")
}



#### CNV within species ------------
#### list gene families with CNV ------------------
du_rm_fl<-list.files("data/DefenseNstress/AT/cnv_within/ranges_dup_rem2",full.names=T,pattern="pseu_dup_rem2.txt")
ogList<-read.csv("data/ogList_with_functions.csv")

out_cnv<-NULL
cnv_count<-NULL
for(i in seq_along(du_rm_fl)){
  og_cnv<-read.table(du_rm_fl[i],h=T)
  og_id<-substr(basename(du_rm_fl[i]),1,9)
  og_line<-ogList[ogList$Family.ID==og_id,]
  at_genes<-unique(stringr::str_trim(stringr::str_split_fixed(unlist(stringr::str_split(og_line[,"at_genes"],"\\,")),"\\.",n=2)[,1]))
  og_gene_count<-length(at_genes)
  assembly<-unique(og_cnv$assembly)
  out<-lapply(assembly,function(x){
    asX<-og_cnv[og_cnv$assembly==x,]
    gtab<-table(asX$gene)
    g_count<-length(gtab)
    gtab<-data.frame(gtab)
    gtab<-gtab[match(at_genes,gtab$Var1),]
    gtab$Var1<-at_genes
    gtab[is.na(gtab)]<-0
    gain<-sum(gtab$Freq>1)
    loss<-sum(gtab$Freq<1)
    gain_loss<-sum(gtab$Freq!=1)
    c(og_id,x,og_gene_count,g_count,gain,loss,gain_loss,sum(gtab$Freq))
  })
  out<-do.call(rbind,out)
  out2<-c(og_id,og_gene_count,sum(as.numeric(out[,7])!=0))
  cnv_count<-rbind(cnv_count,out2)
  out_cnv<-rbind(out_cnv,out)
}
colnames(out_cnv)<-c("family_id","assembly","family_size","gene_count","gain","loss","gain_loss_count","total_gene_count")
out_cnv<-as.data.frame(out_cnv)
write.csv(out_cnv,"data/DefenseNstress/AT/cnv_within/Athal_sig_fam_gain_loss2.2.csv",row.names = F)

## Assembly vs gene count per gene family 
library(tidyr)
library(dplyr)

out_cnv<-read.csv("data/DefenseNstress/AT/cnv_within/Athal_sig_fam_gain_loss2.2.csv")
matrix_out <- out_cnv %>%
  select(family_id, assembly, total_gene_count) %>%
  pivot_wider(names_from = assembly, values_from = total_gene_count)
og_vs_ass<-matrix_out
og_vs_ass[is.na(og_vs_ass)]<-0
write.table(og_vs_ass,"data/DefenseNstress/AT/cnv_within/geneFamily_vs_assembly_cnv2.2.txt",quote=F,row.names = F, sep="\t")


#### histograms of CNV per family -----------
out_cnv<-read.csv("data/DefenseNstress/AT/cnv_within/Athal_sig_fam_gain_loss2.2.csv")
ogs<-unique(out_cnv$family_id)
pdf("plots/AT_gene_counts_per_family_barplots2.2.pdf",h=5,w=8)
par(mfrow=c(5,3))
for(i in seq_along(ogs)){
  tm <- out_cnv[out_cnv$family_id == ogs[i],]
  tm$total_gene_count<-as.numeric(tm$total_gene_count)
  midpoints <- barplot(table(tm$total_gene_count), main = ogs[i], col = i, border = F)
  target_value <- as.numeric(tm$family_size[1])  # The actual numeric value
  closest_bar <- which(names(table(tm$total_gene_count)) == target_value)
  mtext(paste("Reference Family Size:", target_value), side = 3, line = 0.5, cex = 0.8)
  if (length(closest_bar) > 0) {
    lines(x=c(midpoints[closest_bar],midpoints[closest_bar]),y=c(0,max(table(tm$total_gene_count))), col="darkolivegreen",lwd=2)
  }
}
dev.off()

#### CNV per family table ----------------------
library(dplyr)
out_cnv<-read.csv("data/DefenseNstress/AT/cnv_within/Athal_sig_fam_gain_loss2.2.csv")
all_assemblies <- unique(out_cnv$assembly)
ogs_cnv <- out_cnv %>%
  group_by(family_id) %>%
  group_map(~ {
    df <- .x
    family <- .y$family_id
    ref_count <- unique(df$family_size)
    full_df <- data.frame(
      assembly = all_assemblies,
      gene_count = 0
    )
    full_df$gene_count[match(df$assembly, full_df$assembly)] <- df$gene_count
    
    # Calculate metrics
    median_count <- median(full_df$gene_count)
    mean_count <- mean(full_df$gene_count)
    sd_count <- sd(full_df$gene_count)
    cv <- ifelse(mean_count == 0, NA, sd_count / mean_count)
    dsp_ind<-var(full_df$gene_count)/mean_count
    cnv_percentage <- mean(full_df$gene_count != ref_count)
    
    tibble(
      family_id = family,
      refCount = ref_count,
      median_count = median_count,
      cnv_percentage = round(cnv_percentage, 3),
      CV = round(cv, 3),
      dispersion_index=dsp_ind
    )
  }) %>%
  bind_rows()
write.csv(ogs_cnv,"data/DefenseNstress/AT/cnv_within/cnv_percentage_per_family2.csv",row.names=F)


##### >>>>> MCA <<<<<<<< ########
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)

og_vs_ass<-read.table("data/DefenseNstress/AT/cnv_within/geneFamily_vs_assembly_cnv2.2.txt",h=T)
# Convert data frame to factor (cnv states)
pca_input <- as.data.frame(og_vs_ass)

nm <- gsub("[._]genomic","",sapply(colnames(pca_input)[-1],function(x){substr(x,17,nchar(x))}))
colnames(pca_input)[2:ncol(pca_input)] <- nm
rownames(pca_input) <- pca_input$family_id 
pca_input <- pca_input[, -1]  

# Ensure categorical variables are treated as factors
pca_matrix <- as.data.frame(t(pca_input))
pca_matrix[] <- lapply(pca_matrix, as.factor)

# Perform MCA 
mca_result <- MCA(pca_matrix, graph = F)

# Extract variance explained
explained_variance <- mca_result$eig[, 2]  # Percentage of variance per component
pc1_var <- round(explained_variance[1], 1)
pc2_var <- round(explained_variance[2], 1)

# Extract MCA coordinates
mca_df <- as.data.frame(mca_result$ind$coord)
mca_df$Assembly <- rownames(mca_df)
colnames(mca_df)<-c(paste0("Dim",1:5),"Assembly")

# Visualize MCA with variance in axis labels
pdf("plots/AT_MCAdim34_of_defNstres_families_136_assm_v1.pdf",h=8,w=12)
ggplot(mca_df, aes(x = Dim1, y = Dim2, label = Assembly)) +
  geom_point(size = 4, color = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  labs(title = "A.thaliana: MCA of LR-Assemblies Based on CNV in Def.Stress Gene Families",
       x = paste0("Dimension 1 (", pc1_var, "% variance)"),
       y = paste0("Dimension 2 (", pc2_var, "% variance)"))

dev.off()

pops<-data.frame(colnames(og_vs_ass))




### Ecology CNV per family vs  ------------------
# >>> IN THE SCRIPT "ecology.R" <<<<<




### TE ASSESSMENT -------------------------------------
# join all gene ranges
library(dplyr)
library(tidyr)
du_rm_fl<-list.files("data/DefenseNstress/AT/cnv_within/ranges_dup_rem2",full.names=T,pattern = "dup_rem2.txt")
tm_out<-NULL
for(j in seq_along(du_rm_fl)){
  tm2<-read.table(du_rm_fl[j],h=T)
  tm_out<-rbind(tm_out,tm2)
}
write.table(tm_out,"data/DefenseNstress/AT/cnv_within/all_DnS_blast_duprm2_joined.tsv",quote = F, row.names = F)


# search for TE classes in all genes from all assemblies
# Most frequent TEs associated with CNV
te_classes <- c("DNA/Helitron","DNA/MULE-MuDR","LTR/Copia", "LTR/Gypsy", "SINE/tRNA", "SINE/unknown","LINE/L1", "LINE/unknown","LTR/Ty3", "LTR/unknown")
tm_out<-read.table("data/DefenseNstress/AT/cnv_within/all_DnS_blast_duprm2_joined.tsv",h=T)
fls<-list.files("data/DefenseNstress/TE_profiling/RepeatMasker_results/AT",full.names=T,pattern = ".out")
ov_join<-NULL
sum_join<-NULL
pb<-txtProgressBar(max=length(fls),style=3,width=50)
for( i in seq_along(fls)){
  setTxtProgressBar(pb,i)
  tm<-data.table::fread(fls[i],skip=3,fill = 16)
  colnames(tm)[1:15] <- c("SW_score", "perc_div", "perc_del", "perc_ins", 
                          "query_sequence", "query_start", "query_end", "query_left", 
                          "strand", "repeat_name", "repeat_class", 
                          "repeat_start", "repeat_end", "repeat_left", "ID")
  tm<-data.frame(tm)
  te_df <- subset(tm, repeat_class %in% te_classes)
  
  assnam<-substr(basename(fls[i]),1,13)
  tm_ass<-tm_out[grep(assnam,tm_out$assembly),]
  
  tm_ass_buffered <- tm_ass %>%
    mutate(start_buffered = pmax(start - 1000, 1),
           end_buffered = end + 1000)
  
  te_df_clean <- te_df %>%
    rename(scaffold = query_sequence,
           te_start = query_start,
           te_end = query_end)
  
  overlap_df <- inner_join(
    tm_ass_buffered, te_df_clean, 
    by = "scaffold", 
    relationship = "many-to-many"
  ) %>%
    filter(te_start <= end_buffered & te_end >= start_buffered)
  
  summary_df <- overlap_df %>%
    group_by(gene) %>%
    summarise(n_TE_hits = n(),
              repeat_classes = paste(unique(repeat_class), collapse = ", "))
  
  ov_join<-rbind(ov_join,overlap_df)
  sum_join<-rbind(sum_join,data.frame(summary_df))
}
# export the output
write.csv(ov_join,"data/DefenseNstress/TE_profiling/AT_DnS_genes_TE_overlap.csv",row.names=F)
write.csv(sum_join,"data/DefenseNstress/TE_profiling/AT_DnS_genes_TE_overlap_summary.csv")


#### Visualize TE -------
library(ggplot2)
ov_join<-read.csv("data/DefenseNstress/TE_profiling/AT_DnS_genes_TE_overlap.csv")
total_assemblies <- 136
te_summary <- ov_join %>%
  group_by(gene, repeat_class) %>%
  summarise(assemblies = n_distinct(assembly), .groups = 'drop')


pdf("plots/AT_DnS_genes_TE_summary_barplot_v2.pdf",w=9,h=6)
ggplot(te_summary, aes(x = gene, y = assemblies, fill = repeat_class)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = .1,size = 7)) +
  labs(
    x = "Gene Family",
    y = "Number of assemblies with TE near gene members",
    fill = "Repeat class",
    title = "TE class presence in each gene family across assemblies"
  )
ggplot(te_summary, aes(x = gene, y = repeat_class, fill = assemblies)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "dodgerblue", high = 2) +
  theme_minimal() +
  labs(x = "Gene", y = "TE class", fill = "Assemblies", title = "TE presence heatmap") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = .1,size = 7) 
  )

# Summarize and sort data for correct slice order
repeat_summary <- ov_join %>%
  count(repeat_class) %>%
  arrange(desc(repeat_class)) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(round(percentage, 1), "%"),
    ypos = cumsum(percentage) - 0.5 * percentage
  )

# Donut chart with labels
ggplot(repeat_summary, aes(x = 2, y = percentage, fill = repeat_class)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  theme_void() +
  geom_text(aes(y = ypos, label = label), color = "black", size = 8) +
  labs(title = "Distribution of Repeat Classes near Genes", fill = "Repeat class")

dev.off()



### Extract genes from assemblies ---------
library(biomartr)
library(seqinr)
library(Biostrings)
genomes_path<-list.files("data/LR_assemblies/136_LR_fnas",full.names = T)
fam_files<-list.files("data/DefenseNstress/AT/cnv_within/ranges_dup_rem2",full.names=T,pattern="dup_rem2")
for(f in 1:length(fam_files)){ 
  og<-substr(basename(fam_files[f]),1,9)
  genes_sorted_path<-fam_files[f]
  gs_ranges<-read.table(genes_sorted_path,h=T)
  gs_ranges$assembly<-gsub("genomic","",gs_ranges$assembly)
  
  gene<-unique(gs_ranges$gene)
  
  genes<-list()
  for(j in seq_along(gene)){
    tmp<-gs_ranges[gs_ranges$gene==gene[j],]
    assembly<-unique(tmp$assembly)
    
    ass_list<-list()
    for( k in seq_along(assembly)){
      tmp1<-tmp[tmp$assembly==assembly[k],]
      scaffs<-unique(tmp1$scaffold)
      genome_file<-genomes_path[grep(assembly[k],genomes_path)]
      genome<-read_genome(genome_file)
      
      gene_list<-list()
      for(i in seq_along(scaffs)){
        tmp2<-tmp1[tmp1$scaffold==scaffs[i],]
        chr<-genome[[grep(scaffs[i],names(genome))]]
        dnastrin<-DNAStringSet()
        for(m in 1:nrow(tmp2)){
          dnastrin[[m]]<-chr[(tmp2$start[m]-200):(tmp2$end[m]+200)]
          names(dnastrin)[m]<-paste0(assembly[k],"_",gene[j],"_",scaffs[i],"_copy_",m)
        }
        gene_list[[i]]<-dnastrin
        print(dnastrin)
      }
      ass_list[[k]]<-gene_list
      print(assembly[k])
    }
    genes[[j]]<-ass_list
    
  }
  ggs<-unlist(genes)
  ggs<-do.call(c,ggs)
  nm<-sapply(names(ggs),function(x){substr(x,17,nchar(x))})
  write.fasta(as(ggs,"list"),names=nm,file.out = paste0("data/defensNstress/AT/sequences2/",og,"_all_genes_v2.2.fasta"))
}


## Population genetics summaries ---------------
#### HELIXER annotation of the extracted sequences ----------------------------
seq_dir="data/defensNstress/AT/sequences2"
out_dir="data/CNVcals/defensNstress/AT/annot2"

mkdir -p $out_dir
for file in ${seq_dir}/*v2.2.fasta; do
samp=`basename $file | sed 's/.fasta//'`
echo $samp
helixer --lineage land_plant --fasta-path $file --species Arabidopsis_thaliana --gff-output-path ${out_dir}/${samp}.gff3;
done

#### EXTRACTING helixer predicted CDS regions -------------------
## >>>>>>>> for separate genes <<<<<<<<<
library(Biostrings)
col.names<-c("seq","source","type","start","end","score","strand","phase","atributes")
fls<-list.files("data/defensNstress/AT/annot2",full.names = T,pattern="_v2.2.gff3")
seq_fl<-list.files("data/AT/sequences2",full.names=T,pattern="v2.2.fasta")

pb<-txtProgressBar(max=length(fls),style=3,width = 50)
for(k in seq_along(fls)){
  setTxtProgressBar(pb,k)
  og<-substr(basename(fls[k]),1,9)
  helix_out<-read.table(fls[k], col.names = col.names)
  # CDS extraction
  all_seq<-readDNAStringSet(seq_fl[grep(og,seq_fl)])
  genes <- unique(stringr::str_extract(helix_out$seq, "AT\\w{7}"))
  for(j in seq_along(genes)){
    tm<-helix_out[grep(genes[j],helix_out$seq),]
    seqs<-unique(tm$seq)
    #extract each sequence cds and the full gene
    full_genes<-list()
    seq_list<-list()
    for(i in seq_along(seqs)){
      ass<-seqs[i]
      tmp_seq<-all_seq[grep(ass,names(all_seq))]
      tmp_full<-tm[tm$seq==ass,]
      if(sum(tmp_full$type=="gene")>1){tmp_full<-tmp_full[1:which(tmp_full$type=="gene")[2]-1,]}
      gcoord<-tmp_full[tmp_full$type=="gene",c("start","end")]
      fg<-tmp_seq[[1]][gcoord$start:gcoord$end]
      ##---
      if (any(tmp_full$strand == "-")) {
        cds_rows <- tmp_full[tmp_full$type == "CDS", ]
        cds_rows <- cds_rows[order(cds_rows$start), ]
        cds_seq_list <- mapply(function(start, end) {
          subseq(tmp_seq, start = start, end = end)
        }, cds_rows$start, cds_rows$end, SIMPLIFY = FALSE)
        cds_concat <- do.call(xscat, cds_seq_list)
        cds_final <- reverseComplement(cds_concat)
        seq_list[[i]]<-cds_final[[1]]
        #reverse full gene 
        fg<-reverseComplement(fg)
      } else {
        tmp<-tmp_full[tmp_full$type=="CDS",]
        ##---
        seqq<-NULL
        for(l in 1:nrow(tmp)){
          coords<-unlist(tmp[l,4:5])
          tseqq<-tmp_seq[[1]][coords[1]:coords[2]]
          if(l==2){seqq<-c(seqq[[1]],tseqq)}else{seqq<-c(seqq,tseqq)}
        }
        seq_list[[i]]<-seqq
      }
      full_genes[[i]]<-fg
    }
    #save fasta
    ass_nm<-seqs
    names(seq_list)<-ass_nm
    all_cds_aug<-DNAStringSet(unlist(seq_list))
    all_cds_aug<-as(all_cds_aug,"list")
    fam_dir<-"data/defensNstress/AT/cds/sep_genes2"
    if(!dir.exists(fam_dir)){dir.create(fam_dir)}
    seqinr::write.fasta(all_cds_aug,names = ass_nm,file.out=paste0(fam_dir,"/","cds_v5_",og,"_",genes[j],".fasta"))
    full_dir="data/defensNstress/AT/sequences2/aug_ext"
    all_fullg<-DNAStringSet(unlist(full_genes))
    all_fullg<-as(all_fullg,"list")
    seqinr::write.fasta(all_fullg,names = ass_nm,file.out=paste0(full_dir,"/","full_gene_",og,"_",genes[j],"_v5.fasta"))
    
  }
}


### Find truncated/pseudo genes ------------------------
library(Biostrings)
fls<-list.files("data/defensNstress/AT/cds/sep_genes2",full.names=T,pattern="cds_v5")
full_fls<-list.files("data/defensNstress/AT/sequences2/aug_ext",full.names=T,pattern="_v5.fasta")
out_dir<-"data/defensNstress/AT/cds/sep_genes2/fixed/"
out_dfull<-"data/defensNstress/AT/sequences2/aug_ext/fixed/"
if(!dir.exists(out_dir)){dir.create(out_dir)}
if(!dir.exists(out_dfull)){dir.create(out_dfull)}

ll<-lapply(fls,function(x){
  tm<-readDNAStringSet(x)
  lth<-width(tm)
  lmed<-median(lth)
  upl<-lmed*1.2
  lowl<-lmed*0.8
  pseuN<-which(lth < lowl | lth > upl)
  pseu<-names(tm)[pseuN]
  if(length(pseu)<1){pseu<-NA}
  og<-substr(basename(x),8,16)
  gene<-substr(basename(x),18,26)
  out<-cbind(og,gene,pseu)
  
  if(length(pseu)<2){if(!is.na(pseu)){tm<-tm[-pseuN]}}
  if(length(tm)<3){print(paste0(basename(x), " has less than 3 seqs"))}
  seqinr::write.fasta(as(tm,"list"),names=names(tm),file.out=paste0(out_dir,gsub(".fasta","_fixed.fasta",basename(x))))
  
  fulgfile<-full_fls[grep(substr(basename(x),8,26),full_fls)] 
  ful<-readDNAStringSet(fulgfile)
  if(length(pseu)<2){if(!is.na(pseu))ful<-ful[-pseuN]}
  seqinr::write.fasta(as(ful,"list"),names=names(ful),file.out=paste0(out_dfull,gsub(".fasta","_fixed.fasta",basename(fulgfile))))
  return(out)
})
ll<-data.frame(do.call(rbind,ll))
write.csv(ll,"data/defensNstress/AT/stats/AT_putative_pseudo_genes_based_on_cds_length_v5.csv",row.names=F)



### MAFFT ALIGNMENT OF AUGUSTUS PREDICTED FULL GENES ---------
# in fixed folders you find the fixed sequence by removing "pseudo" genes
# "pseudo" genes were determined by checking the cds lengths
input_dir="data/defensNstress/AT/sequences2/aug_ext/fixed"
output_dir="data/defensNstress/AT/align2/full_genes/sep_genes_mafft"

mkdir -p "$output_dir"
for file in "$input_dir"/*.fasta; do
base=$(basename "$file" .fasta)
mafft --auto --adjustdirection --anysymbol "$file" > "$output_dir/${base}_aligned.fasta"
echo "Full gene alignment done for $base"
done

### MACSE ALIGNMENT OF GENE CDS ---------------
# in fixed folders you find the fixed sequence by removing "pseudo" genes
# "pseudo" genes were determined by checking the cds lengths
input_dir="data/defensNstress/AT/cds/sep_genes2/fixed"
output_dir="data/defensNstress/AT/align2/cds/sep_genes_macse"
macse_jar="macse_v2.07.jar"

mkdir -p "$output_dir"
for file in "$input_dir"/*.fasta; do
base=$(basename "$file" .fasta)
java -jar "$macse_jar" -prog alignSequences \
-seq "$file" \
-out_NT "$output_dir/${base}_aligned_NT.fasta" \
-out_AA "$output_dir/${base}_aligned_AA.fasta" 
echo "Codon alignment done for $base"
done



## STATS ----------------------
### Pseudogene table -------------------------------
pse<-read.csv("data/DefenseNstress/AT/stats/AT_putative_pseudo_genes_based_on_cds_length_v5.csv")
library(dplyr)
summary_df <- pse %>%
  filter(!is.na(pseu)) %>%
  group_by(og) %>%
  summarise(
    genes = n_distinct(gene),
    freq_gene = gene %>% 
      table() %>%
      which.max() %>%
      names()
  )
all_ogs <- unique(pse$og)
summary_df <- tibble(og = all_ogs) %>%
  left_join(summary_df, by = "og") %>%
  mutate(
    genes = ifelse(is.na(genes), 0, genes),
    freq_gene = ifelse(is.na(freq_gene), NA, freq_gene)
  )
summary_df<-data.frame(summary_df)
# combine with cnv percentages 
cnv_perc<-read.csv("data/DefenseNstress/AT/cnv_within/cnv_percentage_per_family2.csv",h=T)
summary_df<-summary_df[match(cnv_perc$family,summary_df$og),]
summary_df_joined<-data.frame(cnv_perc,summary_df[,-1])
colnames(summary_df_joined)[6:7]<-c("pseudo_genes","frequent_pseudogene")
summary_df_joined$pse_perc<-summary_df_joined$pseudo_genes/summary_df_joined$refCount
write.csv(summary_df_joined,"data/DefenseNstress/AT/stats/AT_putative_pseudo_genes_counts_with_cnv_v5.csv",row.names = F)




### Pi / Theta --------------------------------------
#calculate pi/theta from mafft alignment of full genes
library(gSoup)
library(Biostrings)
algn_files<-list.files("data/defensNstress/AT/align2/full_genes/sep_genes_mafft",full.names=T,pattern=".fasta")

pb<-txtProgressBar(max=length(algn_files),style=3,width=50)
div_table<-NULL
for(i in seq_along(algn_files)){
  setTxtProgressBar(pb,i)
  algn_seq<-readDNAStringSet(algn_files[i])
  og<-substr(basename(algn_files[i]),11,19)
  gene<-substr(basename(algn_files[i]),21,29)
  if(length(algn_seq)<3){per_gene_dstat<-c(og,gene,rep(NA,7))
  div_table<-rbind(div_table,per_gene_dstat)
  next}
  dstats<-div.stats(algn_seq)
  per_gene_dstat<-c(og,gene,length(algn_seq),mean(width(algn_seq)),unlist(dstats[-5]))
  div_table<-rbind(div_table,per_gene_dstat)
}
div_table<-data.frame(div_table)
colnames(div_table)[1:4]<-c("gene_family","gene","no_seqs","mean_algn_len")
write.csv(div_table,"data/defensNstress/AT/stats/AT_pi_theta_table_v5.csv",row.names=F)

#### plot pi / theta --------------------------
# boxplots 
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 

pi_theta <- read.csv("data/DefenseNstress/AT/stats/AT_pi_theta_table_v5.csv")
#  data
plot_data <- pi_theta %>%
  filter(!is.na(pi), !is.na(theta)) %>%
  pivot_longer(
    cols = c(pi, theta),
    names_to = "statistic", 
    values_to = "value"
  ) %>%
  mutate(statistic = factor(statistic, levels = c("pi", "theta")),
         statistic_label = ifelse(statistic == "pi", "π", "θ"))

# Main boxplot with points
main_plot <- ggplot(plot_data, aes(x = gene_family, y = value)) +
  geom_boxplot(
    aes(fill = statistic),
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = statistic),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.5,
    alpha = 0.6
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("π", "θ")) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("π", "θ")) +
  labs(x = "Gene Family", y = "Diversity Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

# Density plot
density_plot <- ggplot(plot_data, aes(x = value, fill = statistic)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(x = "Diversity Value", y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Combine plots
combined_plot <- main_plot + density_plot + 
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "π and θ Distributions by Gene Family",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  )
print(combined_plot)
# save pdf in 10 to 3.8 inch ratio


### Pi A/S ------------------------------------
# remove "!" from alignments before reading into R << with bash
# MAC
# sed -i.bak '' 's/!/N/g' *.fasta
# # LINUX
# sed -i.bak 's/!/N/g' *.fasta

library(Biostrings)
library(gSoup)
algn_files<-list.files("data/defensNstress/AT/align2/cds/sep_genes_macse",full.names=T,pattern=".fasta")
pb<-txtProgressBar(max=length(algn_files),style=3,width=50)
div_table<-NULL
for(i in seq_along(algn_files)){
  setTxtProgressBar(pb,i)
  algn_seq<-readDNAStringSet(algn_files[i])
  og<-substr(basename(algn_files[i]),9,16)
  gene<-substr(basename(algn_files[i]),18,26)
  if(length(algn_seq)<3){per_gene_dstat<-c(og,gene,rep(NA,7))
  div_table<-rbind(div_table,per_gene_dstat)
  next}
  dstats<-div.stats(algn_seq,pi_AS = T)
  per_gene_dstat<-c(og,gene,length(algn_seq),mean(width(algn_seq)),unlist(dstats[-5]))
  div_table<-rbind(div_table,per_gene_dstat)
}
div_table<-data.frame(div_table)
colnames(div_table)[1:4]<-c("gene_family","gene","no_seqs","mean_algn_len")
write.csv(div_table,"data/defensNstress/AT/stats/AT_pi_AS_table_v5.csv",row.names=F)

#### Pi A/S boxplots ------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 
library(ggside) 

pias_tab<-read.csv("data/DefenseNstress/AT/stats/AT_pi_AS_table_v5.csv")

# Calculate log2(piA/piS) ratio
pias_ratio <- pias_tab %>%
  filter(!is.na(pi_AS.corrected_pi_syn), 
         !is.na(pi_AS.corrected_pi_nonsyn),
         pi_AS.corrected_pi_syn > 0) %>%
  mutate(log2_ratio = log2(pi_AS.corrected_pi_nonsyn / pi_AS.corrected_pi_syn))

# Create the plot density
ggplot(pias_ratio, aes(x = gene_family, y = log2_ratio)) +
  # Main plot elements
  geom_boxplot(fill = "#4E79A7", alpha = 0.7, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "#4E79A7") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  
  # Right-side density plot (single combined density)
  geom_ysidedensity(
    aes(x = after_stat(density)),  
    fill = "#4E79A7", 
    alpha = 0.3,
    position = "identity"
  ) +
  
  # scaling and theming
  scale_ysidex_continuous(expand = expansion(mult = c(0, 0.1)), breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  
  # Labels and theme
  labs(
    x = "Gene Family",
    y = "log2(πN/πS)",
    title = "Nonsynonymous to Synonymous Diversity Ratio",
    subtitle = "Positive values indicate higher nonsynonymous diversity (πN > πS)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    ggside.panel.scale = 0.2
  )

# save pdf in 10 x 3.8 inch ratio




# LYRATA -------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Extract gene sequences and save in FASTA 
al_ids<-def_gene_ids[,c(1,4)]
og_ids<-def_gene_ids$Orthogroup
# genome and the gff files
al_gff <- biomartr::read_gff("data/Alyrata_384_v2.1.gene_exons.gff3")
al_genome <- readDNAStringSet("data/Alyrata_384_v1.fa")

for(j in 1:nrow(al_ids)){
  newIds<-stringr::str_split_fixed(unlist(stringr::str_split(al_ids[j,2],pattern="\\,")),pattern="\\.",n=2)[,1]
  gene_list<-list()
  for(i in seq_along(newIds)){
    g_match<-al_gff[grep(gsub("\\ ","",newIds[i]),al_gff$attribute),]
    rr<-which(g_match$type=="gene")[1]
    coords<-range(g_match[rr,4:5])
    chr<-as.character(g_match[1,1])
    gene<-al_genome[[which(names(al_genome)==chr)]][min(coords):max(coords)]
    gene_list[[i]]<-gene
  }
  seqinr::write.fasta(gene_list,names = newIds,file.out=paste0("data/DefenseNstress/AL/refseq/",unlist(al_ids[j,1]),".fasta"))
}


### find gene lengths --------
fls<-list.files("data/DefenseNstress/AL/refseq",full.names=T)
al_gen_length<-do.call(rbind,lapply(fls,function(x){dd<-readDNAStringSet(x);return(cbind(names(dd),width(dd),gsub(".fasta","",basename(x))))}))
al_gen_length<-data.frame(id=al_gen_length[,1],length=al_gen_length[,2],og=al_gen_length[,3])
al_gen_length$length<-as.numeric(al_gen_length$length)
al_gen_length$id<-as.character(stringr::str_trim(al_gen_length$id))
write.table(al_gen_length,"data/DefenseNstress/AL/all_gene_lengths.txt",row.names = F,quote=F,sep="\t")


## BLAST FOR GENES --------------------------------
reference_dir="data/DefensNstress/AL/refseqs"
seq_dir="data/all_assmbl"
out_dir="DefensNstress/AL/blast/"

for ref_fa in ${reference_dir}/*.fasta; do
ref=`basename $ref_fa | sed 's/.fasta//'`
echo $ref
for file in ${seq_dir}/*.fa; do
samp=`basename $file | sed 's/.fa//'`
echo $samp
blastn -query $ref_fa -subject $file -outfmt '7 std stitle' -out ${out_dir}${ref}${samp}.blastn ;
done;
done


## CHEKC BLAST HITS  --------------------------
b_files<-list.files("data/DefenseNstress/AL/blast/out",full.names = T)
og_ids<-al_ids$Orthogroup
for(j in seq_along(og_ids)){
  og_files<-b_files[grep(og_ids[j],b_files)]
  og_added<-NULL
  for(i in seq_along(og_files)){
    b_hits<-read.table(og_files[i],comment.char = "#",
                       col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","stitle"))
    b_hits$stitle<-stringr::str_split_fixed(gsub(".blastn","",basename(og_files[i])),
                                            pattern="_",n=2)[,2]
    og_added<-rbind(og_added,b_hits)
  }
  write.table(og_added,paste0("data/DefenseNstress/AL/blast/og_sorted/",og_ids[j],"_bhits_v1.txt"),quote=F,row.names=F,sep="\t")
}
gSoup::blast_plot(b_hits,cols=c("stitle","sseqid","sstart","send","qseqid"))


## GROUP THE BLAST HITS  -----------------------------
al_gen_length<-read.table("data/DefenseNstress/AL/all_gene_lengths.txt",h=T)
og_sorted<-list.files("data/DefenseNstress/AL/blast/og_sorted/",full.names = T)
pb<-txtProgressBar(min=0,max=length(og_sorted),style = 3,width = 50)
for(m in seq_along(og_sorted)){
  setTxtProgressBar(pb,m)
  all_blast<-read.table(og_sorted[m],header = T)
  all_blast<-all_blast[all_blast$pident>=90,]
  og<-gsub("bhits_v1.txt","",basename(og_sorted[m]))
  all_blast$assembly<-stringr::str_split_fixed(all_blast$stitle, ".ragtag",n=2)[,1]
  all_blast$assembly<-stringr::str_split_fixed(all_blast$assembly, "_fixed",n=2)[,1]
  assembly<-unique(all_blast$assembly)
  out<-NULL
  for(i in seq_along(assembly)){
    tmp<-all_blast[all_blast$assembly==assembly[i],]
    gene<-unique(tmp$qseqid)
    for(j in seq_along(gene)){
      gen_len<-al_gen_length[al_gen_length$id==gene[j],2]
      max_len<-gen_len*1.25
      tmp1<-tmp[tmp$qseqid==gene[j],]
      chr<-unique(tmp1$sseqid)
      copy_counter <- 0
      for(k in seq_along(chr)){
        tmp2<-tmp1[tmp1$sseqid==chr[k],]
        cls<-gSoup::find_blast_clusters(tmp2[,c("sstart","send")],max_range = max_len)
        rng<-do.call(rbind,lapply(cls,range))
        tmp3<-data.frame(assembly=assembly[i],gene=gene[j],scaffold=chr[k],start=rng[,1],end=rng[,2],length=rng[,2]-rng[,1],max_length=max_len)
        tmp3<-tmp3[tmp3$length>(gen_len*.75),]
        
        if (nrow(tmp3) > 0) {
          tmp3$copy <- seq_len(nrow(tmp3)) + copy_counter
          copy_counter <- copy_counter + nrow(tmp3)  
          out <- rbind(out, tmp3)
        }
      }
    }
  }
  write.table(out,paste0("data/DefenseNstress/AL/blast/length_sorted/",og,"length_sorted_v5.txt"),quote=F,sep="\t",row.names=F)
}


#### CNV within species ------------
bp_fl<-list.files("data/DefenseNstress/AL/blast/length_sorted/",full.names = T,pattern = "_v5.txt")
pb<-txtProgressBar(min=1,max=length(bp_fl),style=3,width = 50)
for(k in seq_along(bp_fl)){
  setTxtProgressBar(pb,k)
  bp_tab<-read.table(bp_fl[k],h=T)
  assembly<-unique(bp_tab$assembly)
  tm0<-NULL
  for(g in seq_along(assembly)){
    df<-bp_tab[bp_tab$assembly==assembly[g],]
    if(nrow(df)>1){
      for (i in 1:(nrow(df))) {
        overlapping_rows<-c(i)
        if(i<nrow(df)){
          for (j in (i + 1):nrow(df)) {
            if (gSoup::overlap_percentage(df$start[i], df$end[i], df$start[j], df$end[j]) > 0.7) {
              overlapping_rows<-c(unique(overlapping_rows),j)
              
            }
          }
        }
        for(rr in seq_along(overlapping_rows)){
          ll<-df[overlapping_rows[rr],2]
          tt<-df[which(df$gene==ll & df$copy==rr),]
          tm0<-rbind(tm0,tt)
        }
        tm0<-tm0[!duplicated(tm0[,4:5]),]
      }
    } else {
      tm0<-rbind(tm0,df)
    }
  }
  write.table(tm0,paste0("data/DefenseNstress/AL/cnv_within/dup_remd/",gsub("length_sorted","pseu_dup_rem",basename(bp_fl[k]))),row.names = F,quote=F,sep="\t")
}

#### list gene families with CNV per species ------------------
du_rm_fl<-list.files("data/DefenseNstress/AL/cnv_within/dup_remd",full.names=T, pattern="_v5.txt")
ogList<-read.csv("data/ogList_with_functions.csv")
def_gene_ids<-data.table::fread("data/defense_and_stress_fast_evol_families_gene_ids_4_species.txt",h=T)
out_cnv<-NULL
cnv_count<-NULL
for(i in seq_along(du_rm_fl)){
  og_cnv<-read.table(du_rm_fl[i],h=T)
  og_id<-substr(basename(du_rm_fl[i]),1,9)
  og_line<-def_gene_ids[def_gene_ids$Orthogroup==og_id,]
  al_genes<-unique(stringr::str_trim(stringr::str_split_fixed(unlist(stringr::str_split(og_line[,4],"\\,")),"\\.",n=2)[,1]))
  og_gene_count<-length(al_genes)
  assembly<-unique(og_cnv$assembly)
  # x<-assembly[2]
  out<-lapply(assembly,function(x){
    asX<-og_cnv[og_cnv$assembly==x,]
    gtab<-table(asX$gene)
    g_count<-length(gtab)
    gtab<-data.frame(gtab)
    gtab<-gtab[match(al_genes,gtab$Var1),]
    gtab$Var1<-al_genes
    gtab[is.na(gtab)]<-0
    gain<-sum(gtab$Freq>1)
    loss<-sum(gtab$Freq<1)
    gain_loss<-sum(gtab$Freq!=1)
    c(og_id,x,og_gene_count,g_count,gain,loss,gain_loss,sum(gtab$Freq))
  })
  out<-do.call(rbind,out)
  out2<-c(og_id,og_gene_count,sum(as.numeric(out[,7])!=0))
  cnv_count<-rbind(cnv_count,out2)
  out_cnv<-rbind(out_cnv,out)
}
colnames(out_cnv)<-c("family_id","assembly","family_size","gene_count","gain","loss","gain_loss_count","total_gene_count")
out_cnv<-as.data.frame(out_cnv)
write.csv(out_cnv,"data/DefenseNstress/AL/cnv_within/Alyrata_sig_fam_gain_loss_v5.csv",row.names = F)



#### histograms of CNV per family -----------
out_cnv<-read.csv("data/DefenseNstress/AL/cnv_within/Alyrata_sig_fam_gain_loss_v5.csv")
ogs<-unique(out_cnv$family_id)
pdf("plots/Alyr_gene_counts_per_family_barplots_v5.2.pdf",h=5,w=8)
par(mfrow=c(5,3))
for(i in seq_along(ogs)){
  tm <- out_cnv[out_cnv$family_id == ogs[i],]
  if(nrow(tm)<1){next}
  tm$total_gene_count<-as.numeric(tm$total_gene_count)
  midpoints <- barplot(table(tm$total_gene_count), main = ogs[i], col = i, border = F)
  target_value <- as.numeric(tm$family_size[1])
  closest_bar <- which(names(table(tm$total_gene_count)) == target_value)
  mtext(paste("Reference Family Size:", target_value), side = 3, line = 0.5, cex = 0.8)
  if (length(closest_bar) > 0) {
    lines(x=c(midpoints[closest_bar],midpoints[closest_bar]),y=c(0,max(table(tm$total_gene_count))), col="darkolivegreen",lwd=2)
  }
}
dev.off()


#### CNV per family table ----------------------
library(dplyr)
out_cnv<-read.csv("data/DefenseNstress/AL/cnv_within/Alyrata_sig_fam_gain_loss_v5.csv")
all_assemblies <- unique(out_cnv$assembly)
ogs_cnv <- out_cnv %>%
  group_by(family_id) %>%
  group_map(~ {
    df <- .x
    family <- .y$family_id
    ref_count <- unique(df$family_size)
    full_df <- data.frame(
      assembly = all_assemblies,
      gene_count = 0
    )
    full_df$gene_count[match(df$assembly, full_df$assembly)] <- df$gene_count
    
    # Calculate metrics
    median_count <- median(full_df$gene_count)
    mean_count <- mean(full_df$gene_count)
    sd_count <- sd(full_df$gene_count)
    cv <- ifelse(mean_count == 0, NA, sd_count / mean_count)
    dsp_ind<-var(full_df$gene_count)/mean_count
    cnv_percentage <- mean(full_df$gene_count != ref_count)
    
    tibble(
      family_id = family,
      refCount = ref_count,
      median_count = median_count,
      cnv_percentage = round(cnv_percentage, 3),
      CV = round(cv, 3),
      dispersion_index=dsp_ind
    )
  }) %>%
  bind_rows()
write.csv(ogs_cnv,"data/DefenseNstress/AL/cnv_within/cnv_percentage_per_family_v5.csv",row.names=F)



### CNV per family vs Ecology ------------------
# ** Run multiple regression on cnv vs env on the def & stress gene families ***
# in the R script ecology.R



### TE ASSESSMENT -------------------------------------
# Find if fast evolving gene family members have TE or located near TE
library(dplyr)
library(tidyr)

te_classes <- c("DNA/Helitron","DNA/MULE-MuDR","LTR/Copia", "LTR/Gypsy", "SINE/tRNA", "SINE/unknown","LINE/L1", "LINE/unknown","LTR/Ty3", "LTR/unknown")

du_rm_fl<-list.files("data/DefenseNstress/AL/cnv_within/dup_remd",full.names=T, pattern="_v5.txt")
tm_out<-NULL
for(j in seq_along(du_rm_fl)){
  tm2<-read.table(du_rm_fl[j],h=T)
  tm_out<-rbind(tm_out,tm2)
}
write.table(tm_out,"data/DefenseNstress/AL/cnv_within/all_DnS_blast_duprm_v5_joined.tsv",quote = F, row.names = F)

tm_out<-read.table("data/DefenseNstress/AL/cnv_within/all_DnS_blast_duprm_v5_joined.tsv",h=T)
fls<-list.files("data/DefenseNstress/TE_profiling/RepeatMasker_results/AL",full.names=T, pattern = ".fa.out")

ov_join<-NULL
sum_join<-NULL
pb<-txtProgressBar(max=length(fls),style=3,width=50)
for( i in seq_along(fls)){
  setTxtProgressBar(pb,i)
  tm<-data.table::fread(fls[i],skip=3,fill = 16)
  colnames(tm)[1:15] <- c("SW_score", "perc_div", "perc_del", "perc_ins", 
                          "query_sequence", "query_start", "query_end", "query_left", 
                          "strand", "repeat_name", "repeat_class", 
                          "repeat_start", "repeat_end", "repeat_left", "ID")
  tm<-data.frame(tm)
  te_df <- subset(tm, repeat_class %in% te_classes)
  
  assnam<-gsub("-","_",gsub(".fa.out","",basename(fls[i])))
  tm_ass<-tm_out[tm_out$assembly==assnam,]
  
  tm_ass_buffered <- tm_ass %>%
    mutate(start_buffered = pmax(start - 1000, 1),
           end_buffered = end + 1000)
  
  te_df_clean <- te_df %>%
    rename(scaffold = query_sequence,
           te_start = query_start,
           te_end = query_end)
  
  overlap_df <- inner_join(
    tm_ass_buffered, te_df_clean, 
    by = "scaffold", 
    relationship = "many-to-many"
  ) %>%
    filter(te_start <= end_buffered & te_end >= start_buffered)
  
  summary_df <- overlap_df %>%
    group_by(gene) %>%
    summarise(n_TE_hits = n(),
              repeat_classes = paste(unique(repeat_class), collapse = ", "))
  
  ov_join<-rbind(ov_join,overlap_df)
  sum_join<-rbind(sum_join,data.frame(summary_df))
}

write.csv(ov_join,"data/DefenseNstress/TE_profiling/AL_DnS_genes_TE_overlap_v5.csv",row.names=F)
write.csv(sum_join,"data/DefenseNstress/TE_profiling/AL_DnS_genes_TE_overlap_summary_v5.csv")


#### Visualize TE -------
library(ggplot2)
ov_join<-read.csv("data/DefenseNstress/TE_profiling/AL_DnS_genes_TE_overlap_v5.csv")
total_assemblies <- 21
te_summary <- ov_join %>%
  group_by(gene, repeat_class) %>%
  summarise(assemblies = n_distinct(assembly), .groups = 'drop')

pdf("plots/AL_DnS_genes_TE_summary_barplot_v5_v2.pdf",w=9,h=6)
ggplot(te_summary, aes(x = gene, y = assemblies, fill = repeat_class)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = .1,size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  labs(
    x = "Gene",
    y = "Number of assemblies with TE near gene",
    fill = "Repeat class",
    title = "TE class presence near each gene across assemblies"
  )
ggplot(te_summary, aes(x = gene, y = repeat_class, fill = assemblies)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "dodgerblue", high = 2) +
  theme_minimal() +
  labs(x = "Gene", y = "TE class", fill = "Assemblies", title = "TE presence heatmap") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = .1,size = 7) 
  )

# Summarize and sort data for correct slice order
repeat_summary <- ov_join %>%
  count(repeat_class) %>%
  arrange(desc(repeat_class)) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(round(percentage, 1), "%"),
    ypos = cumsum(percentage) - 0.5 * percentage
  )

# Donut chart with labels
ggplot(repeat_summary, aes(x = 2, y = percentage, fill = repeat_class)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  theme_void() +
  geom_text(aes(y = ypos, label = label), color = "black", size = 8) +
  labs(title = "Distribution of Repeat Classes near Genes", fill = "Repeat class")

dev.off()




#### Extract genes from assemblies ---------
library(biomartr)
library(seqinr)
library(Biostrings)
genomes_path<-list.files("data/all_assmbl",full.names = T)
fam_files<-list.files("data/defensNstress/AL/cnv_within",full.names=T,pattern = "dup_rem_v5.txt")

for(f in 1:length(fam_files)){
  og<-substr(basename(fam_files[f]),1,9)
  genes_sorted_path<-fam_files[f]
  gs_ranges<-read.table(genes_sorted_path,h=T)
  gene<-unique(gs_ranges$gene)
  
  genes<-list()
  for(j in seq_along(gene)){
    tmp<-gs_ranges[gs_ranges$gene==gene[j],]
    assembly<-unique(tmp$assembly)
    
    ass_list<-list()
    for( k in seq_along(assembly)){
      tmp1<-tmp[tmp$assembly==assembly[k],]
      scaffs<-unique(tmp1$scaffold)
      genome_file<-genomes_path[grep(assembly[k],genomes_path)]
      genome<-read_genome(genome_file)
      
      gene_list<-list()
      for(i in seq_along(scaffs)){
        tmp2<-tmp1[tmp1$scaffold==scaffs[i],]
        chr<-genome[[grep(scaffs[i],names(genome))]]
        dnastrin<-DNAStringSet()
        for(m in 1:nrow(tmp2)){
          dnastrin[[m]]<-chr[(tmp2$start[m]-200):(tmp2$end[m]+200)]
          names(dnastrin)[m]<-paste0(gene[j],"_",assembly[k],"_",substr(scaffs[i],1,10),"_copy_",m)
        }
        gene_list[[i]]<-dnastrin
        print(dnastrin)
      }
      ass_list[[k]]<-gene_list
      print(assembly[k])
    }
    genes[[j]]<-ass_list
  }
  
  ggs<-unlist(genes)
  ggs<-do.call(c,ggs)
  write.fasta(as(ggs,"list"),names=names(ggs),file.out = paste0("DefensNstress/AL/sequences/",og,"_all_genes_v5.fasta"))
}


#### HELIXER annotation of the extracted sequences ----------------------------
seq_dir="DefensNstress/AL/sequences"
out_dir="DefensNstress/AL/annot"

mkdir -p $out_dir
for file in ${seq_dir}/*v5.fasta; do
samp=`basename $file | sed 's/.fasta//'`
echo $samp
helixer --lineage land_plant --fasta-path $file --species Arabidopsis_lyrata --gff-output-path ${out_dir}/${samp}.gff3;
done

#### EXTRACTING helixer predicted CDS regions -------------------
## >>>>>>>> for separate genes <<<<<<<<<
library(Biostrings)
col.names<-c("seq","source","type","start","end","score","strand","phase","atributes")
fls<-list.files("DefensNstress/AL/annot/AUG",full.names = T,pattern="v5.gff3")
seq_fl<-list.files("DefensNstress/AL/sequences",full.names=T,pattern="v5.fasta")

pb<-txtProgressBar(max=length(fls),style=3,width = 50)
for(k in seq_along(fls)){
  setTxtProgressBar(pb,k)
  og<-substr(basename(fls[k]),1,9)
  helix_out<-biomartr::read_gff(fls[k])
  colnames(helix_out)[1]<-"seq"
  # CDS extraction
  all_seq<-readDNAStringSet(seq_fl[grep(og,basename(seq_fl))])
  genes <- unique(substr(names(all_seq),1,9))
  for(j in seq_along(genes)){
    tm<-helix_out[grep(genes[j],helix_out$seq),]
    seqs<-unique(tm$seq)
    #extract each sequence cds and the full gene
    full_genes<-list()
    seq_list<-list()
    for(i in seq_along(seqs)){
      ass<-seqs[i]
      tmp_seq<-all_seq[grep(ass,names(all_seq))]
      tmp_full<-tm[tm$seq==ass,]
      if(sum(tmp_full$type=="gene")>1){tmp_full<-tmp_full[1:which(tmp_full$type=="gene")[2]-1,]}
      gcoord<-tmp_full[tmp_full$type=="gene",c("start","end")]
      fg<-tmp_seq[[1]][gcoord$start:gcoord$end]
      ##---
      if (any(tmp_full$strand == "-")) {
        cds_rows <- tmp_full[tmp_full$type == "CDS", ]
        cds_rows <- cds_rows[order(cds_rows$start), ]
        cds_seq_list <- mapply(function(start, end) {
          subseq(tmp_seq, start = start, end = end)
        }, cds_rows$start, cds_rows$end, SIMPLIFY = FALSE)
        cds_concat <- do.call(xscat, cds_seq_list)
        cds_final <- reverseComplement(cds_concat)
        seq_list[[i]]<-cds_final[[1]]
        #reverse full gene 
        fg<-reverseComplement(fg)
      } else {
        tmp<-tmp_full[tmp_full$type=="CDS",]
        ##---
        seqq<-NULL
        for(l in 1:nrow(tmp)){
          coords<-unlist(tmp[l,4:5])
          tseqq<-tmp_seq[[1]][coords[1]:coords[2]]
          if(l==2){seqq<-c(seqq[[1]],tseqq)}else{seqq<-c(seqq,tseqq)}
        }
        seq_list[[i]]<-seqq
      }
      full_genes[[i]]<-fg
    }
    #save fasta
    ass_nm<-seqs
    names(seq_list)<-ass_nm
    all_cds_aug<-DNAStringSet(unlist(seq_list))
    all_cds_aug<-as(all_cds_aug,"list")
    # cds
    fam_dir<-"DefensNstress/AL/cds/sep_genes/v5.1"
    if(!dir.exists(fam_dir)){dir.create(fam_dir)}
    seqinr::write.fasta(all_cds_aug,names = ass_nm,file.out=paste0(fam_dir,"/","cds_v5_",og,"_",genes[j],".fasta"))
    #full gene
    full_dir="DefensNstress/AL/sequences/aug_ext/v5.1"
    if(!dir.exists(full_dir)){dir.create(full_dir)}
    all_fullg<-DNAStringSet(unlist(full_genes))
    all_fullg<-as(all_fullg,"list")
    seqinr::write.fasta(all_fullg,names = ass_nm,file.out=paste0(full_dir,"/","full_gene_",og,"_",genes[j],".fasta"))
    
  }
}


### Find truncated/pseudo genes ------------------------
library(Biostrings)
fls<-list.files("DefensNstress/AL/cds/sep_genes/v5",full.names=T,pattern=".fasta")
full_fls<-list.files("DefensNstress/AL/sequences/aug_ext/v5",full.names=T,pattern=".fasta")
out_dir<-"DefensNstress/AL/cds/sep_genes/v5/fixed/"
out_dfull<-"DefensNstress/AL/sequences/aug_ext/v5/fixed/"

ll<-lapply(fls,function(x){
  tm<-readDNAStringSet(x)
  lth<-width(tm)
  lmed<-median(lth)
  upl<-lmed*1.25
  lowl<-lmed*0.75
  pseuN<-which(lth < lowl | lth > upl)
  pseu<-names(tm)[pseuN]
  if(length(pseu)<1){pseu<-NA}
  og<-substr(basename(x),8,16)
  gene<-substr(basename(x),18,26)
  out<-cbind(og,gene,pseu)
  if(length(pseu)<2){if(!is.na(pseu)){tm<-tm[-pseuN]}}
  if(length(tm)<3){print(paste0(basename(x), " has less than 3 seqs"))}
  seqinr::write.fasta(as(tm,"list"),names=names(tm),file.out=paste0(out_dir,gsub(".fasta","_fixed.fasta",basename(x))))
  fulgfile<-full_fls[grep(gsub("cds_v5_","",basename(x)),full_fls)]
  ful<-readDNAStringSet(fulgfile)
  if(length(pseu)<2){if(!is.na(pseu))ful<-ful[-pseuN]}
  seqinr::write.fasta(as(ful,"list"),names=names(ful),file.out=paste0(out_dfull,gsub(".fasta","_fixed.fasta",basename(fulgfile))))
  return(out)
})
ll<-data.frame(do.call(rbind,ll))
write.csv(ll,"DefensNstress/AL/stat/AL_putative_pseudo_genes_based_on_cds_length_v5.csv",row.names=F)

### MAFFT ALIGNMENT OF AUGUSTUS PREDICTED FULL GENES ---------
# in fixed folders you find the fixed sequence by removing "pseudo" genes
# "pseudo" genes were determined by checking the cds lengths (see STATS)
input_dir="DefensNstress/AL/sequences/aug_ext/v5/fixed"
output_dir="DefensNstress/AL/align/full_genes/v5/sep_genes_mafft"

mkdir -p "$output_dir"
for file in "$input_dir"/*.fasta; do
base=$(basename "$file" .fasta)
mafft --auto --adjustdirection --anysymbol "$file" > "$output_dir/${base}_aligned_v5.fasta"
echo "Full gene alignment done for $base"
done

### MACSE ALIGNMENT OF GENE CDS ---------------
# in fixed folders you find the fixed sequence by removing "pseudo" genes
# "pseudo" genes were determined by checking the cds lengths (see STATS)
input_dir="DefensNstress/AL/cds/sep_genes/v5/fixed"
output_dir="DefensNstress/AL/align/cds/v5/sep_genes"
macse_jar="macse_v2.07.jar"

mkdir -p "$output_dir"
for file in "$input_dir"/*.fasta; do
base=$(basename "$file" .fasta)
java -jar "$macse_jar" -prog alignSequences \
-seq "$file" \
-out_NT "$output_dir/${base}_v5_aligned_NT.fasta" \
-out_AA "$output_dir/${base}_v5_aligned_AA.fasta" 
echo "Codon alignment done for $base"
done



## STATS ----------------------
### Pseudogene table -------------------------------
library(dplyr)
pse<-read.csv("data/DefenseNstress/AL/stats/AL_putative_pseudo_genes_based_on_cds_length_v5.csv")
summary_df <- pse %>%
  filter(!is.na(pseu)) %>%
  group_by(og) %>%
  summarise(
    genes = n_distinct(gene),
    freq_gene = gene %>% 
      table() %>%
      which.max() %>%
      names()
  )
all_ogs <- unique(pse$og)
summary_df <- tibble(og = all_ogs) %>%
  left_join(summary_df, by = "og") %>%
  mutate(
    genes = ifelse(is.na(genes), 0, genes),
    freq_gene = ifelse(is.na(freq_gene), NA, freq_gene)
  )
summary_df<-data.frame(summary_df)
# combine with cnv percentages 
cnv_perc<-read.csv("data/DefenseNstress/AL/cnv_within/cnv_percentage_per_family_v5.csv")
summary_df<-summary_df[match(cnv_perc$family,summary_df$og),]
summary_df_joined<-data.frame(cnv_perc,summary_df[,-1])
colnames(summary_df_joined)[6:7]<-c("pseudo_genes","frequent_pseudogene")
summary_df_joined$pse_perc<-summary_df_joined$pseudo_genes/summary_df_joined$refCount
write.csv(summary_df_joined,"data/DefenseNstress/AL/stats/AL_putative_pseudo_genes_counts_with_cnv_v5.csv",row.names = F)


### Pi / Theta --------------------------------------
#calculate pi/theta from mafft alignment of full genes
library(gSoup)
library(Biostrings)
algn_files<-list.files("data/DefenseNstress/AL/align/full_gene/v5",full.names=T,pattern=".fasta")

pb<-txtProgressBar(max=length(algn_files),style=3,width=50)
div_table<-NULL
for(i in seq_along(algn_files)){
  setTxtProgressBar(pb,i)
  algn_seq<-readDNAStringSet(algn_files[i])
  og<-substr(basename(algn_files[i]),11,19)
  gene<-substr(basename(algn_files[i]),21,29)
  if(length(algn_seq)<3){per_gene_dstat<-c(og,gene,rep(NA,7))
  div_table<-rbind(div_table,per_gene_dstat)
  next}
  dstats<-div.stats(algn_seq)
  per_gene_dstat<-c(og,gene,length(algn_seq),mean(width(algn_seq)),unlist(dstats[-5]))
  div_table<-rbind(div_table,per_gene_dstat)
}
div_table<-data.frame(div_table)
colnames(div_table)[1:4]<-c("gene_family","gene","no_seqs","mean_algn_len")
write.csv(div_table,"data/DefenseNstress/AL/stats/AL_pi_theta_table_v5.csv",row.names=F)



#### plot pi / theta --------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 

pi_theta <- read.csv("data/DefenseNstress/AL/stats/AL_pi_theta_table_v5.csv")
#  data
plot_data <- pi_theta %>%
  filter(!is.na(pi), !is.na(theta)) %>%
  pivot_longer(
    cols = c(pi, theta),
    names_to = "statistic", 
    values_to = "value"
  ) %>%
  mutate(statistic = factor(statistic, levels = c("pi", "theta")),
         statistic_label = ifelse(statistic == "pi", "π", "θ"))

# Main boxplot with points
main_plot <- ggplot(plot_data, aes(x = gene_family, y = value)) +
  geom_boxplot(
    aes(fill = statistic),
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = statistic),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.5,
    alpha = 0.6
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("π", "θ")) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("π", "θ")) +
  labs(x = "Gene Family", y = "Diversity Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

# Density plot
density_plot <- ggplot(plot_data, aes(x = value, fill = statistic)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(x = "Diversity Value", y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Combine plots
combined_plot <- main_plot + density_plot + 
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "π and θ Distributions by Gene Family",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  )

print(combined_plot)
# save pdf in 10 to 3.8 inch ratio




### Pi A/S -----------------
# remove "!" from alignments before reading into R << with bash
# MAC
# sed -i.bak '' 's/!/N/g' *.fasta
# # LINUX
# sed -i.bak 's/!/N/g' *.fasta

library(Biostrings)
library(gSoup)
algn_files<-list.files("data/DefenseNstress/AL/align/cds/v5",full.names=T,pattern=".fasta")
pb<-txtProgressBar(max=length(algn_files),style=3,width=50)
div_table<-NULL
for(i in seq_along(algn_files)){
  setTxtProgressBar(pb,i)
  algn_seq<-readDNAStringSet(algn_files[i])
  og<-substr(basename(algn_files[i]),9,16)
  gene<-substr(basename(algn_files[i]),18,26)
  if(length(algn_seq)<3){per_gene_dstat<-c(og,gene,rep(NA,13))
  div_table<-rbind(div_table,per_gene_dstat)
  next}
  dstats<-div.stats(algn_seq,pi_AS = T)
  per_gene_dstat<-c(og,gene,length(algn_seq),mean(width(algn_seq)),unlist(dstats[-5]))
  div_table<-rbind(div_table,per_gene_dstat)
}
div_table<-data.frame(div_table)
colnames(div_table)[1:4]<-c("gene_family","gene","no_seqs","mean_algn_len")
write.csv(div_table,"data/DefenseNstress/AL/stats/AL_pi_AS_table_v5.csv",row.names=F)


#### Pi A/S boxplot -----------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 
library(ggside) 

pias_tab<-read.csv("data/DefenseNstress/AL/stats/AL_pi_AS_table_v5.csv")

# Calculate log2(piA/piS) ratio
pias_ratio <- pias_tab %>%
  filter(!is.na(pi_AS.corrected_pi_syn), 
         !is.na(pi_AS.corrected_pi_nonsyn),
         pi_AS.corrected_pi_syn > 0) %>%
  mutate(log2_ratio = log2(pi_AS.corrected_pi_nonsyn / pi_AS.corrected_pi_syn))

# Create the plot with density
ggplot(pias_ratio, aes(x = gene_family, y = log2_ratio)) +
  # Main plot elements
  geom_boxplot(fill = "#4E79A7", alpha = 0.7, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "#4E79A7") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  
  # Right-side density plot (single combined density)
  geom_ysidedensity(
    aes(x = after_stat(density)),  
    fill = "#4E79A7", 
    alpha = 0.3,
    position = "identity"
  ) +
  
  #  scaling and theming
  scale_ysidex_continuous(expand = expansion(mult = c(0, 0.1)), breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  
  # Labels and theme
  labs(
    x = "Gene Family",
    y = "log2(πN/πS)",
    title = "Nonsynonymous to Synonymous Diversity Ratio",
    subtitle = "Positive values indicate higher nonsynonymous diversity (πN > πS)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    ggside.panel.scale = 0.2  
  )

# save pdf in 10 x 3.8 inch ratio

