######### CNV WITH ALL 231 GENE FAMILIES ############
library(Biostrings)

setwd("/Users/piyalkaru/Desktop/DDORF/Ann/231ogs/cnv_geneFam")

# THALIANA -------------------------
# gene ids of gene families
ogList<-data.table::fread("data/ogList_with_functions_v2.csv")
gene_ids<-data.table::fread("data/Orthogroups.tsv",h=T)
gene_ids231<-gene_ids[gene_ids$Orthogroup%in%ogList$Family.ID,]

#drop non-chromosomal gene families
to_remove<-grep("ATM",gene_ids231$Athaliana_Araport11.priTranscripts)
gene_ids231<-gene_ids231[-to_remove,]

## Extract gene sequences and save in FASTA  ------------------------------
library(Biostrings)

out_dir="data/231/AT/refseq/"
dir.create(out_dir,recursive = T)

at_ids<-gene_ids231[,c(1,5)]
# genome and the gff files
at_genome<-readDNAStringSet("data/TAIR10_allChr.fas")
athal_gff <-data.frame(biomartr::read_gff("data/Araport11_GFF3_genes_transposons.Feb2022.gff"))

pb<-txtProgressBar(max=nrow(at_ids),width = 50,style = 3)
for(j in 1:nrow(at_ids)){
  setTxtProgressBar(pb,j)
  newIds<-stringr::str_split_fixed(unlist(stringr::str_split(at_ids[j,2],pattern="\\,")),pattern="\\.",n=2)[,1]
  gene_list<-list()
  for(i in seq_along(newIds)){
    g_match<-athal_gff[grep(gsub("\\ ","",newIds[i]),athal_gff$attribute),] #AT1G70020
    rr<-which(g_match$type=="gene")[1]
    coords<-range(g_match[rr,4:5])
    chr<-g_match[rr,1]
    gene<-at_genome[[grep(chr,names(at_genome))]][min(coords):max(coords)]
    gene_list[[i]]<-gene
  }
  seqinr::write.fasta(gene_list,names = newIds,file.out=paste0(out_dir,unlist(at_ids[j,1]),".fasta"))
}


### find gene lengths --------
fls<-list.files("data/231/AT/refseq/",full.names=T)
at_gen_length<-do.call(rbind,lapply(fls,function(x){dd<-readDNAStringSet(x);return(cbind(names(dd),Biostrings::width(dd)))}))
at_gen_length<-data.frame(id=at_gen_length[,1],length=at_gen_length[,2])
at_gen_length$length<-as.numeric(at_gen_length$length)
at_gen_length$id<-as.character(stringr::str_trim(at_gen_length$id))
write.table(at_gen_length,"data/231/AT/AT_231_gene_lengths.txt",row.names=F, quote=F,sep="\t")

## BLAST FOR GENES --------------------------------
reference_dir="data/231/AT/refseq"
seq_dir="data/136_LR_fnas"
out_dir="data/231/AT/blast1/"

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
b_files<-list.files("data/231/AT/blast1/",full.names=T)
og_ids<-at_ids$Orthogroup
out_dir="data/231/AT/blast1/og_sorted/"
dir.create(out_dir,recursive = T)
col_names<-c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","stitle")

pb<-txtProgressBar(max=length(og_ids),style = 3, width=50)
for(j in seq_along(og_ids)){
  setTxtProgressBar(pb,j)
  og_files<-b_files[grep(og_ids[j],b_files)]
  og_added<-NULL
  for(i in seq_along(og_files)){
    b_hits<-read.table(og_files[i],comment.char = "#",col.names = col_names)
    if(nrow(b_hits)<1) next
    b_hits$stitle<-stringr::str_split_fixed(gsub(".blastn","",basename(og_files[i])),pattern="_",n=2)[,2]
    og_added<-rbind(og_added,b_hits)
  }
  write.table(og_added,paste0(out_dir,og_ids[j],"_bhits_v1.txt"),quote=F,row.names=F,sep="\t")
}

## GROUP THE BLAST HITS  -----------------------------
out_dir<-"data/231/AT/blast1/length_sorted/"
dir.create(out_dir,recursive = T)

at_gen_length<-read.table("data/231/AT_231_gene_lengths.txt",h=T)
og_sorted<-list.files("data/231/AT/blast1/og_sorted/",full.names = T,pattern = ".txt")
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
        tmp3<-tmp3[tmp3$length>(gen_len*.5),] # keep hits bigger than 50% of gene length
        if (nrow(tmp3) > 0) {
          tmp3$copy <- seq_len(nrow(tmp3)) + copy_counter
          copy_counter <- copy_counter + nrow(tmp3)  # Update the copy counter
          out <- rbind(out, tmp3)
        }
      }
    }
  }
  write.table(out,paste0(out_dir,og,"_length_sorted_v1.txt"),quote=F,sep="\t",row.names=F)
}


## CNV within species ------------
# remove pseudo-duplicates >> where the same region is flagged for two different genes
out_dir<-"data/231/AT/blast1/dup_remd/"
dir.create(out_dir,recursive = T)

bp_fl<-list.files("data/231/AT/blast1/length_sorted/",full.names = T,pattern = "length_sorted_v1.txt")
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
  write.table(tm0,paste0(out_dir,gsub("length_sorted","pseu_dup_rem",basename(bp_fl[k]))),row.names = F,quote=F,sep="\t")
}

#### list gene families with CNV ------------------
out_dir<-"data/231/AT/cnv_within/"
dir.create(out_dir,recursive = T)

du_rm_fl<-list.files("data/231/AT/blast1/dup_remd/",full.names=T,pattern="pseu_dup_rem_v1.txt")
ogList<-read.csv("data/ogList_with_functions_v2.csv")

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
write.csv(out_cnv,file.path(out_dir,"Athal231_sig_fam_gain_loss_v1.csv"),row.names = F)


#### Assembly vs gene count per gene family -----------------
library(tidyr)
library(dplyr)
out_dir<-"data/231/AT/cnv_within/"
dir.create(out_dir,recursive = T)

out_cnv<-read.csv("data/231/AT/cnv_within/Athal231_sig_fam_gain_loss_v1.csv")
matrix_out <- out_cnv %>%
  select(family_id, assembly, total_gene_count) %>%
  pivot_wider(names_from = assembly, values_from = total_gene_count)
og_vs_ass<-matrix_out
og_vs_ass[is.na(og_vs_ass)]<-0
write.table(og_vs_ass,"data/231/AT/cnv_within/Athaliana231_geneFamily_vs_assembly_cnv_v1.txt",quote=F,row.names = F, sep="\t")


#### histograms of CNV per family -----------
# **** improved plotting is in script plotting.R ****

### TE ASSESSMENT -------------------------------------
# **** full workflow is in TE_assessment.R script ****

### CNV ASSESSMENT ---------------------------------
# CNV assessment of both within species and between species 
# is in the script cnv_assmnt2.R

### ECOLOGY ----------------------------------------
# CNV-ENV assessment 
# the analysis is in the script cnv-env_test.R



### SEQUENCE ANALYSIS ------------------------------
#### Extract genes from assemblies ---------
library(biomartr)
library(seqinr)
library(Biostrings)

out_dir<-"data/231/AT/blast1/sequences"
if(!dir.exists(out_dir)){dir.create(out_dir,showWarnings = F, recursive = T)}

genomes_path<-list.files("data/AT/136_LR_fnas",full.names = T)
fam_files<-list.files("data/231/AT/blast1/dup_remd/",full.names=T,pattern="dup_rem")

pb<-txtProgressBar(max=length(fam_files),style = 3, width=50)
for(f in 1:length(fam_files)){ 
  setTxtProgressBar(pb,f)
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
      genome_file<-genomes_path[grep(assembly[k],genomes_path)]
      genome<-read_genome(genome_file)
      
      scaffs<-unique(tmp1$scaffold)
      gene_list<-list()
      for(i in seq_along(scaffs)){
        tmp2<-tmp1[tmp1$scaffold==scaffs[i],]
        chr<-genome[[grep(scaffs[i],names(genome))]]
        dnastrin<-DNAStringSet()
        for(m in 1:nrow(tmp2)){
          flank<-tmp2[m,"length"]*.2
          dnastrin[[m]]<-chr[(tmp2$start[m]-flank):(tmp2$end[m]+flank)]
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
  write.fasta(as(ggs,"list"),names=nm,file.out = file.path(out_dir,paste0(og,"_231_all_genes_v1.fasta")))
}

#### PSEUDO & FUNCTIONAL ASSESSMENT ------------------

# step 1 ------------
# builds a master per-copy locus table >> 01_build_copy_locus_db.py
cd ./scripts

tt_tables="data/231/AT/blast1/dup_remd/"
reference_og="data/231/AT/refseq/"
extracted_copy_fastas="data/231/AT/blast1/sequences/"
output_step1="data/231/AT/blast1/pseu_func/step1/"

mkdir -p $output_step1

python 01_build_copy_locus_db.py \
--tt-dir $tt_tables \
--ref-dir $reference_og \
--copy-dir $extracted_copy_fastas \
--out-dir $output_step1

# step 2 ------------
# does a reference-guided alignment of each extracted locus to its matching reference gene sequence, chooses the best strand, and outputs alignment metrics you can use for the next classification step >> 02_align_copies_to_reference.py

# This script:
# •	reads the master table from script 1
# •	loads the copy fasta and reference fasta
# •	finds the matching reference gene sequence for each copy
# •	aligns the extracted copy to the reference on both strands
# •	keeps the better strand
# •	outputs alignment metrics you will need for the next step
cd ./scripts
copy_locus_master="data/231/AT/blast1/pseu_func/step1/copy_locus_master.tsv"
out_dir="data/231/AT/blast1/pseu_func/step2/"
align_out="data/231/AT/blast1/pseu_func/step2/algnments"

python 02_align_copies_to_reference_checkpointed.py \
--master $copy_locus_master \
--out-dir $align_out \
--write-oriented-fasta \
--write-alignment-strings \
--batch-size 200 \
--progress-every 100

# step 3 ---------------
# Compute a structural-integrity profile for every copy
# 
# For each candidate copy, record at least:
#   •	CDS recovered: yes/no
# •	ORF recovered: yes/no
# •	translated protein length
# •	protein length / reference protein length
# •	start codon intact: yes/no
# •	stop codon intact: yes/no
# •	number of internal stop codons
# •	number of frameshifts
# •	exon count compared with reference
# •	alignment coverage of CDS
# •	alignment identity of protein
# •	core domain presence/absence
# •	key motif or catalytic residue retention, if known for that family
# 
# Script 3: re-annotate each extracted locus against the matching reference gene using a spliced protein-to-genome aligner, recover a predicted CDS/protein, and assign an initial structural class
# P = putatively pseudogenized, U0 = intact candidate, A = ambiguous

# python 03_reannotate_and_structural_classify.py \
# --master /path/to/output_step1/copy_locus_master.tsv \
# --step2 /path/to/output_step2/copy_reference_alignment_metrics.tsv \
# --out-dir /path/to/output_step3

#>>>> extract CDS for this step to work <<<<<<

cd ./scripts
gene_fam_ref_dir="data/231/AT/refseq/full_genes"
gff="data/Araport11_GFF3_genes_transposons.Feb2022.gff"
genome="data/TAIR10_allChr.fas"
out_dir="data/231/AT/refseq_cds"

mkdir -p $out_dir

# CDS extraction
# this approach produces the least number of stop codons
cd ./scripts

mkdir -p data/231/AT/refseq_cds_helixer

python 02b_extract_ref_CDS_helixer_v3.py \
--in-dir data/231/AT/refseq/full_genes \
--out-dir data/231/AT/refseq_cds_helixer \
--species Arabidopsis_thaliana

# run step 3
cd ./scripts
copy_locus_master="data/231/AT/blast1/pseu_func/step1/copy_locus_master.tsv"
alignment_metrics="data/231/AT/blast1/pseu_func/step2/algnments/copy_reference_alignment_metrics.tsv"
out_dir="data/231/AT/blast1/pseu_func/step3/"

mkdir -p $out_dir

python 03_reannotate_and_structural_classify_v2.py \
--master $copy_locus_master \
--step2 $alignment_metrics \
--out-dir $out_dir


# step 4 ----------
# among the structurally intact copies, classify
# F = putatively functionally changed vs U = apparently unchanged
# using sequence divergence relative to the matched reference gene across assemblies
# 
# This script uses only the copies that came out as U0 from step 3.
# 
# It is deliberately conservative. It does not claim true neo- or subfunctionalization. It just marks an intact copy as F if it is unusually diverged from the matched reference gene compared with other intact copies of that same gene across assemblies.
# python 04_classify_functional_change.py \
# --step3 /path/to/output_step3/step3_structural_classification.tsv \
# --out-dir /path/to/output_step4

cd ./scripts
structural_classification="data/231/AT/blast1/pseu_func/step3/step3_structural_classification.tsv"
out_dir="data/231/AT/blast1/pseu_func/step4/"

mkdir -p $out_dir

python 04_classify_functional_changes.py \
--step3 $structural_classification \
--out-dir $out_dir


#With explicit thresholds:
python 04_classify_functional_change.py \
--step3 /path/to/output_step3/step3_structural_classification.tsv \
--out-dir /path/to/output_step4 \
--abs-pid-cutoff 90 \
--len-ratio-low 0.70 \
--len-ratio-high 1.30 \
--rz-cutoff 3


# step 5 ------------
# a script that takes the final classification table and writes the filtered FASTA files you need for downstream population-genetic analyses, plus summary tables that quantify pseudogenization and functional-change prevalence per gene and per family.
# 
# How to use the outputs
# 
# For diversity statistics:
#   •	use popgen_input_U_only.tsv or fastas/U_only.fasta for the strictest dataset
# •	use popgen_input_U_plus_F.tsv or fastas/U_plus_F.fasta for all putatively functional copies
# •	use excluding_P if you want to keep ambiguous copies but drop pseudogenes
# 
# For the CNV vs pseudogenization analysis:
#   •	family_compact_summary_for_stats.tsv is the cleanest starting table
# 
# Useful columns there:
#   •	n_total_copies
# •	n_P
# •	n_F
# •	n_U
# •	prop_pseudogenized
# •	prop_functional_nonpseudo
# •	prop_changed_among_nonpseudo
# 
# Save this as 05_filter_and_export_for_popgen.py.
# python 05_filter_and_export_for_popgen.py \
# --step4 /path/to/output_step4/step4_final_copy_classification.tsv \
# --copy-fasta-dir /path/to/extracted_copy_fastas \
# --out-dir /path/to/output_step5
cd ./scripts

step4_final_copy_classification="data/231/AT/blast1/pseu_func/step4/step4_final_copy_classification.tsv"
extracted_copy_fastas="data/231/AT/blast1/sequences"
out_dir="data/231/AT/blast1/pseu_func/step5/"

mkdir -p $out_dir

python 05_filter_and_export_for_popgen.py \
--step4 $step4_final_copy_classification \
--copy-fasta-dir $extracted_copy_fastas \
--out-dir $out_dir
# ***************************************************************
# *** visualization of these results are in plotting.R script *** 
# ***************************************************************

# step 6 -----------------
# It does three things:
#   1.	tests whether families with higher copy number show more pseudogenization or more putative functional change
# 2.	summarizes per-family copy-fate patterns in analysis-ready tables
# 3.	compares population-genetic statistics across datasets if you already have per-family π / θ / Tajima’s D results for:
# •	U_only
# •	U_plus_F
# •	optionally all or excluding_P
# 
# It uses simple, interpretable models:
# •	Spearman correlations
# •	binomial GLMs on counts
# •	optional paired comparisons for popgen outputs
# python 06_analyze_cnv_vs_copy_fate_and_popgen.py \
# --family-summary /path/to/output_step5/family_compact_summary_for_stats.tsv \
# --popgen-u /path/to/popgen_U_only.tsv \
# --popgen-uf /path/to/popgen_U_plus_F.tsv \
# --popgen-excl-p /path/to/popgen_excluding_P.tsv \
# --popgen-all /path/to/popgen_all.tsv \
# --out-dir /path/to/output_step6

cd ./scripts

compact_summary="data/231/AT/blast1/pseu_func/step5/family_compact_summary_for_stats.tsv"
popgen_U_only="data/231/AT/blast1/pseu_func/step5/popgen_input_U_only.tsv"
popgen_U_plus_F="data/231/AT/blast1/pseu_func/step5/popgen_input_U_plus_F.tsv"
popgen_excluding_P="data/231/AT/blast1/pseu_func/step5/popgen_input_excluding_P.tsv"
popgen_input_all="data/231/AT/blast1/pseu_func/step5/popgen_input_all.tsv"
out_dir="data/231/AT/blast1/pseu_func/step6/"

mkdir -p $out_dir

python 06_analyze_cnv_vs_copy_fate_and_popgen.py \
--family-summary $compact_summary \
--popgen-u $popgen_U_only \
--popgen-uf $popgen_U_plus_F \
--popgen-excl-p $popgen_excluding_P \
--popgen-all $popgen_input_all \
--out-dir $out_dir





# LYRATA ---------------------------

library(Biostrings)

setwd("/Users/piyalkaru/Desktop/DDORF/Ann/231ogs/cnv_geneFam")

# gene ids of gene families
ogList<-data.table::fread("data/ogList_with_functions_v2.csv",h=T)
gene_ids<-data.table::fread("data/Orthogroups.tsv",h=T)
gene_ids231<-gene_ids[gene_ids$Orthogroup%in%ogList$Family.ID,]

#drop non-chromosomal gene families
to_remove<-grep("ATM",gene_ids231$Athaliana_Araport11.priTranscripts)
gene_ids231<-gene_ids231[-to_remove,]
write.csv(gene_ids231,"data/all_0.01_ogs_w_o_ATM_genes.csv",row.names = F)

## Extract gene sequences and save in FASTA  ------------------------------

out_dir="data/231/AL/refseq/"
dir.create(out_dir,recursive = T)

## Extract gene sequences and save in FASTA 
al_ids<-gene_ids231[,c(1,4)]
og_ids<-gene_ids231$Orthogroup
# genome and the gff files
al_gff <- biomartr::read_gff("data/Alyrata_384_v2.1.gene.gff3")
al_genome <- readDNAStringSet("data/Alyrata_384_v1.fa")

pb<-txtProgressBar(max=nrow(al_ids),style = 3, width=50)
for(j in 1:nrow(al_ids)){
  setTxtProgressBar(pb,j)
  newIds<-stringr::str_split_fixed(unlist(stringr::str_split(al_ids[j,2],pattern="\\,")),pattern="\\.",n=2)[,1]
  if(length(newIds)<1 & newIds[1]=="") next
  gene_list<-list()
  for(i in seq_along(newIds)){
    g_match<-al_gff[grep(gsub("\\ ","",newIds[i]),al_gff$attribute),]
    rr<-which(g_match$type=="gene")[1]
    coords<-range(g_match[rr,4:5])
    chr<-as.character(g_match[rr,1])
    gene<-al_genome[[which(names(al_genome)==chr)]][min(coords):max(coords)]
    gene_list[[i]]<-gene
  }
  seqinr::write.fasta(gene_list,names = newIds,file.out=paste0(out_dir,unlist(al_ids[j,1]),".fasta"))
}


### find gene lengths --------
fls<-list.files("data/231/AL/refseq/",full.names=T)
al_gen_length<-do.call(rbind,lapply(fls,function(x){dd<-readDNAStringSet(x);return(cbind(names(dd),width(dd),gsub(".fasta","",basename(x))))}))
al_gen_length<-data.frame(id=al_gen_length[,1],length=al_gen_length[,2],og=al_gen_length[,3])
al_gen_length$length<-as.numeric(al_gen_length$length)
al_gen_length$id<-as.character(stringr::str_trim(al_gen_length$id))
write.table(al_gen_length,"data/231/AL/all_gene_lengths.txt",row.names = F,quote=F,sep="\t")


## BLAST FOR GENES --------------------------------
reference_dir="data/231/AL/refseq"
seq_dir="data/Alyrata_assmbl"
out_dir="data/231/AL/blast1/"

mkdir -p $out_dir
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
out_dir<-"data/231/AL/blast1/og_sorted/"
dir.create(out_dir,recursive = T)
col_names <-c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","stitle")

b_files<-list.files("data/231/AL/blast1",full.names = T)
og_ids<-al_ids$Orthogroup
for(j in seq_along(og_ids)){
  og_files<-b_files[grep(og_ids[j],b_files)]
  og_added<-NULL
  for(i in seq_along(og_files)){
    b_hits<-read.table(og_files[i],comment.char = "#",col.names = col_names)
    if(nrow(b_hits)<1) next
    b_hits$stitle<-gsub(paste0(og_ids[j],"|.blastn"),"",basename(og_files[i]))
    og_added<-rbind(og_added,b_hits)
  }
  write.table(og_added,paste0(out_dir,og_ids[j],"_bhits_v1.txt"),quote=F,row.names=F,sep="\t")
}
gSoup::blast_plot(b_hits,cols=c("stitle","sseqid","sstart","send","qseqid"))


## GROUP THE BLAST HITS  -----------------------------
out_dir<-"data/231/AL/blast1/length_sorted/"
dir.create(out_dir,recursive = T)

al_gen_length<-data.frame(data.table::fread("data/231/AL/all_gene_lengths.txt",h=T))
og_sorted<-list.files("data/231/AL/blast1/og_sorted/",full.names = T)
pb<-txtProgressBar(min=0,max=length(og_sorted),style = 3,width = 50)
for(m in seq_along(og_sorted)){
  setTxtProgressBar(pb,m)
  all_blast<-read.table(og_sorted[m],header = T)
  all_blast<-all_blast[all_blast$pident>=90,]
  og<-gsub("_bhits_v1.txt","",basename(og_sorted[m]))
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
        tmp3<-data.frame(assembly=assembly[i],gene=gene[j],scaffold=chr[k],start=rng[,1],end=rng[,2],length=rng[,2]-rng[,1],max_length=unlist(max_len))
        tmp3<-tmp3[tmp3$length>(gen_len*.5),]
        
        if (nrow(tmp3) > 0) {
          tmp3$copy <- seq_len(nrow(tmp3)) + copy_counter
          copy_counter <- copy_counter + nrow(tmp3)  
          out <- rbind(out, tmp3)
        }
      }
    }
  }
  write.table(out,paste0(out_dir,og,"_length_sorted_v1.txt"),quote=F,sep="\t",row.names=F)
}


## CNV within species ------------
# remove pseudo-duplicates >> where the same region is flagged for two different genes
out_dir<-"data/231/AL/blast1/dup_remd/"
dir.create(out_dir,recursive = T)

bp_fl<-list.files("data/231/AL/blast1/length_sorted/",full.names = T,pattern = "_v1.txt")
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
  write.table(tm0,paste0(out_dir,gsub("length_sorted","pseu_dup_rem",basename(bp_fl[k]))),row.names = F,quote=F,sep="\t")
}


#### list gene families with CNV ------------------
out_dir<-"data/231/AL/cnv_within/"
dir.create(out_dir,recursive = T)

du_rm_fl<-list.files("data/231/AL/blast1/dup_remd/",full.names=T, pattern="_v1.txt")
ogList<-read.csv("data/ogList_with_functions.csv")
gene_ids231<-read.csv("data/all_0.01_ogs_w_o_ATM_genes.csv",h=T)
out_cnv<-NULL
cnv_count<-NULL
for(i in seq_along(du_rm_fl)){
  og_cnv<-read.table(du_rm_fl[i],h=T)
  og_id<-substr(basename(du_rm_fl[i]),1,9)
  og_line<-gene_ids231[gene_ids231$Orthogroup==og_id,]
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
write.csv(out_cnv,paste0(out_dir,"Alyrata231_sig_fam_gain_loss_v1.csv"),row.names = F)

#### Assembly vs gene count per gene family -----------
library(tidyr)
library(dplyr)

out_cnv<-read.csv("data/231/AL/cnv_within/Alyrata231_sig_fam_gain_loss_v1.csv")
matrix_out <- out_cnv %>%
  select(family_id, assembly, total_gene_count) %>%
  pivot_wider(names_from = assembly, values_from = total_gene_count)
og_vs_ass<-matrix_out
og_vs_ass[is.na(og_vs_ass)]<-0
write.table(og_vs_ass,"data/231/AL/cnv_within/Alyrata231_geneFamily_vs_assembly_cnv_v1.txt",quote=F,row.names = F, sep="\t")


#### CNV per family table ----------------------
library(dplyr)
out_dir<-"data/231/AL/cnv_within/"
dir.create(out_dir,recursive = T)

out_cnv<-read.csv("data/231/AL/cnv_within/Alyrata231_sig_fam_gain_loss_v1.csv")
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
write.csv(ogs_cnv,paste0(out_dir,"cnv_percentage_per_family_v1.csv"),row.names=F)

#### histograms of CNV per family -----------
# **** improved plotting is in script plotting.R ****


### TE ASSESSMENT -------------------------------------
# **** full workflow is in TE_assessment.R script ****


### CNV ASSESSMENT ---------------------------

# CNV assessment of both within species and between species 
# is in the script cnv_assmnt2.R

### ECOLOGY ----------------------------------------
# CNV-ENV assessment 
# the analysis is in the script cnv-env_test.R


#prepare AL assembly names
al_ogass<-read.table("data/231/AL/cnv_within/Alyrata231_geneFamily_vs_assembly_cnv_v1.txt",h=T)
cnm<-c("family_id","73_3a","AL08","al1","AL27","ALR","BAM12","BOR","F1_14","KN003","KN005","KN006","MN47","NLA","NT1","NT12","NT8","NT9","PU6","TE11","TE4","TE8","WS1")
al_ogass<-al_ogass[,-14] # drop NT1 << reference bias
colnames(al_ogass)<-cnm
write.table(al_ogass,"data/231/AL/cnv_within/Alyrata231_geneFamily_vs_assembly_cnv_v1.1.txt",row.names = F, quote = F, sep = "\t")



### SEQUENCE ANALYSIS ------------------------------
#### Extract blasted genes from assemblies ---------
library(biomartr)
library(seqinr)
library(Biostrings)
genomes_path<-list.files("data/Alyrata_assmbl",full.names = T)
fam_files<-list.files("data/231/AL/blast1/dup_remd",full.names=T,pattern = "dup_rem_v1.txt")

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
          flank<-tmp2[m,"length"]*.2
          dnastrin[[m]]<-chr[(tmp2$start[m]-flank):(tmp2$end[m]+flank)]
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
  write.fasta(as(ggs,"list"),names=names(ggs),file.out = paste0("data/231/AL/blast1/sequences/",og,"_all_genes_v1.fasta"))
}

#### PSEUDO & FUNCTIONAL ASSESSMENT ------------------
# step 1 ------
#builds a master per-copy locus table >> 01_build_copy_locus_db.py
cd data/runs

tt_tables="data/231/AL/blast1/dup_remd/"
reference_og="data/231/AL/refseq/"
extracted_copy_fastas="data/231/AL/blast1/sequences/"
output_step1="data/231/AL/blast1/pseu_func/step1/"

python 01_build_copy_locus_db_alyrata.py \
--tt-dir $tt_tables \
--ref-dir $reference_og \
--copy-dir $extracted_copy_fastas \
--out-dir $output_step1

# Step 2 ----
# does a reference-guided alignment of each extracted locus to its matching reference gene sequence, chooses the best strand, and outputs alignment metrics you can use for the next classification step >> 02_align_copies_to_reference.py

# This script:
# •	reads the master table from script 1
# •	loads the copy fasta and reference fasta
# •	finds the matching reference gene sequence for each copy
# •	aligns the extracted copy to the reference on both strands
# •	keeps the better strand
# •	outputs alignment metrics you will need for the next step
copy_locus_master="data/231/AL/blast1/pseu_func/step1/copy_locus_master.tsv"
out_dir="data/231/AL/blast1/pseu_func/step2/"

python 02_align_copies_to_reference_lyrata.py \
--master $copy_locus_master \
--out-dir $out_dir \
--write-oriented-fasta 


# Step 3 -----
# re-annotate each extracted locus against the matching reference gene using a spliced protein-to-genome aligner, recover a predicted CDS/protein, and assign an initial structural class
# P = putatively pseudogenized, U0 = intact candidate, A = ambiguous

# python 03_reannotate_and_structural_classify_alyrata_compatible.py \
# --master /path/to/copy_locus_master.tsv \
# --step2 /path/to/copy_reference_alignment_metrics.tsv \
# --out-dir /path/to/step3

#>>>> extract CDS for this step to work <<<<<<
cd data/runs
gene_fam_ref_dir="data/231/AL/refseq"
reference_genome="/home/pkaru/hpc-project/CNV/refs/Alyrata_384_v2.1.gene.gff3"
out_dir="data/231/AL/refseq_cds"

mkdir -p $out_dir

# this approach produces the least number of stop codons
python 02b_extract_ref_CDS_helixer_v3.py \
--in-dir data/231/AL/refseq/full_genes \
--out-dir data/231/AL/refseq_cds_helixer \
--species Arabidopsis_lyrata


# run step3
cd data/runs
copy_locus_master="data/231/AL/blast1/pseu_func/step1/copy_locus_master.tsv"
alignment_metrics="data/231/AL/blast1/pseu_func/step2/alignments/copy_reference_alignment_metrics.tsv"
out_dir="data/231/AL/blast1/pseu_func/step3/"

mkdir -p $out_dir

python 03_reannotate_and_structural_classify_alyrata_compatible_v2.py \
--master $copy_locus_master \
--step2 $alignment_metrics \
--out-dir $out_dir


# Script 4: among the structurally intact copies, classify ----
# F = putatively functionally changed vs U = apparently unchanged
# using sequence divergence relative to the matched reference gene across assemblies
# 
# This script uses only the copies that came out as U0 from step 3.
# 
# It is deliberately conservative. It does not claim true neo- or subfunctionalization. It just marks an intact copy as F if it is unusually diverged from the matched reference gene compared with other intact copies of that same gene across assemblies.
# python 04_classify_functional_change.py \
# --step3 /path/to/output_step3/step3_structural_classification.tsv \
# --out-dir /path/to/output_step4

cd data/runs
structural_classification="data/231/AL/blast1/pseu_func/step3/step3_structural_classification.tsv"
out_dir="data/231/AL/blast1/pseu_func/step4/"

mkdir -p $out_dir

python 04_classify_functional_changes.py \
--step3 $structural_classification \
--out-dir $out_dir


# step 5 (5-10) ----
# a script that takes the final classification table and writes the filtered FASTA files you need for downstream population-genetic analyses, plus summary tables that quantify pseudogenization and functional-change prevalence per gene and per family.
# 
# How to use the outputs
# 
# For diversity statistics:
#   •	use popgen_input_U_only.tsv or fastas/U_only.fasta for the strictest dataset
# •	use popgen_input_U_plus_F.tsv or fastas/U_plus_F.fasta for all putatively functional copies
# •	use excluding_P if you want to keep ambiguous copies but drop pseudogenes
# 
# For the CNV vs pseudogenization analysis:
#   •	family_compact_summary_for_stats.tsv is the cleanest starting table
# 
# Useful columns there:
#   •	n_total_copies
# •	n_P
# •	n_F
# •	n_U
# •	prop_pseudogenized
# •	prop_functional_nonpseudo
# •	prop_changed_among_nonpseudo
# 
# Save this as 05_filter_and_export_for_popgen.py.
# python 05_filter_and_export_for_popgen.py \
# --step4 /path/to/output_step4/step4_final_copy_classification.tsv \
# --copy-fasta-dir /path/to/extracted_copy_fastas \
# --out-dir /path/to/output_step5
cd data/runs

step4_final_copy_classification="data/231/AL/blast1/pseu_func/step4/step4_final_copy_classification.tsv"
extracted_copy_fastas="data/231/AL/blast1/sequences"
out_dir="data/231/AL/blast1/pseu_func/step5/"

mkdir -p $out_dir

python 05_filter_and_export_for_popgen.py \
--step4 $step4_final_copy_classification \
--copy-fasta-dir $extracted_copy_fastas \
--out-dir $out_dir


# step 6 -----------------
# It does three things:
#   1.	tests whether families with higher copy number show more pseudogenization or more putative functional change
# 2.	summarizes per-family copy-fate patterns in analysis-ready tables
# 3.	compares population-genetic statistics across datasets if you already have per-family π / θ / Tajima’s D results for:
# •	U_only
# •	U_plus_F
# •	optionally all or excluding_P
# 
# It uses simple, interpretable models:
# •	Spearman correlations
# •	binomial GLMs on counts
# •	optional paired comparisons for popgen outputs
# python 06_analyze_cnv_vs_copy_fate_and_popgen.py \
# --family-summary /path/to/output_step5/family_compact_summary_for_stats.tsv \
# --popgen-u /path/to/popgen_U_only.tsv \
# --popgen-uf /path/to/popgen_U_plus_F.tsv \
# --popgen-excl-p /path/to/popgen_excluding_P.tsv \
# --popgen-all /path/to/popgen_all.tsv \
# --out-dir /path/to/output_step6

cd data/runs

compact_summary="data/231/AL/blast1/pseu_func/step5/family_compact_summary_for_stats.tsv"
popgen_U_only="data/231/AL/blast1/pseu_func/step5/popgen_input_U_only.tsv"
popgen_U_plus_F="data/231/AL/blast1/pseu_func/step5/popgen_input_U_plus_F.tsv"
popgen_excluding_P="data/231/AL/blast1/pseu_func/step5/popgen_input_excluding_P.tsv"
popgen_input_all="data/231/AL/blast1/pseu_func/step5/popgen_input_all.tsv"
out_dir="data/231/AL/blast1/pseu_func/step6/"

mkdir -p $out_dir

python 06_analyze_cnv_vs_copy_fate_and_popgen.py \
--family-summary $compact_summary \
--popgen-u $popgen_U_only \
--popgen-uf $popgen_U_plus_F \
--popgen-excl-p $popgen_excluding_P \
--popgen-all $popgen_input_all \
--out-dir $out_dir

