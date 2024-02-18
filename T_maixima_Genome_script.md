# T. maxima Genome Script

Li, R., Leiva, C., Lemer, S., Kirkendale, L., & Li, J.  (n.d.). Photosymbiosis profoundly influenced animal genome architecture and gene evolution – Insights from the giant clam genome. Manuscript submitted for publication.


## 1. Gene Family Analyses

OrthoDB Pipeline

```bash
docker pull ezlabgva/orthologer:v3.0.0
docker run -u $(id -u) -v $(pwd)/odb:/odbwork ezlabgva/orthologer:v3.0.0 orthomapper -c run -p Tmax -f Tmax_ReAnn.fasta -n 6447


# CAFE
# Run with Universal lambda
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 1 -p -c 6 -o k1p_1lambda&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 2 -p -c 6 -o k2p_1lambda&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 3 -p -c 6 -o k3p_1lambda&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 4 -p -c 6 -o k4p_1lambda&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 5 -p -c 6 -o k5p_1lambda&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -k 6 -p -c 6 -o k6p_1lambda&

# Run with 2 seperate lambda for Bivalvia and Gastropoda
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 1 -p -c 6 -o k1p_2lambdas&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 2 -p -c 6 -o k2p_2lambdas&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 3 -p -c 6 -o k3p_2lambdas&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 4 -p -c 6 -o k4p_2lambdas&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 5 -p -c 6 -o k5p_2lambdas&
nohup cafe5 -i BivGas_GeneCount.tsv -t BivGas_ultrametric.tre -y Biv_Gas_separate_lambda.txt -k 6 -p -c 6 -o k6p_2lambdas&

# Results Mining
# Only viewing the trees with significant changes
echo $'#nexus\nbegin trees;'>Significant_trees.tre
grep "*" Gamma_asr.tre >>Significant_trees.tre
echo "end;">>Significant_trees.tre # 6794 trees

## Only viewing the trees significant changes in Tmax
echo $'#nexus\nbegin trees;'>Tmax_Significant_trees.tre
grep -F "<17>*" Gamma_asr.tre >>Tmax_Significant_trees.tre
echo "end;">>Tmax_Significant_trees.tre # 207 trees
# Get the gene family IDs that have significant changes in Tmax
grep -o "\w*at6447\w*" Tmax_Significant_trees.tre > Tmax_Significant_trees_IDs #207

# Get the gene family expaned in Tmax
cat Gamma_change.tab |cut -f1,19|grep "+[1-9]" >Tmax.expanded # 2783
# Get the gene family extracted in Tmax
cat Gamma_change.tab |cut -f1,19|grep "-" >Tmax.contracted # 3490

############ Significant Branch on the Tree ########################
#  Get the gene family SIGNIFICANTLY expaned in Tmax
grep -Fwf Tmax_Significant_trees_IDs Tmax.expanded | cut -f1 > Tmax.expanded.Sig.treebranch #122
#  Get the gene family SIGNIFICANTLY contracted in Tmax
grep -Fwf Tmax_Significant_trees_IDs Tmax.contracted | cut -f1 > Tmax.contracted.Sig.treebranch  #85

# Functional Enrichment
# Annotation of Tmax with emapper
nohup emapper.py -m diamond --output_dir $outputdir -o Tmax_ReAnn --cpu 40 -i $input --data_dir $db&
cat Tmax_ReAnn.emapper.annotations | sed '/^#/d' > Tmax_ReAnn.clean.annotations

#get the triMax most_expanded gene IDs
cat Tmax_most.expanded.OG2gene.tab | cut -f2 > Tmax_most.expanded_IDs

# assign 1 when most_expanded gene IDs are presented
# 114 genes with GO annotation from 2156 genes (while there are 13539 genes with GO Annotation in the total )
grep -Fwf Tmax_most.expanded_IDs triMax_gene2GO_IDs  > triMax_most_expanded_GOs_IDs

# Then run the GO_MWU script
```

OrthoFinder PipeLine

Data


| Species Code | Species Name          | Original File Name                           | New File Name     |
|--------------|-----------------------|---------------------------------------------|-------------------|
| 6454_0       | Haliotis rufescens    | GCF_003343065.1_H.ruf_v1.0_protein.faa.gz | halRuf.faa.gz     |
| 6500_0       | Aplysia californica   | GCF_000002075.1_AplCal3.0_protein.faa.gz   | aplCal.faa.gz     |
| 6526_0       | Biomphalaria glabrata | GCF_000457365.1_ASM45736v1_protein.faa.gz  | bioGla.faa.gz     |
| 6565_0       | Crassostrea virginica | GCF_002022765.2_C_virginica-3.0_protein.faa.gz | craVir.faa.gz |
| 6573_0       | Mizuhopecten yessoensis | GCF_002113885.1_ASM211388v2_protein.faa.gz | mizYes.faa.gz     |
| 6579_0       | Pecten maximus        | GCF_902652985.1_xPecMax1.1_protein.faa.gz  | pecMax.faa.gz     |
| 6596_0       | Mercenaria mercenaria | GCF_021730395.1_MADL_Memer_1_protein.faa.gz | merMer.faa.gz     |
| 29159_0      | Crassostrea gigas     | GCF_902806645.1_cgigas_uk_roslin_v1_protein.faa.gz | craGig.faa.gz |
| 36100_0      | Haliotis rubra        | GCF_003918875.1_ASM391887v1_protein.faa.gz  | halRub.faa.gz     |
| 37623_0      | Ostrea edulis         | GCF_947568905.1_xbOstEdul1.1_protein.faa.gz | ostEdu.faa.gz     |
| 37653_0      | Octopus bimaculoides  | GCF_001194135.2_ASM119413v2_protein.faa.gz | octBim.faa.gz     |
| 225164_0     | Lottia gigantea       | GCF_000327385.1_Helro1_protein.faa.gz     | lotGig.faa.gz     |
| 400727_0     | Pomacea canaliculata  | GCF_003073045.1_ASM307304v1_protein.faa.gz | pomCan.faa.gz     |
| 1735272_0    | Gigantopelta aegis    | GCF_016097555.1_Gae_host_genome_protein.faa.gz | gigAeg.faa.gz |
| 2607531_0    | Octopus sinensis      | GCF_006345805.1_ASM634580v1_protein.faa.gz | octSin.faa.gz     |

```bash
# run OrthoFinder
nohup orthofinder -f Separate_Genomes_NCBI &
# Find Orthogroup IDs specific to triMax and save them to a file
# /home/ruiqi/ruiqi_data/GiantClamGenome/OrthoDB/Separate_Genomes_NCBI/OrthoFinder/Results_Jan06/Orthogroups
awk 'NR > 1 && $NF == $15 {print $1}' Orthogroups.GeneCount.tsv > triMax_assigned_specific_Ortho_ID

#combine "unassigned genes from triMax
sed '1d' Orthogroups_UnassignedGenes.tsv | cut -f1,15 | grep "Tmax"  > triMax_Unassigned_OG2Gene
cat triMax_Unassigned_OG2Gene triMax_assigned_specific_OG2Gene > triMax_specific_OG2Gene.txt
sed -i 's/Tmax_Tmax/Tmax/g' triMax_specific_OG2Gene.txt
```

```{r}
# Annotation and Enrichment
############################################################################################################
# Part 1 - add annotation to all Tmax specific gene families
############################################################################################################
OG2gene <- read.table(file="triMax_specific_OG2Gene.txt", sep = "\t")

colnames(OG2gene) <- c("familyID","Gene")
Annot <- read.delim(file = "Tmax_ReAnn.clean.annotations", sep = "\t")
colnames(Annot) <- c("Gene","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","best_tax_level","Preferred_name",
                     "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","taxonomic scope","eggNOG OGs",
                     "best eggNOG OG","COG","eggNOG free text desc.")

Tmax_annot <- merge(x=OG2gene,y=Annot,
                        by="Gene", all.x=TRUE)
sorted_Tmax_annot <- Tmax_annot[order(Tmax_annot$familyID), ]

write.table(sorted_Tmax_annot, file = "Tmax_specific_OG_Annot.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
############################################################################################################
# Part 2 - plot bar graph for description for each gene family, mark the description that exceeds 50% in red
############################################################################################################

sorted_Tmax_annot<-read.delim(file = "triMax_assigned_specific_OG_Annot.tsv", sep = "\t", header = TRUE)
unique_og_values <- unique(sorted_Tmax_annot$familyID)
# Replace empty strings with NA in the "eggNOG free text desc." column
sorted_Tmax_annot$eggNOG.free.text.desc.[sorted_Tmax_annot$eggNOG.free.text.desc. == ""] <- NA

# Reform if the description is too long (more than 10 words, then only keep 5 and add ...)
# Function to truncate and add "..."
truncate_and_ellipsis <- function(text, max_words = 10) {
  words <- unlist(strsplit(text, " "))
  if (length(words) > max_words) {
    truncated_text <- paste(words[1:5], collapse = " ")
    return(paste(truncated_text, "..."))
  }
  return(text)
}

# Apply the function to the specified column
sorted_Tmax_annot$eggNOG.free.text.desc. <- sapply(
  sorted_Tmax_annot$eggNOG.free.text.desc.,
  truncate_and_ellipsis
)

# Plot functions
for (og_value in unique_og_values) {
  # Filter the dataframe for the specific familyID
  family_df <- sorted_Tmax_annot %>%
    filter(familyID == og_value)
  # Count the frequency of each description
  description_counts <- family_df %>%
    group_by(eggNOG.free.text.desc.) %>%
    summarize(count = n())
  # Calculate the percentage for each description
  description_counts <- description_counts %>%
    mutate(percentage = (count / sum(count)) * 100)
  # Create a horizontal bar graph with red highlights for descriptions > 50%
  gg <- ggplot(description_counts, aes(x = percentage, y = reorder(eggNOG.free.text.desc., -percentage), fill = percentage > 50)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("lightgrey", "lightcoral"), guide = FALSE) +
    labs(title = paste("Frequency of Descriptions for  Gene Family", og_value),
         x = "Percentage", y = "Description") +
    theme(axis.text.y = element_text(hjust = 1))
  # Save the plot with the name "[OG value]_box_plot.png"
  ggsave(paste0("./GeneFamily_description_bargraph/", og_value, "_description_bargraph.png"), gg, width = 16, height = 9, units = "in")
}

############################################################################################################
# Part 3 - summarize description for each gene family, use the second most frequent description if the most one is NA
############################################################################################################
sorted_Tmax_annot<-read.delim(file = "triMax_assigned_specific_OG_Annot.tsv", sep = "\t", header = TRUE)
sorted_Tmax_annot <- sorted_Tmax_annot %>% select(Gene,familyID,eggNOG.free.text.desc.)
# change space to NA
sorted_Tmax_annot <- sorted_Tmax_annot %>% mutate_all(na_if,"")

# Find the most common value in function desc.
# If the most common one is NA, use the second common one
# If all the values are NA, just use NA
sorted_Tmax_annot <- sorted_Tmax_annot %>%
  group_by(familyID) %>%
  mutate(most_common = {
    # Count the occurrences of each value, excluding NA
    counts <- table(na.omit(eggNOG.free.text.desc.))

    if (length(counts) == 0) {
      # If all values are NA, return NA
      NA
    } else {
      # Get the name of the value with maximum count
      names(sort(counts, decreasing = TRUE))[1]
    }
  }) %>%
  ungroup()

# Create Tmax_summary with familyID, Genes, and most_common
Tmax_summary <- sorted_Tmax_annot %>%
  group_by(familyID) %>%
  summarise(
    Genes = n(),  # Count of genes in each familyID
    most_common = first(most_common)  # The most common value from 'most_common' column
  )

# Calculate Frequency for each familyID independently
Tmax_summary <- Tmax_summary %>%
  rowwise() %>%
  mutate(Frequency = {
    common_value <- most_common
    total_genes <- Genes
    if (is.na(common_value)) {
      NA_real_  # If most_common is NA, set Frequency as NA
    } else {
      common_count <- sum(sorted_Tmax_annot$familyID == familyID & sorted_Tmax_annot$eggNOG.free.text.desc. == common_value, na.rm = TRUE)
      (common_count / total_genes) * 100
    }
  })
############################################################################################################
# Part 4 - summarize function of un-assigned genes
############################################################################################################
OG2gene <- read.table(file="triMax_Unassigned_OG2Gene", sep = "\t")

colnames(OG2gene) <- c("familyID","Gene")
Annot <- read.delim(file = "Tmax_ReAnn.clean.annotations", sep = "\t")
colnames(Annot) <- c("Gene","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","best_tax_level","Preferred_name",
                     "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","taxonomic scope","eggNOG OGs",
                     "best eggNOG OG","COG","eggNOG free text desc.")

Tmax_unasnd_annot <- merge(x=OG2gene,y=Annot,
                    by="Gene", all.x=TRUE)
sorted_Tmax_unasnd_annot <- Tmax_unasnd_annot[order(Tmax_unasnd_annot$familyID), ]
sorted_Tmax_unasnd_annot <- sorted_Tmax_unasnd_annot %>% mutate_if(is.character, na_if, "")
#write.table(sorted_Tmax_unasnd_annot, file = "Tmax_unasnd_OG_Annot.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

################### summarize function of un-assigned genes
sorted_Tmax_unasnd_annot <-read.delim(file = "Tmax_unasnd_OG_Annot.tsv", sep = "\t", header = TRUE)

Tmax_unasnd_summary <- sorted_Tmax_unasnd_annot %>%
  group_by(eggNOG.free.text.desc.) %>%
  summarise(
    Count = n()
  ) %>%
  ungroup() %>%
  arrange(desc(Count))
# Non-TE and Non-NA
Tmax_unasnd_summary_nonTE_nonNA  <- Tmax_unasnd_summary %>%
  filter(!eggNOG.free.text.desc. %in% exclude_values, !is.na(eggNOG.free.text.desc.))

write.table(Tmax_unasnd_summary, file = "Tmax_UNassigned_Annot_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

############################################################################################################
# Part 5 - GO MWU enrichment
############################################################################################################
input="triMax_specific_GOs.csv"
goAnnotations="triMax_gene2GO.txt"
goDatabase="go.obo"
source("gomwu.functions.R")
go=c("MF", "CC", "BP")
for (goDivision in go){
	gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl",
               largest=0.1,
               smallest=5,
               clusterCutHeight=0.25,
    )
}
```


## 2. Phylogenomics

```bash
# grep each OG and save to gene Id to a file
while read line
 do
 grep -w $line Tmax_Mollusca_OG2gene.tab | cut -f2 > ./SingCopyOG_14_genelist/${line}_genelist.txt
 done < Mollusca_SingCopyOG_14_OG
# Retrieve sequences
# Retrieve mollusca specific genes
for file in odb11v_cleanheader.part*.fasta.cidx
do
  cat Mollusca_gene_ID | cdbyank $file > ${file}_Mollusca_gene_seq.fasta
done
cat *_Mollusca_gene_seq.fasta > Mollusca_gene_OrthoDB.fasta
# Retrieve seqs from each OG and save to separate fasta files
for file in ./SingCopyOG_14_genelist/*_genelist.txt
do
  filename=`basename $file`
  name=${filename%_*}
  echo $name
  cat $file | cdbyank Tmax_Mollusca.fasta.cidx > ./SingCopyOG_14_seqs/${name}.fasta
done
# Alignment, Trim, Phylogeny Inference
align each gene separately,
while read name;
do mafft ${name}.fasta > $name.output
 trimal -in $name.output -out ${name}_trimmed.output -automated1
done < Mollusca_SingCopyOG_14_OG

./catfasta2phyml/catfasta2phyml.pl --concatenate --verbose --fasta *_trimmed.output > Tmax_Mollusca_SingCopyOG_14.fasta
echo "alignment done!"

# raxml
run1=all1
run2=all2
run=Tmax_Mollusca_SingCopyOG_14
raxmlHPC-PTHREADS -T 36 -p 12345 -N 20 -s Tmax_Mollusca_SingCopyOG_14.fasta -m PROTGAMMAWAG -n $run1
raxmlHPC-PTHREADS -T 36 -p 12345 -b 12345 -s Tmax_Mollusca_SingCopyOG_14.fasta -N 100 -m PROTGAMMAWAG -n $run2
raxmlHPC-PTHREADS -T 36 -f b -t RAxML_bestTree.${run1}  -z RAxML_bootstrap.${run2} -m PROTGAMMAWAG -n $run

```
## 3. Repeat Annotation

Pipeline demo with Cerastoderma edule genome

```bash
# # RepeatModeler
BuildDatabase -name cerEdu C_edule.fna

# Annotate Repeats
git clone https://github.com/darencard/GenomeAnnotation.git

./GenomeAnnotation/repclassifier -t 36 -d Mollusca -u cerEdu_families.prefix.fa.unknown \
-k cerEdu_families.prefix.fa.known -a cerEdu_families.prefix.fa.known \
-o round-1_RepbaseMollusca-Self
#986 unknowns after round 1

# round 2 with know library

./GenomeAnnotation/repclassifier -t 36 -u round-1_RepbaseMollusca-Self/round-1_RepbaseMollusca-Self.unknown \
-k round-1_RepbaseMollusca-Self/round-1_RepbaseMollusca-Self.known \
-a round-1_RepbaseMollusca-Self/round-1_RepbaseMollusca-Self.known -o round-2_Self
# 825 unknowns after round 2

# round 3 with know library

./GenomeAnnotation/repclassifier -t 36 -u round-2_Self/round-2_Self.unknown \
-k round-2_Self/round-2_Self.known \
-a round-2_Self/round-2_Self.known -o round-3_Self
# 765 unknowns after round 3

# round 4 with know library

./GenomeAnnotation/repclassifier -t 36 -u round-3_Self/round-3_Self.unknown \
-k round-3_Self/round-3_Self.known \
-a round-3_Self/round-3_Self.known -o round-4_Self
# 742 unknowns after round 4

# round 5 with know library

./GenomeAnnotation/repclassifier -t 36 -u round-4_Self/round-4_Self.unknown \
-k round-4_Self/round-4_Self.known \
-a round-4_Self/round-4_Self.known -o round-5_Self
# 736 unknowns after round 5

cat round-5_Self.unknown round-5_Self.known > C_edule_denovo_annot.fa

# round 1: annotate/mask simple repeats
nohup /home/ruiqi/ruiqi_data/2023GiantClam/RepeatMaskerInstall/RepeatMasker/RepeatMasker -pa 48 -a -e ncbi -dir 01_simple_out -noint -xsmall C_edule.fna &

mv nohup.out logs/01_simplemask.log

# rename outputs for 01_simple_out
for file in *; do mv -v "$file" "${file/C_edule/cerEdu_simple_mask}"; done
mv cerEdu_simple_mask.fna.masked cerEdu_simple_mask.fna


# round 2: annotate/mask Tetrapoda elements sourced from Repbase using output from 1st round of RepeatMasker
nohup /home/ruiqi/ruiqi_data/2023GiantClam/RepeatMaskerInstall/RepeatMasker/RepeatMasker -pa 48 -a -e ncbi -dir 02_mollusca_out -nolow -species Mollusca 01_simple_out/cerEdu_simple_mask.fna &

mv nohup.out logs/02_molluscaMask.log

# round 2:rename outputs
for file in *; do mv -v "$file" "${file/simple_mask/mollusca_mask}"; done
mv cerEdu_mollusca_mask.fna.masked cerEdu_mollusca_mask.masked.fasta

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
nohup /home/ruiqi/ruiqi_data/2023GiantClam/RepeatMaskerInstall/RepeatMasker/RepeatMasker -pa 48 -a -e ncbi -dir 03_known_out -nolow -lib round-5_Self.known \
02_mollusca_out/cerEdu_mollusca_mask.masked.fasta &

# round 3:rename outputs
for file in *; do mv -v "$file" "${file/mollusca_mask/known_mask}"; done
mv cerEdu_known_mask.masked.fasta.masked cerEdu_known_mask.masked.fasta

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output from 3nd round of RepeatMasker
nohup /home/ruiqi/ruiqi_data/2023GiantClam/RepeatMaskerInstall/RepeatMasker/RepeatMasker -pa 48 -a -e ncbi -dir 04_unknown_out -nolow -lib round-5_Self.unknown \
03_known_out/cerEdu_known_mask.masked.fasta &
# round 4:rename outputs
for file in *; do mv -v "$file" "${file/known_mask/unknown_mask}"; done
mv cerEdu_unknown_mask.masked.fasta.masked cerEdu_unknown_mask.masked.fasta

# create directory for full results
mkdir -p 05_full_out

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/cerEdu_simple_mask.fna.cat.gz 02_mollusca_out/cerEdu_mollusca_mask.fna.cat.gz 03_known_out/cerEdu_known_mask.masked.fasta.cat.gz 04_unknown_out/cerEdu_unknown_mask.masked.fasta.cat.gz

cat 01_simple_out/cerEdu_simple_mask.fna.cat.gz 02_mollusca_out/cerEdu_mollusca_mask.fna.cat.gz 03_known_out/cerEdu_known_mask.masked.fasta.cat.gz 04_unknown_out/cerEdu_unknown_mask.masked.fasta.cat.gz > 05_full_out/cerEdu_full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/cerEdu_simple_mask.fna.out \
<(cat 02_mollusca_out/cerEdu_mollusca_mask.fna.out | tail -n +4) \
<(cat 03_known_out/cerEdu_known_mask.masked.fasta.out | tail -n +4) \
<(cat 04_unknown_out/cerEdu_unknown_mask.masked.fasta.out | tail -n +4) \
> 05_full_out/cerEdu.full_mask.out


# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/cerEdu_simple_mask.fna.align \
02_mollusca_out/cerEdu_mollusca_mask.fna.align \
03_known_out/cerEdu_known_mask.masked.fasta.align \
04_unknown_out/cerEdu_unknown_mask.masked.fasta.align \
> 05_full_out/cerEdu.full_mask.align

# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
nohup /home/ruiqi/ruiqi_data/2023GiantClam/RepeatMaskerInstall/RepeatMasker/ProcessRepeats -a -species Mollusca 05_full_out/cerEdu_full_mask.cat.gz &

mv nohup.out logs/05_fullmask.log

# Insertion time
nohup ../../Parsing-RepeatMasker-Outputs/parseRM.pl -i ../cerEdu_full_mask.align -l 120,5 -m 0.00109 -v&
```

## 4. Circos Plot

```bash
# calculate gene density
cat Tmax_merged_longestisoform.gff | cut -f 1,4,5 > genes.bed
cut -d ' ' -f 3,6 chromosome.txt | tr ' ' '\t' > chromosome.genome
bedtools makewindows -g chromosome.genome -w 1000000 > chromosome.windows
bedtools coverage -a chromosome.windows -b genes.bed | cut -f 1-4 > genes_num.txt

# calculate repeat density
cat PO1429_Tridacna_maxima.RepeatMasked.gff | cut -f 1,4,5 > repeats.bed
bedtools coverage -a chromosome.windows -b repeats.bed | cut -f 1-4 > repeats_num.txt

# calculate GC content
bedtools nuc -fi $genome -bed chromosome.windows > gc_raw.txt

# Configuration file
karyotype = chromosome.txt
chromosomes_units=1000000
<ideogram>

<spacing>
default = 0.005r
</spacing>

radius    = 0.75r
thickness = 25p
fill      = yes

show_bands = yes
fill_bands = yes

stroke_color     = dgrey
stroke_thickness = 2p

show_label       = yes
label_font       = default
label_radius     = 1.05r
label_size       = 30
label_parallel   = yes


show_ticks          = yes
show_tick_labels    = yes

</ideogram>

show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d

<tick>
spacing        = 10u
size           = 10p
show_label = yes
label_offset = 8p
label_size = 15p
format = %d
</tick>

<tick>
spacing        = 5u
size           = 5p
show_label = no
</tick>

</ticks>

<plots>

<plot>
type    = heatmap
file    = GCrate.txt
color   = purd-9-seq
r1      = 0.85r
r0      = 0.80r
</plot>

<plot>
type      = heatmap
file = repeats_num.txt
color   = blues-9-seq
r1      = 0.91r
r0      = 0.86r
</plot>
<plot>
type      = line
thickness = 1
max_gap = 1u
file = repeats_num.txt
color   = black
r1      = 0.91r
r0      = 0.86r
</plot>

<plot>
type = heatmap
file = genes_num.txt
color = purples-9-seq
r1   = 0.99r
r0   = 0.92r
</plot>

<plot>
type = line
thickness = 1
max_gap = 1u
file = genes_num.txt
color = black
r1   = 0.99r
r0   = 0.92r
</plot>


</plots>

<image>

# Run the program
circos -conf circos.conf
```
# 5. Targeted Gene Finder

Targeted Gene Finder demo with Paradox CDX gene

Find sequence: [Nucula tumidula homeobox parahox cdx (Cdx) mRNA, complete cds](https://www.ncbi.nlm.nih.gov/nuccore/KX365143.1)

```bash
# base Environment
# makedb
makeblastdb -in C_edule.fna -dbtype nucl -parse_seqids -out C_edule
makeblastdb -in Tmax_ReAnn_dna.fasta -dbtype nucl -parse_seqids -out T_maxima_Annot
makeblastdb -in GCF_902652985.1_xPecMax1.1_cds_from_genomic.fna -dbtype nucl -parse_seqids -out PecMax

## tblastn
tblastx -query Paradox_cdx.fasta \
    -db C_edule \
    -outfmt 6 -evalue 1e-5 -num_threads 30 >  Paradox_cdx_C_edule.outfmt6

tblastx -query Paradox_cdx.fasta \
    -db T_maxima_Annot \
    -outfmt 6 -evalue 1e-5 -num_threads 30 > Paradox_cdx_T_maxima_Annot.outfmt6
```

# 6. Transcriptome

One tissue each (mantle,  adductor  muscle,  go-nads,  gills,  byssus,  visceral  mass,  and  kidney) from [this study](https://www.ncbi.nlm.nih.gov/nuccore/1897846671)

Adductor muscle 2 (TM_A_N02) is from Li et al. 2020.

Mantle mRNA available 29/53/65 day at control (29.2°C; ambient (condition N) and elevated (30.7°C) (condition E) temperature

| Tissue | Sample Name | Original File Names | SRR number |
|--------|------------|----------------------|------------|
|Adductor Muscle|TM_A_N01|Tmax_B_S4_R1_001_val_1.fq.gz, Tmax_B_S4_R2_001_val_2.fq.gz|NA|
|Adductor Muscle|TM_A_N02|SRR8217859_1.fastq,SRR8217859_2.fastq|SRR8217859|
|Adductor Muscle|TM_A_N03|SRR12486110_1.fastq, SRR12486110_2.fastq|SRR12486110|
|Mantle|TM_M_N01|SRR12486156_1.fastq, SRR12486156_2.fastq|SRR12486156|
|Mantle|TM_M_N02|SRR12486155_1.fastq, SRR12486155_2.fastq|SRR12486155|
|Mantle|TM_M_N03|SRR12486141_1.fastq, SRR12486141_2.fastq|SRR12486141|
|Mantle|TM_M_N04|SRR12486109_1.fastq, SRR12486109_2.fastq|SRR12486109|
|Mantle|TM_M_E01|SRR12486154_1.fastq, SRR12486154_2.fastq|SRR12486154|
|Mantle|TM_M_E02|SRR12486147_1.fastq, SRR12486147_2.fastq|SRR12486147|
|Mantle|TM_M_E03|SRR12486140_1.fastq, SRR12486140_2.fastq|SRR12486140|

```bash
# Trim the reads following default setting in Trinity

trimmomatic PE -threads 20 TM_A_N01_1.fastq.gz TM_A_N01_2.fastq.gz \
TM_A_N01_1.trim.fq.gz TM_A_N01_1.untrim.fq.gz \
TM_A_N01_2.trim.fq.gz TM_A_N01_2.untrim.fq.gz \
ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25

#Use a loop
#!/bin/bash
for file in $(ls *_1.fastq.gz)
do
  echo $file
  file="${file%%_1.*}"
  echo $file
  trimmomatic PE -threads 36 ${file}_1.fastq.gz ${file}_2.fastq.gz \
  TrimmedReads/${file}_1.trim.fq.gz Untrimmed/${file}_1.untrim.fq.gz \
  TrimmedReads/${file}_2.trim.fq.gz Untrimmed/${file}_2.untrim.fq.gz \
  ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
done
```

### 3. QC and Quantification

```bash
# Prepare the sample file
for file in *_1.trim.fq.gz; do
    prefix="${file%%_1.trim.fq.gz}"
    path="/home/ruiqi/ruiqi_data/GiantClamGenome/Transcriptome_GiantClam/TrimmedReads/"
    col1="${prefix:3:1}${prefix:5:1}"
    col2="${prefix:0:8}"
    col3="$path${file}"
    col4="$path${prefix}_2.trim.fq.gz"
    echo -e "$col1\t$col2\t$col3\t$col4"
done > TM_samples.txt

# env: trinity
# Use align_and_estimate_abundance.pl to quantify unigenes
nohup /home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/align_and_estimate_abundance.pl \
--transcripts Tmax_ReAnn_dna.fasta \
--seqType fq \
--samples_file TM_samples.txt \
--est_method kallisto \
--prep_reference \
--kallisto_add_opts "-t 36" &

# get the abundace file list
find . -name "abundance.tsv" | sort | tee TM_abundance_files.list

# abuandance to matrix
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/abundance_estimates_to_matrix.pl --est_method kallisto \
--gene_trans_map none \
--out_prefix TM_all \
--name_sample_by_basedir \
--quant_files TM_abundance_files.list

# Counting Numbers of Expressed Transcripts or Genes
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
TM_all.isoform.TPM.not_cross_norm | tee TM_all.isoform_matrix.TPM.not_cross_norm.counts_by_min_TPM

# DE
# Use TPM > 5 as cut off for DE analyses
#Retained 16119 / 46469 = 34.69% of total transcripts.
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/filter_low_expr_transcripts.pl --matrix TM_all_afterQC.isoform.TPM.not_cross_norm --transcripts Tmax_ReAnn_dna.fasta --min_expr_any 5 > TPM5_Tmax_ReAnn_dna.fasta
# clean header
sed -i '/^>/ s/ .*//' TPM5_Tmax_ReAnn_dna.fasta
#get the gene list (TPM>5)
grep ">" TPM5_Tmax_ReAnn_dna.fasta | sed 's/>//' > TPM5_Tmax.geneID
#get the headers
head -1 TM_all_afterQC.isoform.counts.matrix >TPM5_TM_all.isoform.counts.matrix
# get the count matrix on the isoform list
grep -wFf TPM5_Tmax.geneID TM_all_afterQC.isoform.counts.matrix >> TPM5_TM_all.isoform.counts.matrix

## Step 1: Identifying DEs
# TPM can't be in the file name, rename it to
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/run_DE_analysis.pl \
--matrix expr5_TM_all.isoform.counts.matrix \
--samples_file TM_DE_samples.txt \
--contrasts TM_contrasts.txt \
--output . \
--method DESeq2

## Step 2: Extracting and Clustering DEs
## please note that we are using the TMM.EXPR.matrix
# filtering it first
#get the headers
head -1 TM_all_afterQC.isoform.TMM.EXPR.matrix > expr5_TM_all_afterQC.isoform.TMM.EXPR.matrix
# get the count matrix on the isoform list
grep -wFf TPM5_Tmax.geneID TM_all_afterQC.isoform.TMM.EXPR.matrix >> expr5_TM_all_afterQC.isoform.TMM.EXPR.matrix

/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/analyze_diff_expr.pl \
--matrix expr5_TM_all_afterQC.isoform.TMM.EXPR.matrix \
--samples TM_DE_samples.txt \
--output TM_P${Pvalue}_C${C} \
--P 5e-2 \
--C 1
#${prefix}.sampleA_vs_sampleB.voom.DE_results.P0.001_C2.sampleA-UP.subset : the expression matrix subset for features up-regulated in sampleA
# ${prefix}.sampleA_vs_sampleB.voom.DE_results.P0.001_C2.sampleB-UP.subset : the expression matrix subset for features up-regulated in sampleB

# Enrichment

setwd("/home/ruiqi/ruiqi_data/GiantClamGenome/Transcriptome_GiantClam/GO_MWU")

goDatabase="go.obo"
source("gomwu.functions.R")
###########
# Tmax
###########
goAnnotations="expr5_Tmax_gene2go.tab"
# Input
go=c("MF", "CC", "BP")
input_files=c("Tmax_MN_vs_ME_lfc.csv","Tmax_MN_vs_AN_lfc.csv")
# Calculating stats
for (goDivision in go){
  for (input in input_files){
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl",
               largest=0.1,
               smallest=5,
               clusterCutHeight=0.25,
    )
  }
}

library(tidyverse)
 go=c("MF", "CC", "BP")
 input_files=c("Tmax_MN_vs_AN_lfc.csv","Tmax_MN_vs_ME_lfc.csv")
 # Calculating stats)
 for (goDivision in go){
   for (input in input_files){
     mwu_data <- read.table(paste0("MWU_",goDivision,"_",input), sep = " ", header=TRUE) %>%
       filter(p.adj<0.05)
     mwu_data$change[mwu_data$delta.rank > 0] <- "Up"
     mwu_data$change[mwu_data$delta.rank < 0] <- "Down"
     mwu_up <- mwu_data %>% filter (change == "Up") %>% select(c("term", "name"))
     write.table(mwu_up, file=paste0("./GOs_GO_MWU/GOs_UP_",substring(input, 1, 13), "_", goDivision, ".tsv"), quote = FALSE, sep="\t",row.names = FALSE, col.names = TRUE)
     mwu_down <- mwu_data %>% filter (change == "Down") %>% select(c("term", "name"))
     write.table(mwu_down, file=paste0("./GOs_GO_MWU/GOs_DOWN_",substring(input, 1, 13), "_", goDivision, ".tsv"), quote = FALSE, sep="\t",row.names = FALSE, col.names = TRUE)
   }
 }
```

# from GO_MWU to genes
```bash
# get all the microtubule/cilium/ciliary/cilia related GO terms in the down-regulated genes in the MN
for file in GOs_DOWN_Tmax*; do
  echo $file
  grep -E 'microtubule|cilium|ciliary|cilia' $file >> cilia_DOWN_MN.txt
done

# if there are multiple GO terms in one line, split into several lines
awk -F'\t' '{
  split($1, terms, ";")
  for (i in terms) {
    print terms[i] "\t" $2
  }
}' cilia_DOWN_MN.txt > cilia_DOWN_MN_1.txt
# Keep the unique GO terms
cat cilia_DOWN_MN_1.txt | sort -u > cilia_DOWN_MN_unique.txt
# get the GO terms
cut -f1 cilia_DOWN_MN_unique.txt > cilia_DOWN_MN_GO_unique.txt



# get the list of Down-Regulated gene in MN
cat expr5_TM_all.isoform.counts.matrix.MN_vs_AN.DESeq2.DE_results.P5e-2_C1.AN-UP.subset | sed '1d' | cut -f1 > MN_vs_AN_Down_genelist

cat expr5_TM_all.isoform.counts.matrix.MN_vs_ME.DESeq2.DE_results.P5e-2_C1.ME-UP.subset | sed '1d' | cut -f1 >> MN_vs_AN_Down_genelist
# get the unique genes
sort -u  MN_vs_AN_Down_genelist > MN_vs_AN_Down_unique_genelist


###  get the list of genes that annotated with the cilia GO terms, and down regulated in MN
grep -Fwf cilia_DOWN_MN_GO_unique.txt Tmax_ReAnn.clean.annotations | grep -Fwf TM_DE_subset/MN_vs_AN_Down_unique_genelist > cilia_DOWN_MN_genes.txt

cut -f1,6,22 cilia_DOWN_MN_genes.txt cilia_DOWN_MN_genes_summary.txt
```
