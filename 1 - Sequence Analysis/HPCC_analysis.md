## Analysis of 16S rRNA gene Illumina sequencing data for 
## Mimulus Reciprocal Transplant Project.


###############################################################

## RAW DATA FOR RECIPROCAL TRANSPLANT CAN BE FOUND AT:
/mnt/research/ShadeLab/Sequence/raw_sequence/PRI_Kearns/mimulus

## COPY THE RAW RECIPROCAL TRANSPLANT DATA INTO WORKING SPACE
cp *.fastq /mnt/research/ShadeLab/WorkingSpace/Bowsher/Mimulus_RecipTrans_16S/

###############################################################

## RAW DATA FOR THE BULK SOILS CAN BE FOUND AT
/mnt/research/ShadeLab/Sequence/raw_sequence/Bowsher/20180525_16S-V4_PE250

## COPY THE RAW BULK SOIL DATA INTO WORKING SPACE
cp Mim* /mnt/research/ShadeLab/WorkingSpace/Bowsher/Mimulus_RecipTrans_16S/

###############################################################

## MOVE TO THE WORKING SPACE AND GUNZIP ALL SEQUENCE FILES
cd /mnt/research/ShadeLab/WorkingSpace/Bowsher/Mimulus_RecipTrans_16S/

gunzip *.gz

## REMOVE SAMPLE 201, which clearly had sequencing errors.
rm Mim201*

# CHECK THAT CORRECT NUMBER OF FILES ARE PRESENT (SHOULD BE 98)
ls -1 | wc -l

# MOVE THIS RAW INPUT DATA TO A NEW DIRECTORY
mkdir InputData
mv *.fastq InputData/

## MERGE PAIRED END READS
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs InputData/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -threads 1 -fastq_minmergelen 250 -fastq_maxmergelen 300

## TRIM PRIMER-BINDING REGIONS FROM SEQUENCES
## THIS TRIMS PRIMER1, AND REVERSE COMPLEMENT OF PRIMER2.
~/.local/bin/cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o cut_merged.fastq mergedfastq/merged.fq > cut_adpt_results.txt

## FILTER THE MERGED AND TRIMMED SEQUENCES
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter cut_merged.fastq -fastq_maxee 1 -fastqout filtered_merged.fq -threads 1

## DEREPLICATE SEQUENCES (FIND UNIQUE SEQUENCES)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_merged.fq -threads 1 -fastqout uniques_merged.fastq -sizeout

## REMOVE SINGLETONS
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize uniques_merged.fastq -fastqout nosigs_uniques_merged.fastq -minsize 2

## CLUSTER DEREPLICATED SEQS INTO ZOTUS (ZERO-OTUS)
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 nosigs_uniques_merged.fastq -zotus rep_set.fna -tabbedout rep_set_report.txt

#############################################################

#####WARNING: SHIFTED SEQUENCES DETECTED.
#####Not a problem though. I checked it using Edgar's suggested #####'usearch_global' method: #####https://www.drive5.com/usearch/manual/otu_qc_offset.html 

#############################################################

## CAPITALIZE ALL Zotu NAMES TO ZOTU (NECESSARY DOWNSTREAM)
sed -i 's/Zotu/ZOTU/g' rep_set.fna

## MAP PRE-DEREPLICATED SEQUENCES TO ZOTUS (MAKE ZOTU TABLE).
## SWITCH TO USEARCH9 FOR DOWNSTREAM COMPATIBILITY WITH QIIME.
 
/mnt/research/rdp/public/thirdParty/usearch9.2.64_i86linux64 -usearch_global cut_merged.fastq -threads 1 -db rep_set.fna -strand plus -id 0.97 -uc zOTU_map.uc -otutabout zOTU_table.txt -biomout zOTU_jsn.biom

###############################################################

## ACTIVATE QIIME ENVIRONMENT (THIS IS QIIME1.9.0 DESPITE CODE)
source activate qiime1.8.0

###############################################################

## ALIGN ZOTU REPRESENTATIVE SEQUENCES

#DOWNLOAD MUSCLE 3.8.1 (https://www.drive5.com/muscle/downloads.htm

#DECOMPRESS AND EXTRACT THE MUSCLE PROGRAM
tar -xvzf muscle3.8.31_i86linux64.tar.gz

#FINALLY, PERFORM THE ALIGNMENT
~/muscle3.8.31_i86linux64 -in rep_set.fna -out alignment.fasta -maxiters 2 -diags1

###############################################################

## ASSIGN TAXONOMY TO REPRESENTATIVE SEQS USING SILVA DATABASE
assign_taxonomy.py -i rep_set.fna -t /mnt/research/ShadeLab/SharedResources/SILVA123_QIIME_release/taxonomy/16S_only/97/majority_taxonomy_7_levels.txt -r /mnt/research/ShadeLab/SharedResources/SILVA123_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -o taxonomy

## ADD TAXONOMY ASSIGNMENTS TO THE .BIOM ZOTU TABLE
biom add-metadata -i zOTU_jsn.biom -o zotu_table_tax.biom --observation-metadata-fp=taxonomy/rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy

## FILTER PLASTIDS, MITO, AND UNASSIGNED FROM .BIOM ZOTU TABLE
filter_taxa_from_otu_table.py -i zotu_table_tax.biom -o zotu_table_tax_rmCM.biom -n D_2__Chloroplast,D_4__Mitochondria,Unassigned

## FILTER PLASTIDS, MITO, AND UNASSIGNED FROM REP SEQS FILE
filter_fasta.py -f rep_set.fna -o rep_set_rmCM.fa -b zotu_table_tax_rmCM.biom

## CHECK .BIOM ZOTU TABLE THAT PLASTIDS ARE REMOVED
grep -c D_2__Chloroplast zotu_table_tax.biom
grep -c D_2__Chloroplast zotu_table_tax_rmCM.biom

## CHECK .BIOM ZOTU TABLE THAT MITOCHONDRIA ARE REMOVED
grep -c D_4__Mitochondria zotu_table_tax.biom
grep -c D_4__Mitochondria zotu_table_tax_rmCM.biom

## CHECK .BIOM ZOTU TABLE THAT UNASSIGNED READS ARE REMOVED
grep -c Unassigned zotu_table_tax.biom
grep -c Unassigned zotu_table_tax_rmCM.biom

## CHECK THAT SEQS HAVE BEEN REMOVED FROM REP SEQS FILE
grep ">" rep_set.fna | wc -l  
grep ">" rep_set_rmCM.fa | wc -l

############################################################

## MAKE PHYLOGENY WITH FASTTREE
make_phylogeny.py -i alignment.fasta -o rep_set.tre

## SUMMARIZE THE ZOTU TABLE
biom summarize-table -i zotu_table_tax_rmCM.biom -o zotu_table_summary.txt

## RAREFY ZOTU TABLE TO LOWEST SEQUENCING DEPTH
single_rarefaction.py -d 22354 -o single_rare.biom -i zotu_table_tax_rmCM.biom

## CONVERT RAREFIED ZOTU TABLE TO TXT FORMAT
biom convert -i single_rare.biom -o zotu_table_forR.txt --to-tsv

############################################################

## CALCULATE ALPHA DIVERSITY
alpha_diversity.py -m PD_whole_tree,shannon -i single_rare.biom -o alpha -t rep_set.tre

## MERGE ALPHA DIVERSITY OUTPUT WITH METADATA
add_alpha_to_mapping_file.py -i alpha -m Mimulus_metadata.txt -o Mimulus_metadata_alpha.txt

## CALCULATE PAIRWISE WEIGHTED UNIFRAC DISTANCE (BETA DIVERSITY)
beta_diversity.py -m weighted_unifrac -i single_rare.biom -o beta_div -t rep_set.tre

## CALCULATE AXIS SCORES FOR PCOA BASED ON WEIGHTED UNIFRAC
principal_coordinates.py -i beta_div -o coords
#Note: axis1 explains 45.8%, axis2 explains 13.9%


###############################################################

## ADD ALPHA DIVERSITY VALUES AND WU AXES TO METADATA

# FIRST, SORT THE METADATA NUMERICALLY
tail -n+2 Mimulus_metadata_alpha.txt | sort -k1 >> Mimulus_metadata_alpha_sorted.txt

# REMOVE HEADER INFORMATION FROM WEIGHTED UNIFRAC DATAFILE
grep "Mim" coords/pcoa_weighted_unifrac_single_rare.txt > coords/WUvalues.txt

# SELECT THE FIRST TWO AXES OF THE WEIGHTED UNIFRAC DATAFILE
cut -f1,2,3 coords/WUvalues.txt > coords/WUaxes1and2.txt

# ADD TAB-DELIMITED HEADER TO FILE CONTAINING FIRST TWO AXES
{ printf 'Sample\tWU_Axis1\tWU_Axis2\n'; cat coords/WUaxes1and2.txt; } > coords/WUaxes1and2_named.txt

# VERIFY THAT METADATA AND WU AXIS1AND2 FILE ARE IN SAME ORDER
cut -f1 Mimulus_metadata_alpha_sorted.txt >new1
cut -f1 coords/WUaxes1and2_named.txt >new2  
paste new1 new2

# PASTE WU AXIS1 AND 2 VALUES INTO METADATA FILE
cut -f2,3 coords/WUaxes1and2_named.txt > b.tmp 

paste Mimulus_metadata_alpha_sorted.txt b.tmp > Mimulus_metadata_FULL.txt

###############################################################

## SUMMARIZE TAXONOMIC DATA BY EACH TAXONOMIC LEVEL

summarize_taxa.py -i single_rare.biom -o taxa_sum

################################################################ 

## EXPORT THE FOLLOWING FILES TO R FOR ECOLOGICAL ANALYSIS:

zotu_table_forR.txt

Mimulus_metadata_FULL.txt

taxa_sum/single_rare_L3.txt

beta_div/weighted_unifrac_single_rare.txt

#############################################################











































