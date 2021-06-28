## Single-Cell RNA-seq 
#### Reproduced from Alex K. Shalek, Nature, 2013
### 1. Prepare Files
##### Build STAR index
```shell
#!/usr/bin/bash
#SBATCH -J starindex
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /data/user_05/projects/shalek_nature_2013/data/starindex_mm10 \
--genomeFastaFiles /data/user_05/projects/shalek_nature_2013/data/reference/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile /data/user_05/projects/shalek_nature_2013/data/reference/gencode.vM27.annotation.gtf \
--sjdbOverhang 100
```
(omit job scheduling lines in the following scripts)
##### Build RSEM index
```shell
rsem-prepare-reference --gtf /data/user_05/projects/shalek_nature_2013/data/reference/gencode.vM27.annotation.gtf \
--bowtie \
/data/user_05/projects/shalek_nature_2013/data/reference/GRCm39.primary_assembly.genome.fa \
/data/user_05/projects/shalek_nature_2013/data/ref_mouse
```
##### Build BLAST database
```shell
#go to your working directory
cd '/data/user_05/projects/shalek_nature_2013/data/blast'
#mkdir ref_prok_rep_genomes
cd ref_prok_rep_genomes_v5

#download prokaryotic database
#wget https://ftp.ncbi.nlm.nih.gov/blast/db/ref_prok_rep_genomes.??.tar.gz
#tar zxvf ref_prok_rep_genomes.??.tar.gz

#extract database
#Could play with -outfmt to get easier parsing for a tax map
blastdbcmd -db ref_prok_rep_genomes -entry all > ref_prok_genomes_v5.fna

cd '/data/user_05/projects/shalek_nature_2013/data/blast'
makeblastdb -in /data/user_05/projects/shalek_nature_2013/data/blast/ref_prok_rep_genomes_v5/ref_prok_genomes_v5.fna \
-dbtype nucl -parse_seqids \
-out prok_db_v5/ref_prok_v5
```
### 2. Fastqc
```shell
# first build a sample list
list_sample='/data/user_05/projects/shalek_nature_2013/srrlisttest.txt'

# specify output directory
# fastq file like sample_R1.fastq.gz/sample_R2.fastq.gz
fastq_dir='/data/user_05/projects/shalek_nature_2013/data/fastq'
output_dir='/data/user_05/projects/shalek_nature_2013/output'

## specify Org, option available at present: mm10, hg38, others get error ! and will exit
## have been replaced with 2019 latest version GTF

#Org='mm10'

#RRNA_FA='dir'
#STAR_INDEX='/data/user_05/projects/shalek_nature_2013/data/starindex_mm10'
#RSEM_TRANS='dir'

#test
for sample in $(cat ${list_sample})
do
#### pre-process
echo "create alignment directory"

align_dir=$(echo ${output_dir}/${sample}/align )
qc_dir=$(echo ${output_dir}/${sample}/qc )
mkdir -p ${align_dir}
mkdir -p ${qc_dir}

echo "create symbolic link"
ln -s ${fastq_dir}/${sample}_1.fastq.gz ${align_dir}/${sample}_R1.fastq
ln -s ${fastq_dir}/${sample}_2.fastq.gz ${align_dir}/${sample}_R2.fastq

## fastqc
echo "fastqc"
fastqc \
-t 6 \
${fastq_dir}/${sample}_1.fastq.gz \
${fastq_dir}/${sample}_2.fastq.gz \
-o ${qc_dir} \
-d ${qc_dir}/
done
```
### 3. STAR align
```shell

list_sample='/data/user_05/projects/shalek_nature_2013/srrlisttest.txt'
STAR_INDEX='/data/user_05/projects/shalek_nature_2013/data/starindex_mm10'
output_dir='/data/user_05/projects/shalek_nature_2013/output'
##star-rsem alignment
for sample in $(cat ${list_sample})
do
align_dir=$(echo ${output_dir}/${sample}/align )
echo "start star alignment"
mkdir -p ${align_dir}/STAR
STAR --genomeDir ${STAR_INDEX} \
--runThreadN 8 \
--readFilesIn \
${align_dir}/${sample}_R1.fastq \
${align_dir}/${sample}_R2.fastq \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--limitBAMsortRAM 30000000000 \
--outFileNamePrefix ${align_dir}/STAR/ \
--outSAMstrandField intronMotif \
--quantMode TranscriptomeSAM
done
```
### 4. BamCoverage
```shell
list_sample='/data/user_05/projects/shalek_nature_2013/srrlisttest.txt'

output_dir='/data/user_05/projects/shalek_nature_2013/output'
#bam file process
#index $ bigwig coverage for igv visualization

for sample in $(cat ${list_sample})
do

align_dir=$(echo ${output_dir}/${sample}/align )

samtools index -@ 6	${align_dir}/STAR/Aligned.sortedByCoord.out.bam

echo "Start bamCoverage"
bamCoverage \
-p 6 \
-b $(echo ${output_dir}/${sample}/align )/STAR/Aligned.sortedByCoord.out.bam \
-o $(echo ${output_dir}/${sample}/align )/STAR/Aligned.sortedByCoord.out.bw
done
```
### 5. RSEM
```shell
list_sample='/data/user_05/projects/shalek_nature_2013/srrlisttest.txt'
RSEM_TRANS='/data/user_05/projects/shalek_nature_2013/data/rsem_trans/ref_mouse'
output_dir='/data/user_05/projects/shalek_nature_2013/output'

for sample in $(cat ${list_sample})
do

align_dir=$(echo ${output_dir}/${sample}/align )

#rsem
mkdir -p ${align_dir}/RSEM_STAR
rsem-calculate-expression -p 8 \
--no-bam-output \
--paired-end \
--alignments ${align_dir}/STAR/Aligned.toTranscriptome.out.bam \
${RSEM_TRANS} \
${align_dir}/RSEM_STAR/RSEM_STAR_${sample}
done
```

### 6. file summary
```shell
list_sample='/data/user_05/projects/shalek_nature_2013/SRR_Acc_List.txt'
output_dir='/data/user_05/projects/shalek_nature_2013/output'

###summary
cd ${output_dir}

mkdir bigwig_file
mkdir bam_file
mkdir gene_result

# put all bigwig (for IGV), gene_result from RSEM (keep col 1,5,6,7) into a new dir
# gene_id transcript_id(s) length effective_length expected_count TPM FPKM
for line in $(cat ${list_sample})
do
cp ./${line}/align/STAR/Aligned.sortedByCoord.out.bw ./bigwig_file/${line}.bw
cat ./${line}/align/RSEM_STAR/*.genes.results |gawk 'OFS="\t"{print $1,$5,$6,$7 }' \
> ./gene_result/${line}.genes.results
done

# create link to bamfile of each sample (sort by genome coordinates)
for line in $(cat ${list_sample})
do
ln -s ${output_dir}/${line}/align/STAR/Aligned.sortedByCoord.out.bam bam_file/${line}.bam
ln -s ${output_dir}/${line}/align/STAR/Aligned.sortedByCoord.out.bam.bai \
bam_file/${line}.bam.bai
done
```

### 7. Summary.csv
```shell
list_sample='/data/user_05/projects/shalek_nature_2013/SRR_Acc_List.txt'
output_dir='/data/user_05/projects/shalek_nature_2013/output'
# summary.csv
# sample, total reads, unique mapping reads/ratio, multiple mapping reads/ratio, rRNA, rRNA_r= rRNA/(unique+multiple)
# add tpm>0, tpm>2 gene count
cd ${output_dir}
#echo -e "Sample,TotalReads,Unique,UniqueRatio,Multiple,MultipleRatio,rRNA,rRNARatio,genes(TPM>0),genes(TPM>2)" \
#> summary.csv
echo -e "Sample,TotalReads,Unique,UniqueRatio,Multiple,MultipleRatio,genes(TPM>0),genes(TPM>2)" \
> summary0.csv

for line in $(cat ${list_sample})
do
Reads=$(cat ./${line}/align/STAR/Log.final.out |grep "Number of input reads" |cut -f 2)
Unique=$(cat ./${line}/align/STAR/Log.final.out |grep "Uniquely mapped reads number" |cut -f 2)
Unique_r=$(cat ./${line}/align/STAR/Log.final.out |grep "Uniquely mapped reads %" |cut -f 2)
Multiple=$(cat ./${line}/align/STAR/Log.final.out |grep "Number of reads mapped to multiple loci" |cut -f 2)
Multiple_r=$(cat ./${line}/align/STAR/Log.final.out |grep "% of reads mapped to multiple loci" |cut -f 2)
#rRNA=$(cat ./${line}/qc/RRNA/rRNA.sam.frag_count )
#rRNA_r=$(echo "scale=2;${rRNA}*100/(${Unique}+${Multiple})" |bc )"%"

genes_1=$(cat gene_result/${line}.genes.results |cut -f 3 |tail -n +2 |gawk '$1 >0' |wc -l )
genes_2=$(cat gene_result/${line}.genes.results |cut -f 3 |tail -n +2 |gawk '$1 >2' |wc -l )

#echo ${line},${Reads},${Unique},${Unique_r},${Multiple},${Multiple_r},${rRNA},${rRNA_r},${genes_1},${genes_2} \
#>> summary.csv
echo ${line},${Reads},${Unique},${Unique_r},${Multiple},${Multiple_r},${genes_1},${genes_2} \
>> summary0.csv
done
```

### 8. Expression matrix
```shell
list_sample='/data/user_05/projects/shalek_nature_2013/SRR_Acc_List.txt'
output_dir='/data/user_05/projects/shalek_nature_2013/output'
mat_name='shalek'

##Matix
#gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
#gene_id 	expected_count	TPM	FPKM

# TPM
cd ${output_dir}
idx=$(cat ${list_sample} |head -n 1)
cat gene_result/${idx}.genes.results | cut -f 1 > RNAseq.${mat_name}.tpm.matrix
for line in $(cat ${list_sample})
do
echo ${line} > ${line}.tmp
cat gene_result/${line}.genes.results |cut -f 3 |tail -n +2 >> ${line}.tmp
paste -d "\t" RNAseq.${mat_name}.tpm.matrix ${line}.tmp > RNAseq.${mat_name}.tpm.matrix.tmp
rm RNAseq.${mat_name}.tpm.matrix
mv RNAseq.${mat_name}.tpm.matrix.tmp RNAseq.${mat_name}.tpm.matrix
rm ${line}.tmp
done
```

### 9. Check IS primer
```shell
output_dir='/data/user_05/projects/shalek_nature_2013/output'
r1='/data/user_05/projects/shalek_nature_2013/data/fastq/SRR578551_1.fastq.gz'

cd ${output_dir}
echo '' > read1_tso.txt
zcat $r1 |grep AAGCAGTGGTATCAACGCAGAG >> read1_tso.txt
zcat $r1 |grep CTCTGCGTTGATACCACTGCTT >> read1_tso.txt
zcat $r1 |grep GAGACGCAACTATGGTGACGAA >> read1_tso.txt
zcat $r1 |grep TTCGTCACCATAGTTGCGTCTC >> read1_tso.txt
wc -l read1_tso.txt
```

### 10. BLAST unmapped reads
```shell
#Get unmapped reads
#STAR alignment to genome
output_dir='/data/user_05/projects/shalek_nature_2013/output'
STAR_INDEX='/data/user_05/projects/shalek_nature_2013/data/starindex_mm10'

r1='/data/user_05/projects/shalek_nature_2013/data/fastq/SRR578551_1.fastq.gz'
r2='/data/user_05/projects/shalek_nature_2013/data/fastq/SRR578551_2.fastq.gz'

STAR --genomeDir ${STAR_INDEX} \
--runThreadN 1 --readFilesIn $r1 $r2 \
--readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within KeepPairs \
--twopassMode Basic --limitBAMsortRAM 300000000000 \
--outFileNamePrefix ${output_dir}/analysis_mapping/ \
--outSAMstrandField intronMotif

cd ${output_dir}/analysis_mapping
#filter unmapped reads
samtools view -b -f 4 Aligned.sortedByCoord.out.bam > unmapped.bam
samtools bam2fq unmapped.bam > unmapped.fastq
cat unmapped.fastq |grep '^@.*/1$' -A 3 --no-group-separator > unmapped.r1.fastq
cat unmapped.fastq |grep '^@.*/2$' -A 3 --no-group-separator > unmapped.r2.fastq

#blast unmapped reads
#convert fastq to fasta
cat unmapped.r1.fastq |grep '^@.*/1$' -A 1 --no-group-separator > unmapped.r1.fasta
sed -i -e 's/^\@/\>\@/' unmapped.r1.fasta
cat unmapped.r2.fastq |grep '^@.*/2$' -A 1 --no-group-separator > unmapped.r2.fasta
sed -i -e 's/^\@/\>\@/' unmapped.r2.fasta
head -n 1000 unmapped.r1.fasta > unmapped.r2.1000.fasta
head -n 1000 unmapped.r2.fasta > unmapped.r2.1000.fasta

#BLAST
blastn -db /data/user_05/projects/shalek_nature_2013/data/blast/prok_db_v5/ref_prok_v5 \
-query /data/user_05/projects/shalek_nature_2013/output/analysis_mapping/unmapped.r1.1000.fasta \
-evalue 1e-10 -word_size 28 \
-outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 10 > unmapped.r1.blast.tsv
```

### Supplement: convert chromosome names
For reference from different sources
```shell
sed 's/chr/chromosome /g' /data/user_05/projects/shalek_nature_2013/data/reference/gencode.vM27.annotation.gtf > /data/user_05/projects/shalek_nature_2013/data/reference/new.gtf
```

### Reproduce fig. 2b
```R
library(dplyr)
####prepare data####
df = as.data.frame(read.table('RNAseq.shalek.tpm.txt'))
colnames(df) <- df[1,]
df[1:nrow(df)-1,] <- df[2:nrow(df),]
df <- df[1:nrow(df)-1,]

#log transform TPM
log_df <- mutate_at(.tbl = df,.vars = vars(2:ncol(df)),.fun = as.numeric)
log_df[,2:ncol(df)] <- log(log_df[,2:ncol(df)]+1)

#read 522 genes from supplementary table 3
library(readxl)
sup <- read_excel("Supplementary_Table3.xls",sheet=1)

#convert gene name to ensembl id
library(biomaRt)
mart <- useMart("ensembl","mmusculus_gene_ensembl")
gene_name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="external_gene_name",values = sup[1:522,1],mart = mart) 
log_df$ensembl_gene_id=gsub("\\..*","",df[,1]) #convert to standard ensembl id
names(gene_name) <- c("ensembl_gene_id","Gene")
library(h2o)
gene_name$Gene <- toupper(gene_name$Gene)
sup_merge <- left_join(sup,gene_name,by="Gene")
sup_merge <- as.data.frame(sup_merge)
#write.csv(sup_merge,file='522genes.csv')


####cell count####
countdf <- merge(sup_merge,log_df,all=FALSE)
cellcount <- countdf
#only count single cell samples
#SRR578551,554,557,563,566,569,572,575,578,581,584,587,590,593,596,599

for (N in 0:9) {
  
  for (i in 1:nrow(countdf)){
    cellcount[i,N+8] <- length(countdf[i,8:23][countdf[i,8:23]>=N & countdf[i,8:23]<N+1])
  }
}

colnames(cellcount)[8] <- "0-1"
colnames(cellcount)[9] <- "1-2"
colnames(cellcount)[10] <- "2-3"
colnames(cellcount)[11] <- "3-4"
colnames(cellcount)[12] <- "4-5"
colnames(cellcount)[13] <- "5-6"
colnames(cellcount)[14] <- "6-7"
colnames(cellcount)[15] <- "7-8"
colnames(cellcount)[16] <- "8-9"
colnames(cellcount)[17] <- "9-10"

cellcount <- cellcount[,1:17]
cellcount_o <- cellcount[order(cellcount$`Fano factor`),] #CV low to high, top to bottom
cellcount_o <- as.matrix(cellcount_o[,8:17])


####plot heatmap####
library(pheatmap)
pheatmap(cellcount_o,border_color=NA,scale='none',legend_breaks=c(0,16),legend_labels=c("0","16"),cluster_rows=FALSE,cluster_cols=FALSE,angle_col=0,labels_row="",colorRampPalette(c("black","yellow"))(50))
```
