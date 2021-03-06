## Homework for 4.1 Sequence Motif
#### 1. Explain the concatenation step
The step “concatenate sequences of the same 3'UTR” picks out all sequence reads of the same 3'UTR and joins the sequences together.
Take transcript ENST00000612051.1 as an example:
In the fasta file before concatenation, it has several separate reads:
```
>ENST00000612051.1(+)
GCACCTCTGGAAGAGCCAACTGTGTGAGATGGTGCAGCCCAGTGGTGGCCCGGCAGCAGATCAGGACGTACTGGGCGAAGAGTCTCCTCTGGGGAAGCCAGCCATGCTGCACCTGCCTTCAGAACAGGGCGCTCCTGAGACCCTCCAGCGCTGCCTGGAGGAGAATCAAGAGCTCCGAG
>ENST00000612051.1(+)
ATGCCATCCGGCAGAGCAACCAGATTCTGCGGGAGCGCTGCGAGGAGCTTCTGCATTTCCAAGCCAGCCAGAGGGAGGAGAAGGAGTTCCTCATGTGCAAGTTCCAGGAGGCCAGGAAACTGGTGGAGAGACTCGGCCTGGAGAAGCTCGATCTGAAGAGGCAGAAGGAGCAGGCTCTGCGGGAGGTGGAGCACCTGAAGAGATGCCAGCAG
--
>ENST00000612051.1(+)
CAGATGGCTGAGGACAAGGCCTCTGTGAAAGCCCAGGTGACGTCCTTGCTCGGGGAGCTGCAGGAGAGCCAGAGTCGCTTGGAGGCTGCCACTAAGGAATGCCAGGCTCTGGAGGGTCG
>ENST00000612051.1(+)
GGCCCGGGCGGCCAGCGAGCAGGCGCGGCAGCTGGAGAGTGAGCGCGAGGCGCTGCAGCAGCAGCACAGCGTGCAGGTGGACCAGCTGCGCATGCAGGGCCAGAGCGTGGAGGCCGCGCTCCGCATGGAGCGCCAGGCCGCCTCGGAGGAGAA
>ENST00000612051.1(+)
GAGGAAGCTGGCCCAGTTGCAGGTGGCCTATCACCAGCTCTTCCAAGAATACGACAACCACATCAAGAGCAGCGTGGTGGGCAGTGAGCGGAAGCGA
>ENST00000612051.1(+)
GGAATGCAGCTGGAAGATCTCAAACAGCAGCTCCAGCAGGCCGAGGAGGCCCTGGTGGCCAAACAGGAGGTGATCGATAAGCTGAAGGAGGAGGCCGAGCAGCACAAGATTGTGATGGAGACCGTTCCGGTGCTGAAGGCCCAG
>ENST00000612051.1(+)
GCGGATATCTACAAGGCGGACTTCCAGGCTGAGAGGCAGGCCCGGGAGAAGCTGGCCGAGAAGAAGGAGCTCCTGCAGGAGCAGCTGGAGCAGCTGCAGAGGGAGTACAGCAAACTGAAGGCCAGCTGTCAGGAGTCGGCCAG
>ENST00000612051.1(+)
GATCGAGGACATGAGGAAGCGGCATGTCGAGGTCTCCCAGGCCCCCTTGCCCCCCGCCCCTG
>ENST00000612051.1(+)
CCTACCTCTCCTCTCCCCTGGCCCTGCCCAGCCAGAGGAGGAGCCCCCCCGAGGAGCCACCTGACTTCTGCTGTCCCAAGTGCCAGTATCAGGCCCCTGATATGGACACCCTGCAGATACATGTCATGGAGTGCATTGAGTAGGGCCGGCCAGTGCAAGGCCACTGCCTGCCGAGGACGTGCCCGGGACCGTGCAGTCTGCGCTTTCCTCTCCCGCCTGCCTAGCCCAGGATGAAGGGCTGGGTGGCCACAACTGGGATGCCACCTGGAGCCCCACCCAGGAGCTGGCCGCGGCACCTTACGCTTCAGCTGTTGATCCGCTGGTCCCCTCTTTTGGGGTAGATGCGGCCCCGATCAGGCCTGACTCGCTGCTCTTTTTGTTCCCTTCTGTCTGCTCGAACCACTTGCCTCGGGCTAATCCCTCCCTCTTCCTCCACCCGGCACTGGGGAAGTCAAGAATGGGGCCTGGGGCTCTCAGGGAGAACTGCTTCCCCTGGCAGAGCTGGGTGGCCGCTCTTCCTCCCACCGGACACCGACCCGCCCGCCGCTGTGCCCTGGGAGTGCTGCCCTCTTACCATGCACACGGGTGCTCTCCTTTTGGGCTGCATGCTATTCCATTTTGCAGCCAGACCGATGTGTATTTAACCAGTCACTATTGATGGACATTTGGGTTGTTTCCCATCTTTTTGTTACCATAAATAATGGCATAGTAAAAA
```
After concatenation, it has only one concatenated read:
```
>ENST00000612051.1
GCACCTCTGGAAGAGCCAACTGTGTGAGATGGTGCAGCCCAGTGGTGGCCCGGCAGCAGATCAGGACGTACTGGGCGAAGAGTCTCCTCTGGGGAAGCCAGCCATGCTGCACCTGCCTTCAGAACAGGGCGCTCCTGAGACCCTCCAGCGCTGCCTGGAGGAGAATCAAGAGCTCCGAGATGCCATCCGGCAGAGCAACCAGATTCTGCGGGAGCGCTGCGAGGAGCTTCTGCATTTCCAAGCCAGCCAGAGGGAGGAGAAGGAGTTCCTCATGTGCAAGTTCCAGGAGGCCAGGAAACTGGTGGAGAGACTCGGCCTGGAGAAGCTCGATCTGAAGAGGCAGAAGGAGCAGGCTCTGCGGGAGGTGGAGCACCTGAAGAGATGCCAGCAGCAGATGGCTGAGGACAAGGCCTCTGTGAAAGCCCAGGTGACGTCCTTGCTCGGGGAGCTGCAGGAGAGCCAGAGTCGCTTGGAGGCTGCCACTAAGGAATGCCAGGCTCTGGAGGGTCGGGCCCGGGCGGCCAGCGAGCAGGCGCGGCAGCTGGAGAGTGAGCGCGAGGCGCTGCAGCAGCAGCACAGCGTGCAGGTGGACCAGCTGCGCATGCAGGGCCAGAGCGTGGAGGCCGCGCTCCGCATGGAGCGCCAGGCCGCCTCGGAGGAGAAGAGGAAGCTGGCCCAGTTGCAGGTGGCCTATCACCAGCTCTTCCAAGAATACGACAACCACATCAAGAGCAGCGTGGTGGGCAGTGAGCGGAAGCGAGGAATGCAGCTGGAAGATCTCAAACAGCAGCTCCAGCAGGCCGAGGAGGCCCTGGTGGCCAAACAGGAGGTGATCGATAAGCTGAAGGAGGAGGCCGAGCAGCACAAGATTGTGATGGAGACCGTTCCGGTGCTGAAGGCCCAGGCGGATATCTACAAGGCGGACTTCCAGGCTGAGAGGCAGGCCCGGGAGAAGCTGGCCGAGAAGAAGGAGCTCCTGCAGGAGCAGCTGGAGCAGCTGCAGAGGGAGTACAGCAAACTGAAGGCCAGCTGTCAGGAGTCGGCCAGGATCGAGGACATGAGGAAGCGGCATGTCGAGGTCTCCCAGGCCCCCTTGCCCCCCGCCCCTGCCTACCTCTCCTCTCCCCTGGCCCTGCCCAGCCAGAGGAGGAGCCCCCCCGAGGAGCCACCTGACTTCTGCTGTCCCAAGTGCCAGTATCAGGCCCCTGATATGGACACCCTGCAGATACATGTCATGGAGTGCATTGAGTAGGGCCGGCCAGTGCAAGGCCACTGCCTGCCGAGGACGTGCCCGGGACCGTGCAGTCTGCGCTTTCCTCTCCCGCCTGCCTAGCCCAGGATGAAGGGCTGGGTGGCCACAACTGGGATGCCACCTGGAGCCCCACCCAGGAGCTGGCCGCGGCACCTTACGCTTCAGCTGTTGATCCGCTGGTCCCCTCTTTTGGGGTAGATGCGGCCCCGATCAGGCCTGACTCGCTGCTCTTTTTGTTCCCTTCTGTCTGCTCGAACCACTTGCCTCGGGCTAATCCCTCCCTCTTCCTCCACCCGGCACTGGGGAAGTCAAGAATGGGGCCTGGGGCTCTCAGGGAGAACTGCTTCCCCTGGCAGAGCTGGGTGGCCGCTCTTCCTCCCACCGGACACCGACCCGCCCGCCGCTGTGCCCTGGGAGTGCTGCCCTCTTACCATGCACACGGGTGCTCTCCTTTTGGGCTGCATGCTATTCCATTTTGCAGCCAGACCGATGTGTATTTAACCAGTCACTATTGATGGACATTTGGGTTGTTTCCCATCTTTTTGTTACCATAAATAATGGCATAGTAAAAA
```

#### 2. Do concatenation with my own script
R script
```R
########### concatenate with multiple transcripts###########
before.fa <- read.table("before_test.txt") #fasta file before concatenation
#after.fa <- read.table("after_test.txt")
out <- file("output.txt",open="wt") #output path

#define concatenate function
concatenate.func <- function(df){
  writeLines(df$transcript.id.trim[1],con=out)
  df.cont <- paste(df$sequence,collapse="")
  writeLines(df.cont,con=out)
}

#prepare dataframe from fasta
library(stringr)
transcript.id <- before.fa[grep("^>ENST",before.fa[,1]),1]
transcript.id.trim <- substr(transcript.id,1,18)
sequence <- before.fa[grep("G|A",before.fa[,1]),1]
before.df <- as.data.frame(cbind(transcript.id,transcript.id.trim,sequence))

#find sequence of the same transcript
group <- split(before.df,before.df$transcript.id.trim)

for (i in 1:length(group)){
  concatenate.func(group[[i]])
}

close(out)
```
Use the following test script with two transcripts as input:
```
>ENST00000000442.10(+)
GGCAAGGGGTGGGACTGGTGGGGGTTCTGGCAGGACCTGCCTAGCATGGGGTCAGCCCCAAGGGCTGGGGCGGAGCTGGGGTCTGGGCAGTGCCACAGCCTGCTGGCAGGGCCAGGGCAATGCCATCAGCCCCTGGGAACAGGCCCCACGCCCTCTCCTCCCCCTCCTAGGGGGTGTCAGAAGCTGGGAACGTGTGTCCAGGCTCTGGGCACAGTGCTGCCCCTTGCAAGCCATAACGTGCCCCCAGAGTGTAGGGGGCCTTGCGGAAGCCATAGGGGGCTGCACGGGATGCGTGGGAGGCAGAAACCTATCTCAGGGAGGGAAGGGGATGGAGGCCAGAGT
>ENST00000000442.10(+)
CTCCCAGTGGGTGATGCTTTTGCTGCTGCTTAATCCTACCCCCTCTTCAAAGCAGAGTGGGACTTGGAGAGCAAAGGCCCATGCCCCCTTCGCTCCTCCTCTCATCATTTGCATTGGGCATTAGTGTCCCCCCTTGAAGCAATAACTCCAAGCAGACTCCAGCCCCTGGACCCCTGGGGTGGCCAGGGCTTCCCCATCAGCTCCCAACGAGCCTCCTCAGGGGGTAGGAGAGCACTGCCTCTATGCCCTGCAGAGCAATAACACTATATTTATTTTTGGGTTTGGCCAGGGAGGCGCAGGGACATGGGGCAAGCCAGGGCCCAGAGCCCTTGGCTGTACAGAGACTCTATTTTAATGTATATTTGCTGCAAAGAGAAACCGCTTTTGGTTTTAAACCTTTAATGAGAAAAAAATATATAATACCGAGCTC
>ENST00000612051.1(+)
GCACCTCTGGAAGAGCCAACTGTGTGAGATGGTGCAGCCCAGTGGTGGCCCGGCAGCAGATCAGGACGTACTGGGCGAAGAGTCTCCTCTGGGGAAGCCAGCCATGCTGCACCTGCCTTCAGAACAGGGCGCTCCTGAGACCCTCCAGCGCTGCCTGGAGGAGAATCAAGAGCTCCGAG
>ENST00000612051.1(+)
ATGCCATCCGGCAGAGCAACCAGATTCTGCGGGAGCGCTGCGAGGAGCTTCTGCATTTCCAAGCCAGCCAGAGGGAGGAGAAGGAGTTCCTCATGTGCAAGTTCCAGGAGGCCAGGAAACTGGTGGAGAGACTCGGCCTGGAGAAGCTCGATCTGAAGAGGCAGAAGGAGCAGGCTCTGCGGGAGGTGGAGCACCTGAAGAGATGCCAGCAG
--
>ENST00000612051.1(+)
CAGATGGCTGAGGACAAGGCCTCTGTGAAAGCCCAGGTGACGTCCTTGCTCGGGGAGCTGCAGGAGAGCCAGAGTCGCTTGGAGGCTGCCACTAAGGAATGCCAGGCTCTGGAGGGTCG
>ENST00000612051.1(+)
GGCCCGGGCGGCCAGCGAGCAGGCGCGGCAGCTGGAGAGTGAGCGCGAGGCGCTGCAGCAGCAGCACAGCGTGCAGGTGGACCAGCTGCGCATGCAGGGCCAGAGCGTGGAGGCCGCGCTCCGCATGGAGCGCCAGGCCGCCTCGGAGGAGAA
>ENST00000612051.1(+)
GAGGAAGCTGGCCCAGTTGCAGGTGGCCTATCACCAGCTCTTCCAAGAATACGACAACCACATCAAGAGCAGCGTGGTGGGCAGTGAGCGGAAGCGA
>ENST00000612051.1(+)
GGAATGCAGCTGGAAGATCTCAAACAGCAGCTCCAGCAGGCCGAGGAGGCCCTGGTGGCCAAACAGGAGGTGATCGATAAGCTGAAGGAGGAGGCCGAGCAGCACAAGATTGTGATGGAGACCGTTCCGGTGCTGAAGGCCCAG
>ENST00000612051.1(+)
GCGGATATCTACAAGGCGGACTTCCAGGCTGAGAGGCAGGCCCGGGAGAAGCTGGCCGAGAAGAAGGAGCTCCTGCAGGAGCAGCTGGAGCAGCTGCAGAGGGAGTACAGCAAACTGAAGGCCAGCTGTCAGGAGTCGGCCAG
>ENST00000612051.1(+)
GATCGAGGACATGAGGAAGCGGCATGTCGAGGTCTCCCAGGCCCCCTTGCCCCCCGCCCCTG
>ENST00000612051.1(+)
CCTACCTCTCCTCTCCCCTGGCCCTGCCCAGCCAGAGGAGGAGCCCCCCCGAGGAGCCACCTGACTTCTGCTGTCCCAAGTGCCAGTATCAGGCCCCTGATATGGACACCCTGCAGATACATGTCATGGAGTGCATTGAGTAGGGCCGGCCAGTGCAAGGCCACTGCCTGCCGAGGACGTGCCCGGGACCGTGCAGTCTGCGCTTTCCTCTCCCGCCTGCCTAGCCCAGGATGAAGGGCTGGGTGGCCACAACTGGGATGCCACCTGGAGCCCCACCCAGGAGCTGGCCGCGGCACCTTACGCTTCAGCTGTTGATCCGCTGGTCCCCTCTTTTGGGGTAGATGCGGCCCCGATCAGGCCTGACTCGCTGCTCTTTTTGTTCCCTTCTGTCTGCTCGAACCACTTGCCTCGGGCTAATCCCTCCCTCTTCCTCCACCCGGCACTGGGGAAGTCAAGAATGGGGCCTGGGGCTCTCAGGGAGAACTGCTTCCCCTGGCAGAGCTGGGTGGCCGCTCTTCCTCCCACCGGACACCGACCCGCCCGCCGCTGTGCCCTGGGAGTGCTGCCCTCTTACCATGCACACGGGTGCTCTCCTTTTGGGCTGCATGCTATTCCATTTTGCAGCCAGACCGATGTGTATTTAACCAGTCACTATTGATGGACATTTGGGTTGTTTCCCATCTTTTTGTTACCATAAATAATGGCATAGTAAAAA
```
The output concatenated file is
```
>ENST00000000442.1
GGCAAGGGGTGGGACTGGTGGGGGTTCTGGCAGGACCTGCCTAGCATGGGGTCAGCCCCAAGGGCTGGGGCGGAGCTGGGGTCTGGGCAGTGCCACAGCCTGCTGGCAGGGCCAGGGCAATGCCATCAGCCCCTGGGAACAGGCCCCACGCCCTCTCCTCCCCCTCCTAGGGGGTGTCAGAAGCTGGGAACGTGTGTCCAGGCTCTGGGCACAGTGCTGCCCCTTGCAAGCCATAACGTGCCCCCAGAGTGTAGGGGGCCTTGCGGAAGCCATAGGGGGCTGCACGGGATGCGTGGGAGGCAGAAACCTATCTCAGGGAGGGAAGGGGATGGAGGCCAGAGTCTCCCAGTGGGTGATGCTTTTGCTGCTGCTTAATCCTACCCCCTCTTCAAAGCAGAGTGGGACTTGGAGAGCAAAGGCCCATGCCCCCTTCGCTCCTCCTCTCATCATTTGCATTGGGCATTAGTGTCCCCCCTTGAAGCAATAACTCCAAGCAGACTCCAGCCCCTGGACCCCTGGGGTGGCCAGGGCTTCCCCATCAGCTCCCAACGAGCCTCCTCAGGGGGTAGGAGAGCACTGCCTCTATGCCCTGCAGAGCAATAACACTATATTTATTTTTGGGTTTGGCCAGGGAGGCGCAGGGACATGGGGCAAGCCAGGGCCCAGAGCCCTTGGCTGTACAGAGACTCTATTTTAATGTATATTTGCTGCAAAGAGAAACCGCTTTTGGTTTTAAACCTTTAATGAGAAAAAAATATATAATACCGAGCTC
>ENST00000612051.1
GCACCTCTGGAAGAGCCAACTGTGTGAGATGGTGCAGCCCAGTGGTGGCCCGGCAGCAGATCAGGACGTACTGGGCGAAGAGTCTCCTCTGGGGAAGCCAGCCATGCTGCACCTGCCTTCAGAACAGGGCGCTCCTGAGACCCTCCAGCGCTGCCTGGAGGAGAATCAAGAGCTCCGAGATGCCATCCGGCAGAGCAACCAGATTCTGCGGGAGCGCTGCGAGGAGCTTCTGCATTTCCAAGCCAGCCAGAGGGAGGAGAAGGAGTTCCTCATGTGCAAGTTCCAGGAGGCCAGGAAACTGGTGGAGAGACTCGGCCTGGAGAAGCTCGATCTGAAGAGGCAGAAGGAGCAGGCTCTGCGGGAGGTGGAGCACCTGAAGAGATGCCAGCAGCAGATGGCTGAGGACAAGGCCTCTGTGAAAGCCCAGGTGACGTCCTTGCTCGGGGAGCTGCAGGAGAGCCAGAGTCGCTTGGAGGCTGCCACTAAGGAATGCCAGGCTCTGGAGGGTCGGGCCCGGGCGGCCAGCGAGCAGGCGCGGCAGCTGGAGAGTGAGCGCGAGGCGCTGCAGCAGCAGCACAGCGTGCAGGTGGACCAGCTGCGCATGCAGGGCCAGAGCGTGGAGGCCGCGCTCCGCATGGAGCGCCAGGCCGCCTCGGAGGAGAAGAGGAAGCTGGCCCAGTTGCAGGTGGCCTATCACCAGCTCTTCCAAGAATACGACAACCACATCAAGAGCAGCGTGGTGGGCAGTGAGCGGAAGCGAGGAATGCAGCTGGAAGATCTCAAACAGCAGCTCCAGCAGGCCGAGGAGGCCCTGGTGGCCAAACAGGAGGTGATCGATAAGCTGAAGGAGGAGGCCGAGCAGCACAAGATTGTGATGGAGACCGTTCCGGTGCTGAAGGCCCAGGCGGATATCTACAAGGCGGACTTCCAGGCTGAGAGGCAGGCCCGGGAGAAGCTGGCCGAGAAGAAGGAGCTCCTGCAGGAGCAGCTGGAGCAGCTGCAGAGGGAGTACAGCAAACTGAAGGCCAGCTGTCAGGAGTCGGCCAGGATCGAGGACATGAGGAAGCGGCATGTCGAGGTCTCCCAGGCCCCCTTGCCCCCCGCCCCTGCCTACCTCTCCTCTCCCCTGGCCCTGCCCAGCCAGAGGAGGAGCCCCCCCGAGGAGCCACCTGACTTCTGCTGTCCCAAGTGCCAGTATCAGGCCCCTGATATGGACACCCTGCAGATACATGTCATGGAGTGCATTGAGTAGGGCCGGCCAGTGCAAGGCCACTGCCTGCCGAGGACGTGCCCGGGACCGTGCAGTCTGCGCTTTCCTCTCCCGCCTGCCTAGCCCAGGATGAAGGGCTGGGTGGCCACAACTGGGATGCCACCTGGAGCCCCACCCAGGAGCTGGCCGCGGCACCTTACGCTTCAGCTGTTGATCCGCTGGTCCCCTCTTTTGGGGTAGATGCGGCCCCGATCAGGCCTGACTCGCTGCTCTTTTTGTTCCCTTCTGTCTGCTCGAACCACTTGCCTCGGGCTAATCCCTCCCTCTTCCTCCACCCGGCACTGGGGAAGTCAAGAATGGGGCCTGGGGCTCTCAGGGAGAACTGCTTCCCCTGGCAGAGCTGGGTGGCCGCTCTTCCTCCCACCGGACACCGACCCGCCCGCCGCTGTGCCCTGGGAGTGCTGCCCTCTTACCATGCACACGGGTGCTCTCCTTTTGGGCTGCATGCTATTCCATTTTGCAGCCAGACCGATGTGTATTTAACCAGTCACTATTGATGGACATTTGGGTTGTTTCCCATCTTTTTGTTACCATAAATAATGGCATAGTAAAAA
```
