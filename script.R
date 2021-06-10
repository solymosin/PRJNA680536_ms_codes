
# export PATH=$PATH:/data/tools/FastQC
# 
# for f in *fastq.gz
# do 
#   fastqc -t 38 $f -o qc
# done 
# 
# multiqc qc
# 
# TRIM='/data/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
# 
# for f in *.fastq.gz
# do 
#   o=${f/'.fastq.gz'/'_trimmed.fastq'}
#   java -jar $TRIM SE -phred33 -threads 38 \
#     $f $o \
#     LEADING:3 \
#     TRAILING:3 \
#     SLIDINGWINDOW:4:30 \
#     MINLEN:50
# done
# 
# export PATH=/data/tools/STAR-2.7.3a/bin/Linux_x86_64:$PATH
# export PATH=/data/tools/subread-2.0.0-source/bin:$PATH
# 
# cd /data/dbs/STAR/GRCh38_99
# 
# STAR -- runMode genomeGenerate \
#      -- genomeDir . \
#      -- genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#      -- sjdbGTFfile Homo_sapiens.GRCh38.99.gtf \
#      -- runThreadN 14 
#  
# idx='/data/dbs/STAR/GRCh38_99'
# 
# for f in *.fastq
# do 
#   root=${f/'.fastq'/''}
#   STAR -- genomeDir $idx \
#        -- readFilesIn $f \
#        -- outFileNamePrefix 'GRCh38_99_'$root \
#        -- outFilterMultimapNmax 1 \
#        -- outReadsUnmapped Fastx \
#        -- outSAMtype BAM SortedByCoordinate \
#        -- twopassMode Basic \
#        -- runThreadN 14  
# done    
# 
# featureCounts \
#   -a $idx/Homo_sapiens.GRCh38.99.gtf \
#   -o ms_featureCounts_GRCh38_99.txt \
#   GRCh38_99*bam

library(Biostrings)
library(tidyverse)
library(ShortRead)

a = readFastq('A-CHAT-2_S18.fastq.gz')
at = readFastq('A-CHAT-2_S18_trimmed.fastq')
b = readFastq('B-spiny-2_S19.fastq.gz')
bt = readFastq('B-spiny-2_S19_trimmed.fastq')
c = readFastq('C-CHAT-833-2_S20.fastq.gz')
ct = readFastq('C-CHAT-833-2_S20_trimmed.fastq')
d = readFastq('D-medium-spiny-833-2_S21.fastq.gz')
dt = readFastq('D-medium-spiny-833-2_S21_trimmed.fastq')

read.counts = read_delim('featureCounts_GRCh38_99.txt', col_names=T, skip=1, delim='\t')
colnames(read.counts) = data.frame(col= gsub('_trimmedAligned.sortedByCoord.out.bam', '', gsub('GRCh38_99_', '', colnames(read.counts)))) %>%
separate(col, 'eleje', sep='_') %>% 
pull(eleje)

cl = read.counts %>% select(-c(1:6)) %>% colSums()

raw = c(length(a), length(b), length(c), length(d))
trim = c(length(at), length(bt), length(ct), length(dt))

mapped = c(23361624, 17541978, 1748441, 18394200)
  
feats = c(41.2, 24.5, 20.1, 35.7)

# QoRTs='/data/tools/QoRTs/QoRTs.jar'
# 
# mkdir outputQC_A
# mkdir outputQC_B
# mkdir outputQC_C
# mkdir outputQC_D
# 
# java -Xmx128G -jar $QoRTs QC --maxPhredScore 45 --singleEnded \
#   GRCh38_99_A-CHAT-2_S18_trimmedAligned.sortedByCoord.out.bam \
#   $idx/Homo_sapiens.GRCh38.99.gtf \
#   outputQC_A
#   
# java -Xmx128G -jar $QoRTs QC --maxPhredScore 45 --singleEnded \
#   GRCh38_99_B-spiny-2_S19_trimmedAligned.sortedByCoord.out.bam \
#   $idx/Homo_sapiens.GRCh38.99.gtf \
#   outputQC_B
#     
# java -Xmx128G -jar $QoRTs QC --maxPhredScore 45 --singleEnded \
#   GRCh38_99_C-CHAT-833-2_S20_trimmedAligned.sortedByCoord.out.bam \
#   $idx/Homo_sapiens.GRCh38.99.gtf \
#   outputQC_C
#   
# java -Xmx128G -jar $QoRTs QC --maxPhredScore 45 --singleEnded \
#   GRCh38_99_D-medium-spiny-833-2_S21_trimmedAligned.sortedByCoord.out.bam \
#   $idx/Homo_sapiens.GRCh38.99.gtf \
#   outputQC_D
    

library(QoRTs)

decode = data.frame(
unique.ID=c('A', 'B', 'C', 'D'), 
qc.data.dir=c('outputQC_A', 'outputQC_B', 'outputQC_C', 'outputQC_D')
)

res = read.qc.results.data('./', decoder=decode, calc.DESeq2 = TRUE, calc.edgeR = TRUE)

tab = get.summary.table(res)

tab[grep('ReadPairs', rownames(tab)),]

rtab = tab[grep('ReadPairs', rownames(tab)),]

# ReadPairs_UniqueGene == cl

ptab = rbind(raw, trim, mapped)
colnames(ptab) = c('A', 'B', 'C', 'D')

rt = tibble(counts=gsub('ReadPairs_', '', rownames(rtab)), rtab)
pt = tibble(counts=rownames(ptab), as.data.frame(ptab))

#######################################################################
# library(AnnotationHub)
# 
# ah = AnnotationHub()
# 
# ahDb = query(ah, pattern = c("Homo sapiens", "EnsDb", 99))
# ahEdb = ahDb[['AH78783']]
# 
# makeEnsembldbPackage(
#  ensdb = dbfile(dbconn(ahEdb)), 
#  version = '0.01',
#  maintainer = 'Norbert Solymosi <solymosi.norbert@gmail.com>',
#  author = 'N Solymosi'
# )
# 
# R CMD build EnsDb.Hsapiens.v99
# 
# R CMD INSTALL EnsDb.Hsapiens.v99_0.01.tar.gz
# 

# R

options(java.parameters = "- Xmx204800m")

library(EnsDb.Hsapiens.v99)
library(org.Hs.eg.db)
library(GO.db)
library(PANTHER.db)
library(KEGGREST)
library(edgeR)
library(DESeq2)
library(WriteXLS)
library(xlsx)
library(tidyverse)


org.db = org.Hs.eg.db
edb = EnsDb.Hsapiens.v99
pthOrganisms(PANTHER.db) = 'HUMAN'

# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
# https://www.uniprot.org/help/accession_numbers

uniprot.db_sel = read_tsv('HUMAN_9606_idmapping_selected.tab', col_names=c('UniProtKB_AC', 'UniProtKB_ID', 'GeneID_EntrezGene', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI_taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed')) 

# wget ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR15.0_human

PTHR = read_tsv(
    'PTHR15.0_human', 
    col_names=c('GeneIdentifier', 'ProteinID', 'SFID', 'FamilyName', 'SubfamilyName', 'MolecularFunction', 'BiologicalProcess', 'CellularComponents', 'ProteinClass', 'Pathway'), 
    quote='') %>% 
  separate(GeneIdentifier, c('GeneIdentifier', 'UniProt'), 'UniProtKB=')
  
# wget ftp://ftp.ensembl.org/pub/release-99/tsv/homo_sapiens/Homo_sapiens.GRCh38.99.uniprot.tsv.gz

ens2uniprot = read_tsv('Homo_sapiens.GRCh38.99.uniprot.tsv', col_names=T) %>%
  rename(UniProt=xref, ens=gene_stable_id) %>%
  select(ens, UniProt) %>% 
  unique()

# wget http://geneontology.org/gene-associations/goa_human.gaf.gz
# wget http://geneontology.org/gene-associations/goa_human_isoform.gaf.gz
# 
# gunzip goa_*.gz

go = read_tsv('goa_human.gaf', skip=31, col_names=F, 
  col_types = cols(.default = "c")) %>%
  select(2,5,7,9,10,11) %>%
  mutate(X9=case_when(X9=='C' ~ 'CC', X9=='P' ~ 'BP', X9=='F' ~ 'MF')) %>%
  rename(UniProt=1, GOid=2, evidence=3, ontology=4, GeneName=5, symbol=6)

read.counts = read_delim('featureCounts_GRCh38_99.txt', col_names=T, skip=1, delim='\t')
colnames(read.counts) = data.frame(col= gsub('_trimmedAligned.sortedByCoord.out.bam', '', gsub('GRCh38_99_', '', colnames(read.counts)))) %>%
separate(col, 'eleje', sep='_') %>% 
pull(eleje)

read.counts = as.data.frame(read.counts)
rownames(read.counts) = read.counts$Geneid
colnames(read.counts) = gsub('-', '_', colnames(read.counts))

readcounts = read.counts[,-c(1,2,3,4,5,6)]
sample_info = data.frame(smpl=colnames(readcounts))
rownames(sample_info) = sample_info$smpl

sample_info$grp = 'grpA'
sample_info$grp[2:3] = 'grpB'
sample_info$grp = factor(sample_info$grp)

dds = DESeqDataSetFromMatrix(
    countData = readcounts, 
    colData = sample_info, 
    design = ~ grp 
)

dds = estimateSizeFactors(dds)

counts_raw = readcounts
colnames(counts_raw) = paste0('raw_', colnames(counts_raw))
counts_raw$ens = rownames(counts_raw)
counts_raw = as_tibble(counts_raw)

counts_normalized = as.data.frame(counts(dds, normalized=T))
colnames(counts_normalized) = paste0('norm_', colnames(counts_normalized))
counts_normalized$ens = rownames(counts_normalized)
counts_normalized = as_tibble(counts_normalized)

m = as.matrix(readcounts)
counts_cpm = as.data.frame(cpm(m))
colnames(counts_cpm) = paste0('cpm_', colnames(counts_cpm))
counts_cpm$ens = rownames(counts_cpm)
counts_cpm = as_tibble(counts_cpm)

# ids = rownames(readcounts)
# n = 1
# tmp = tibble(.rows=0, ens='', UniProt='')
# for(ens in ids){
#   tmp = rbind(tmp, 
#     tibble(ens, UniProt= uniprot.db_sel %>% 
#     filter(str_detect(Ensembl, ens)) %>%
#     pull(UniProtKB_AC)  
#     )
#   )
#   n=n+1
#   print(n)
# }
# save(tmp, file='ens_tmp.RData')

load('ens_tmp.RData')

ens2uniprot = rbind(ens2uniprot, tmp) %>% 
  unique()

ens_uniprot_pthr = inner_join(ens2uniprot, PTHR)

tib = left_join(tibble(ens=rownames(readcounts)), ens_uniprot_pthr) %>% 
  select(ens, UniProt)
  
tib = inner_join(tib, counts_raw)
tib = inner_join(tib, counts_normalized)
tib = inner_join(tib, counts_cpm)

dds.dif = DESeq(dds)

res = results(dds.dif, contrast=c('grp', 'grpA', 'grpB'))
fix = res[tib$ens,]
tib$log2FC = fix$log2FoldChange

annot = left_join(
  tib,
  PTHR %>% 
    select(UniProt, FamilyName, SubfamilyName, ProteinClass, BiologicalProcess, CellularComponents, MolecularFunction)
)

evs = sort(unique(go$evidence))
onts = sort(unique(go$ontology))

# ev = evs[1]
# ont = onts[1]
for(ont in onts){
  for(ev in evs){
    tmp = go %>% 
      filter(ontology==ont, evidence==ev) %>% 
      select(UniProt, GeneName)
    if(dim(tmp)[1]>0){ 
      lst = split(tmp$GeneName, tmp$UniProt) %>% 
        lapply(unique) %>% 
        lapply(sort) %>%  
        lapply(paste, collapse='\n')   
      annot = left_join(
        annot, 
        tibble(UniProt=names(lst), tmp=as.character(lst)) %>%
          rename_at(vars(tmp), ~ paste(ont, ev, sep='_'))
        )
    }
  }
}

lst = keggList('pathway', 'hsa') 
paths = tibble(
  PATH=gsub('path:hsa', '', names(lst)), 
  pathway=gsub(' - Homo sapiens \\(human\\)', '', as.character(lst))
) 

i = 1
query = keggGet(paste0('hsa', paths$PATH[i]))
kegg = as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
  rename(GeneID=1, descr=2) %>%
  mutate(PATH=paths$PATH[i])
  
for(i in 2:dim(paths)[1]){
  query = keggGet(paste0('hsa', paths$PATH[i]))
  if(!is.null(query[[1]]$GENE)){
    kegg = rbind(kegg, 
      as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
      rename(GeneID=1, descr=2) %>%
      mutate(PATH=paths$PATH[i])
    )
  }
#   print(length(query))
}

ens2entrez = read_tsv('Homo_sapiens.GRCh38.99.entrez.tsv', 
  col_types = cols(.default = "c"))

tmp = inner_join(
  inner_join(kegg,  
    inner_join(
      tibble(gene_stable_id=rownames(readcounts)),
      ens2entrez
    ) %>% 
    select(gene_stable_id, xref) %>% 
    unique() %>% 
    rename(ens=1, GeneID=2)
  ) %>% 
  select(ens, PATH),
  paths
)

lst = split(tmp$pathway, tmp$ens) %>% 
  lapply(unique) %>% 
  lapply(sort) %>%  
  lapply(paste, collapse='\n') 
kegg_res = tibble(ens=names(lst), KEGG=as.character(lst))  

annot = left_join(annot, kegg_res) %>% 
  rename(Ensembl=1)

s = (dim(tib)[2]+1):dim(annot)[2]
wb = createWorkbook() 
cs = CellStyle(wb, alignment = Alignment(wrapText = TRUE))
colsty = rep(list(cs), length(s))
names(colsty) = s

sheet = createSheet(wb, sheetName='Chat_GAD_with_annotation')
addDataFrame(annot, sheet, colStyle=colsty)

saveWorkbook(wb, 'results.xlsx')

