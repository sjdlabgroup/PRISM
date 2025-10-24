PRISM: PRecise Identification of Species of the Microbiome
================
Bassel Ghaddar
2025-10-15

## Introduction

<div style="float: right; text-align: center;">

<img src="Figure 1.png" width="250px" style="float: right; margin-left: 40px;" />

<div style="text-align: center;">

<strong>Figure 1:</strong> Overview of PRISM.

</div>

</div>

PRISM identifies truly present microbial taxa in genomic sequencing data
and eliminates falsely positive artifacts and contaminants. It
identifies uniquely identifiable taxa using full-length mapping of a
representative subsample of sequencing reads and taxa first identified
with fast, k-mer-based taxonomic classification. It then employs a
machine learning model to predict tissue-present microbes
vs. contaminants based on multiple features engineered from read mapping
and gene expression statistics (Fig. 1). PRISM is compatible with
RNA-seq, WGS, 16S-seq, and scRNA-seq.

Please see the reference below for more information.

Please contact Bassel Ghaddar (<bassel.ghaddar@gmail.com>) for any
questions.

## Setup

PRISM requires the following dependencies:

1.  Kraken2: <https://github.com/DerrickWood/kraken2>
2.  Minimap2:<https://github.com/lh3/minimap2>
3.  STAR: <https://github.com/alexdobin/STAR>
4.  BLAST:
    <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
5.  Seqkit: <https://bioinf.shenwei.me/seqkit/>
6.  R packages: optparse, ShortRead, tidyverse, furrr, data.table, vegan
7.  Processed genbank files:
    <https://drive.google.com/file/d/1nkdQ3GyG_FfQQI7nRWpzcQDJOTvULuGV/view?usp=sharing>

After downloading this PRISM package and genbank folder, make sure to
unzip the genbank folder and place it in the PRISM package directory.

There is a required BLAST database mapping file that must be created
once:

One-time setup: create accession→taxid map for PRISM

1.  Extract taxid-accession pairs from your BLAST database (edit path as
    needed)
    `blastdbcmd -db /path/to/blast/db/eg/core_nt \   -entry all -outfmt "%T %a" > /path/to/PRISM/taxid_accession_map.tsv`

2.  Sort by taxid to create final map
    `sort -k1,1 /path/to/PRISM/taxid_accession_map.tsv > /path/to/PRISM/sorted_accession_map.txt`

(Re-run only if BLAST database is updated)

PRISM is invoked with a single command-line R function `PRISM.R`:

## Running PRISM

### Input

`PRISM.R`

#### Required inputs / flags

- `--sample` sample name (e.g., ABC for ABC_1.fastq)
- `--data_path` path to FASTQ/FASTA files
- `--kraken_path` path to Kraken2 executable
- `--kraken_db_path` path to Kraken2 database
- `--seqkit_path` path to SeqKit executable
- `--minimap2_path` path to Minimap2 executable
- `--minimap2_index` path to Minimap2 host .mmi index
- `--star_path` path to STAR executable
- `--star_genome_dir` path to STAR host genome index (for –genomeDir)
- `--blast_path` path to BLAST+ executables directory
- `--blast_db_path` path prefix to BLAST database (no extension)
- `--prism_path` path to PRISM repository (root folder)

#### Common / optional flags (defaults shown)

- `--model_org_taxids` path to text file of model-organism taxids
  (default: prism_path/NA → bundled file)
- `--paired` TRUE\|FALSE (default TRUE)
- `--fq1_end` suffix for read1 (default “\_1.fastq”)
- `--fq2_end` suffix for read2 (default “\_2.fastq”)
- `--barcode_only` TRUE\|FALSE; if TRUE, skip BLAST of read1 (default
  FALSE)
- `--max_sample` maximum reads per taxon in subsample (default 100)
- `--min_read_per` minimum reads per X to analyze (default 1e4)
- `--min_uniq_frac` minimum unique k-mer ratio (default 5)
- `--threads` number of threads for external tools (default 1)
- `--min_qcovs` minimum BLAST query coverage percent (default 80)
- `--out_path` output directory (default: {data_path}/{sample}\_prism/)
- `--use_custom_db` TRUE\|FALSE; build temporary subsetted BLAST DB
  (default TRUE)
- `--custom_db_path` directory for temporary subsetted BLAST DB
  (default: {out_path}/data/customdb)
- `--fasta_size_threshold_custom_db` minimum FASTA size (MB) to trigger
  custom DB (default 5)

## Output

PRISM creates a directory for each sample containing the main subfolder
`data/` with all intermediate files and logs.  
The key final outputs are:

- **`X-results.csv`** — Final per-read table after all filtering and
  scoring  
- **`X-counts.csv`** — Per-species summary of detected taxa  
- **`X_1.fa` / `X_2.fa`** — Final PRISM FASTA(s) of retained microbial
  reads (headers annotated with taxid and accession)

#### `X-results.csv`

This is the final BLAST alignment table after: 1. Identifying uniquely
mappable microbial species, 2. Removing human, model organism, and
vector sequences, 3. Resolving multi-mapping reads, and 4. Appending
GenBank annotations and PRISM classification scores.

Each row corresponds to a single sequencing read.  
Columns:

| Column | Description |
|:---|:---|
| **id** | Sequence identifier |
| **read** | Read number (`1` or `2`) for paired-end samples |
| **staxids** | NCBI taxon ID assigned |
| **tax_name** | Scientific name of the taxon |
| **rank** | Highest resolved phylogenetic rank (`k,p,c,o,f,g,s`) |
| **sacc** | BLAST subject accession |
| **pos** | Subject start position of alignment |
| **qcovs** | BLAST query coverage percent (0–100) |
| **pident** | Percent sequence identity |
| **bitscore** | BLAST bit score |
| **gene** | GenBank gene annotation (if available) |
| **product** | GenBank product annotation (if available) |
| **pred** | PRISM score (0 = likely contaminant, 1 = likely truly present) |

#### `X-counts.csv`

Aggregated species-level summary of reads and PRISM scores:

| Column       | Description                |
|:-------------|:---------------------------|
| **tax_name** | Scientific name            |
| **staxids**  | NCBI taxon ID              |
| **n**        | Number of reads assigned   |
| **pred**     | Mean PRISM score per taxon |

#### `X_1.fa` / `X_2.fa`

Final FASTA files containing only PRISM-retained reads.  
Each header is annotated as: read_id \| PRISM \| staxids:{taxid}
sacc:{accession}

## Example

#### RNA-seq data from pancreatic cancer, using f942e6d6-f697-4141-8f1c-58933ca81751_1.fastq and f942e6d6-f697-4141-8f1c-58933ca81751_2.fastq from the CPTAC project.

``` bash
Rscript \
/path/to/PRISM/PRISM.R \
--sample f942e6d6-f697-4141-8f1c-58933ca81751 \
--data_path /path/to/data/ \
--kraken_path /path/to/kraken2-master/kraken2 \
--kraken_db_path /path/to/kraken_db \
--seqkit_path /path/to/seqkit \
--minimap2_path /path/to/minimap2 \
--minimap2_index /path/to/host_hg38.mmi \
--star_path /path/to/STAR \
--star_genome_dir /path/to/star_index \
--model_org_taxids /path/to/PRISM/model_org_taxids.txt \
--blast_path /path/to/blast/2.16/bin/ \
--blast_db_path /path/to/BLAST/db \
--prism_path /path/to/PRISM/ \
--fq1_end _1.fastq \
--fq2_end _2.fastq \
--paired T
```

The files `0ac20066-3954-4a36-8ab3-c62a0c32d988-results.csv` and
`0ac20066-3954-4a36-8ab3-c62a0c32d988-counts.csv` contain the main
results.

Species counts and PRISM scores:

``` r
library(tidyverse)
res = read.csv('./test data/0ac20066-3954-4a36-8ab3-c62a0c32d988-counts.csv') 
tibble(res)
```

    ## # A tibble: 13 × 4
    ##    tax_name                staxids     n  pred
    ##    <chr>                     <int> <int> <dbl>
    ##  1 Enterococcus faecium       1352  5521 0.973
    ##  2 Escherichia coli            562  1531 0.105
    ##  3 Pseudomonas aeruginosa      287   277 0.009
    ##  4 Staphylococcus aureus      1280   262 0.005
    ##  5 Bacteroides fragilis        817   228 0.617
    ##  6 Cupriavidus taiwanensis  164546   199 0.002
    ##  7 Paucibacter sediminis   3019553   175 0    
    ##  8 Pseudomonas tolaasii      29442   161 0.002
    ##  9 Cutibacterium acnes        1747   120 0.006
    ## 10 Pseudomonas yamanorum    515393    99 0.004
    ## 11 Parvimonas micra          33033    43 0.006
    ## 12 Prevotella nigrescens     28133    43 0.002
    ## 13 Malassezia restricta      76775    24 0.001

Examining read-level microbial gene/product data:

``` r
library(tidyverse)
res = read.csv('./test data/0ac20066-3954-4a36-8ab3-c62a0c32d988-results.csv') 
tibble(res)
```

    ## # A tibble: 8,683 × 13
    ##    id        read staxids tax_name rank  sacc  gene  product    pos qcovs pident
    ##    <chr>    <int>   <int> <chr>    <chr> <chr> <chr> <chr>    <int> <int>  <dbl>
    ##  1 K00270:…     2     562 Escheri… S     S429… ""    ""         105   100  100  
    ##  2 K00270:…     2    1352 Enteroc… S     MK33… ""    "16S r…    624   100  100  
    ##  3 K00270:…     1    1352 Enteroc… S     LC56… ""    ""         124    99   90.8
    ##  4 K00270:…     1    1352 Enteroc… S     LC56… ""    ""         124   100   92.2
    ##  5 K00270:…     1    1352 Enteroc… S     CP13… ""    ""         263    99  100  
    ##  6 K00270:…     1 3019553 Pauciba… S     CP11… ""    "23S r… 899421   100  100  
    ##  7 K00270:…     2     562 Escheri… S     AP02… ""    "23S r… 457497   100  100  
    ##  8 K00270:…     1 3019553 Pauciba… S     CP11… ""    "23S r… 899421   100  100  
    ##  9 K00270:…     2     562 Escheri… S     CP09… ""    "23S r… 468371   100  100  
    ## 10 K00270:…     1    1352 Enteroc… S     LC56… ""    ""         124   100   89.6
    ## # ℹ 8,673 more rows
    ## # ℹ 2 more variables: bitscore <dbl>, pred <dbl>

## Reference

Ghaddar B, Blaser M, De S. Revisiting the cancer microbiome using PRISM.
bioRxiv 2025

<https://www.biorxiv.org/content/10.1101/2025.01.21.634087v1>
