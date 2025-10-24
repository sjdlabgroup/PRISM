suppressMessages(library(optparse))
suppressMessages(library(ShortRead))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
suppressMessages(library(data.table))

# ---------- Options ----------
option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--data_path"), action="store", help = "path to fastq files"),
  make_option(c("--kraken_path"), action="store", help = "path to Kraken2 executable"),
  make_option(c("--kraken_db_path"), action="store", help = "path to Kraken2 reference database"),
  make_option(c("--seqkit_path"), action="store", help = "path to SeqKit"),
  make_option(c("--star_path"), action="store", help = "path to STAR"),
  make_option(c("--star_genome_dir"), action="store", help = "path to STAR genome index"),
  make_option(c("--model_org_taxids"), action="store", default = 'NA', help = "path to model organism taxids"),
  make_option(c("--blast_path"), action="store", help = "path to blastn"),
  make_option(c("--blast_db_path"), action="store", help = "path to blast reference database"),
  make_option(c("--prism_path"), action="store", help = "path to PRISM directory"),
  make_option(c("--minimap2_path"), action="store", help = "path to Minimap2 executable"),
  make_option(c("--minimap2_index"), action="store", help = "path to Minimap2 index"),
  make_option(c("--min_read_per"), action="store", default = 10^4, help = "minimum reads per X to analyze"),
  make_option(c("--min_uniq_frac"), action="store", default = 5, help = "minimum unique k-mer ratio"),
  make_option(c("--paired"), action="store", default = TRUE, help = "paired-end (T/F)"),
  make_option(c("--max_sample"), action="store", default = 100, help = "max reads to sample"),
  make_option(c("--barcode_only"), action="store", default = FALSE, help = "barcode-only read 1 (T/F)"),
  make_option(c("--fq1_end"), action="store", default = "_1.fastq", help = "FASTQ file suffix for read 1"),
  make_option(c("--fq2_end"), action="store", default = "_2.fastq", help = "FASTQ file suffix for read 2"),
  make_option(c("--threads"), action="store", default = 1, help = "threads"),
  make_option(c("--min_qcovs"), action="store", default = 80, help = "minimum BLAST qcovs"),
  make_option(c("--out_path"), action="store", default = NULL, help = "output directory"),
  make_option(c("--use_custom_db"), action="store", default = TRUE, help = "Create temporary subsetted BLAST database"),
  make_option(c("--custom_db_path"), action="store", default = NULL, help = "temporary subsetted BLAST database directory"),
  make_option(c("--fasta_size_threshold_custom_db"), action="store", default = 5, help = "Minimum query FASTA file size in megabytes for which to create a temporary subsetted BLAST database")
)
opt = parse_args(OptionParser(option_list = option_list))


# ---------- Parameter Validation ----------
required_params <- c("sample", "data_path", "kraken_path", "kraken_db_path",
                     "seqkit_path", "star_path", "star_genome_dir", 
                     "blast_path", "blast_db_path", "prism_path",
                     "minimap2_path", "minimap2_index")

missing_params <- required_params[sapply(required_params, function(p) is.null(opt[[p]]))]
if (length(missing_params) > 0) {
  stop(paste("Missing required parameters:", paste(missing_params, collapse = ", ")))
}

# ---------- Path / Parameter Setup ----------
sample = opt$sample
data_path = opt$data_path
kraken_path = opt$kraken_path
kraken_db_path = opt$kraken_db_path
seqkit_path = opt$seqkit_path
star_path = opt$star_path
star_genome_dir = opt$star_genome_dir
minimap2_path = opt$minimap2_path
minimap2_index = opt$minimap2_index
blast_path = opt$blast_path
blast_db_path = opt$blast_db_path
prism_path = ifelse(str_detect(opt$prism_path, '/$'), opt$prism_path, paste0(opt$prism_path, '/'))
model_org_taxids = ifelse(opt$model_org_taxids == 'NA', paste0(prism_path, opt$model_org_taxids), opt$model_org_taxids)
min_read_per = opt$min_read_per
min_uniq_frac = opt$min_uniq_frac
paired = opt$paired %>% as.logical()
barcode_only = opt$barcode_only %>% as.logical()
max_sample = opt$max_sample
fq1_end = opt$fq1_end
fq2_end = opt$fq2_end
threads = opt$threads
min_qcovs = opt$min_qcovs
use_custom_db = as.logical(opt$use_custom_db)
custom_db_path = opt$custom_db_path
fasta_size_threshold_custom_db = opt$fasta_size_threshold_custom_db
if(is.null(opt$out_path)){
  out_path = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/data/')
  out_path_final = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/')
} else {
  out_path_final = file.path(opt$out_path)
  out_path = paste0(out_path_final, "/data/")
}

# ---------- Create output directory ----------
#if(dir.exists(out_path_final)){system(paste('rm -r', out_path_final))}
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# ---------- Load Functions ----------
source(paste0(prism_path, 'functions.R'))

# ---------- Logging Setup ----------
log_file_path <- paste0(out_path, sample, "_PRISM.log")
log_file <- file(log_file_path, open = "wt")
on.exit({ close(log_file) }, add = TRUE)
log_message("Starting PRISM pipeline")


# ------- Step 1 ------- 
# ------- Initial Kraken2 host depletion, microbial surveillance, and putative microbial read extraction
if (!file.exists(file.path(out_path, paste0(sample, ".kraken.report.txt")))) {
  out <- run_step("Kraken2 analysis", quote(
    prism_kraken(sample, data_path, kraken_path, kraken_db_path, seqkit_path,
                 out_path, min_read_per, min_uniq_frac, paired, 
                 fq1_end, fq2_end, threads)
  ))
  kr_report = out$kr_report; mpa = out$mpa; sp_taxid = out$taxid
} else {
  log_message("Skipping Kraken2 – output exists.")
  
  kr_report <- read.delim(paste0(out_path, sample, '.kraken.report.txt'), header = F, skip = 2) %>% mutate(V8=str_squish(V8))
  mpa <- read.table(paste0(out_path, sample, '.mpa.prism.txt'))
  sp_taxid <- kr_report$V7[str_detect(kr_report$V6, 'S')]
}

# ------- Step 2 ------- 
# ------- Second host depletion with Minimap2 and microbial read extraction
if (!file.exists(file.path(out_path, paste0(sample, "-minimap.sam")))) {
  out <- run_step("Minimap2 host depletion", quote(
    prism_minimap(sample, out_path, minimap2_path, minimap2_index, threads, paired)
  ))
  fa1 = out$fa1; fa2 = out$fa2
} else {
  log_message("Skipping Minimap – output exists.")
  fa1 = ShortRead::readFasta(file.path(out_path, paste0(sample, "_1.fa")))
  if(paired == T){fa2 = ShortRead::readFasta(file.path(out_path, paste0(sample, "_2.fa")))}
}

# ------- Step 3 ------- 
# ------- Third host depletion with STAR and putative microbial read extraction
if (!file.exists(file.path(out_path, paste0('star/', sample, "_Aligned.out.sam")))) {
  out <- run_step("STAR host depletion", quote(
    prism_star(sample, out_path, star_path, star_genome_dir, paired, id_df, fa1, fa2)
  ))
  fa1 = out$fa1; fa2 = out$fa2 
} else {
  log_message("Skipping STAR – output exists.")
  fa1 = ShortRead::readFasta(file.path(out_path, paste0(sample, "_1.fa")))
  if(paired == T){fa2 = ShortRead::readFasta(file.path(out_path, paste0(sample, "_2.fa")))}
}

# ------- Step 4 ------- 
# ------- Sub-sample microbial reads by species taxids
if (!file.exists(file.path(out_path, paste0(sample, "_sub_fa1")))) {
  id_df <- run_step("Sub-sampling by taxid", quote(
    prism_subsample(sample, out_path, fa1, fa2, sp_taxid, mpa, max_sample, paired)
  ))
} else {
  log_message("Skipping subsampling – output exists.")
  sub_fa1 = ShortRead::readFasta(file.path(out_path, paste0(sample, "_sub_fa1")))
  ids = ShortRead::id(sub_fa1) %>% as.character() 
  tax = ids %>% str_remove('.*taxid\\|') %>% as.numeric()
  id_df = data.frame(taxid = tax, id = ids %>% str_remove('\\s.*'))
}

# ------- Step 5 ------- 
# ------- BLAST on sub-sampled putative microbial reads
if (!file.exists(file.path(out_path, paste0(sample, "_sub_fa1-blast.csv")))) {
  run_step("BLAST (subsample)", quote(
    prism_blast(sample, out_path, blast_path, blast_db_path, model_org_taxids, 
                id_df, mode = 'subsample', uit = NULL, kr_report, barcode_only, paired, threads,
                use_custom_db = use_custom_db, custom_db_path = custom_db_path, fasta_size_threshold_custom_db = fasta_size_threshold_custom_db, 
                remove_custom_db = F, prism_path)
  ))
} else {
  log_message("Skipping subsample BLAST – output exists.")
}

# ------- Step 6 ------- 
# ------- Find species with uniquely mappable reads: BLAST multi-mapping analysis 
if (!file.exists(file.path(out_path, paste0(sample, "-final-blast1.csv")))) {
  if((paired == T & !file.exists(file.path(out_path, paste0(sample, "-final-blast2.csv")))) | paired == F){
    uit <- run_step("Multimapping analysis", quote(
      prism_multimapping(out_path, sample, kr_report, mpa, paired, 
                         blast_file_pattern = "-blast.csv", nreads = 10, min_qcovs)
    ))
  }
}

# ------- Step 7 ------- 
# ------- BLAST full data set against human, model organisms, identified species
if (!file.exists(file.path(out_path, paste0(sample, "-final-blast1.csv")))) {
  if((paired == T & !file.exists(file.path(out_path, paste0(sample, "-final-blast2.csv")))) | paired == F){
    run_step("BLAST (full)", quote(
      prism_blast(sample, out_path, blast_path, blast_db_path, model_org_taxids, 
                  id_df, mode = 'full', uit = uit, kr_report, barcode_only, paired, threads,
                  use_custom_db = use_custom_db, custom_db_path = custom_db_path, fasta_size_threshold_custom_db = fasta_size_threshold_custom_db, 
                  remove_custom_db = T, prism_path)
    ))
  }
} else {
  log_message("Skipping full BLAST – output exists.")
}

# ------- Step 8 ------- 
# ------- Filter non-microbial reads 
out <- run_step("Filtering non-microbial reads", quote(
  prism_filter(out_path, sample, kr_report, mpa, paired, blast_file_pattern = 'final-blast', min_qcovs)
))
blast = out$micro; unmapped = out$unmapped

# ------- Step 9 ------- 
# ------- Merge with Genbank gene/product data
prod <- run_step("GenBank annotation", quote(
  prism_genbank(blast)
))

# ------- Step 10 ------- 
# ------- PRISM counts and score prediction
out <- run_step("PRISM profiling", quote(
  prism_profile(out_path, sample, kr_report, blast, prod)
))
counts = out$counts; prod = out$prod; xgmat = out$xgmat

# ------- Step 11 ------- 
# ------- Save results 
run_step("Saving outputs", quote({
  fwrite(prod, file = paste0(file.path(out_path_final, sample), '-results.csv'))
  fwrite(counts, file = paste0(file.path(out_path_final, sample), '-counts.csv'))
  fwrite(xgmat, file = paste0(file.path(out_path, sample), '-xgmat.csv'))
  fwrite(blast, file = paste0(file.path(out_path, sample), '-filtered-blast.csv'))
  fwrite(unmapped, file = paste0(file.path(out_path, sample), '-unmapped.csv'))
}))

# ------- Step 12 ------- 
# ------- Make PRISM FASTA files
run_step("Generating FASTA", quote(
  prism_FASTA(fa1, fa2, blast, out_path, out_path_final, sample, paired)
))
log_message("PRISM pipeline completed successfully.")
