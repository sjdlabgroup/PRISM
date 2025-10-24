# ---------- Logging ----------
log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  if (exists("log_file") && !is.null(log_file)) {
    writeLines(msg, log_file)
    flush(log_file)
  }
}

log_error <- function(...) {
  log_message("ERROR:", ...)
  stop(paste(...))
}

# ---------- Pipeline Execution ----------
run_step <- function(step_name, expr) {
  log_message("Starting:", step_name)
  tryCatch({
    result <- eval(expr)
    log_message("Completed:", step_name)
    return(result)
  }, error = function(e) {
    log_error("Step failed -", step_name, ":", e$message)
  })
}

# prism_kraken.R
#' Run Kraken2 classification and extract microbial reads.
#'
#' This function invokes Kraken2 on paired- or single-end FASTQ(s), generates standard
#' and MPA-style reports, and extracts microbial reads into FASTA outputs.
#'
#' @param sample         Sample name (string)
#' @param data_path      Directory containing sample FASTQ files
#' @param kraken_path    Path to Kraken2 executable
#' @param kraken_db_path Path to Kraken2 database directory
#' @param seqkit_path    Path to SeqKit executable
#' @param out_path       Directory to write outputs (will be created if missing)
#' @param min_read_per   Minimum reads per X to analyze (default 1e4)
#' @param min_uniq_frac  Minimum ratio of unique k-mer to Kraken reads (default 5)
#' @param paired         Logical: TRUE for paired-end, FALSE for single-end
#' @param fq1_end        Suffix for read1 FASTQ (default "_1.fastq")
#' @param fq2_end        Suffix for read2 FASTQ (default "_2.fastq")
#' @param threads        Number of threads to pass to Kraken2 (default 24)

prism_kraken <- function(sample,
                         data_path,
                         kraken_path,
                         kraken_db_path,
                         seqkit_path,
                         out_path,
                         min_read_per,
                         min_uniq_frac,
                         paired,
                         fq1_end,
                         fq2_end,
                         threads) {
  
  #----  Validate executables & paths ----
  # Check that the Kraken2 binary exists and is executable
  if (!file.exists(kraken_path) || !file.access(kraken_path, 1) == 0) {
    stop("kraken2 executable not found or not executable: ", kraken_path)
  }
  # Check that the Kraken2 database directory exists
  if (!file.exists(kraken_db_path)) {
    stop("Kraken2 DB not found: ", kraken_db_path)
  }
  # Check that the SeqKit binary exists and is executable
  if (!file.exists(seqkit_path) || !file.access(seqkit_path, 1) == 0) {
    stop("SeqKit executable not found or not executable: ", seqkit_path)
  }
  
  #----  Validate input FASTQ/FASTA files ----
  fq1 <- file.path(data_path, paste0(sample, fq1_end))
  if (!file.exists(fq1)) {
    stop("Input file #1 not found: ", fq1)
  }
  if (paired) {
    fq2 <- file.path(data_path, paste0(sample, fq2_end))
    if (!file.exists(fq2)) {
      stop("Input file #2 not found: ", fq2)
    }
  }
  
  
  #---- build Kraken2 command ----
  # Construct the base sample path and output prefix
  base <- file.path(data_path, sample)
  out_prefix <- file.path(out_path, sample)
  
  # If paired-end, include both R1 and R2 inputs and --paired flag
  if (paired) {
    cmd <- paste(
      kraken_path,
      "--db", kraken_db_path,
      "--threads", threads,
      "--paired",
      "--use-names",
      "--report-minimizer-data",
      "--classified-out", paste0(out_prefix, "#.fq"),
      "--output", paste0(out_prefix, ".kraken.output.txt"),
      "--report", paste0(out_prefix, ".kraken.report.txt"),
      paste0(base, fq1_end),
      paste0(base, fq2_end)
    )
  } else {
    # Single-end: only R1 and no --paired flag
    cmd <- paste(
      kraken_path,
      "--db", kraken_db_path,
      "--threads", threads,
      "--use-names",
      "--report-minimizer-data",
      "--classified-out", paste0(out_prefix, ".fq"),
      "--output", paste0(out_prefix, ".kraken.output.txt"),
      "--report", paste0(out_prefix, ".kraken.report.txt"),
      paste0(base, fq1_end)
    )
  }
  # Run the Kraken2 classification command
  code <- system(cmd)
  if (code != 0L) stop("Kraken2 failed (exit code ", code, ")")
  
  #---- generate standard & MPA reports ----
  # Create a minimal “standard” report with only key columns
  std_report <- paste0(out_prefix, ".kraken.report.std.txt")
  std_cmd <- paste(
    "cut -f1-3,6-8", paste0(out_prefix, ".kraken.report.txt"),
    ">", std_report
  )
  code <- system(std_cmd)
  if (code != 0L) stop("cut failed (exit code ", code, ")")
  
  # Convert the standard report into MPA-style format
  mpa_cmd <- paste(
    paste0(prism_path, 'kreport2mpa.py'), "-r", std_report,
    "-o", paste0(out_prefix, ".kraken.report.mpa.txt"),
    "--intermediate-ranks"
  )
  code <- system(mpa_cmd)
  if (code != 0L) stop("kreport2mpa.py failed (exit code ", code, ")")
  
  #---- extract microbial reads ----
  # Path to the original Kraken2 output
  kraken_out   <- paste0(out_prefix, ".kraken.output.txt")
  
  #---- read back reports & determine which species to keep ----
  kr = paste0(out_path, sample, '.kraken.report.txt')
  mpa = paste0(out_path, sample, '.kraken.report.mpa.txt')
  
  # Load report and drop non-microbial entries
  kr = read.delim(kr, header = F) %>% mutate(V8 = trimws(V8))
  kr = kr[-which(kr$V8 %in% c('root', 'unclassified')), ] %>%
    mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
  
  # Load MPA report and initialize taxid column
  mpa = read.delim(mpa, header = F)
  mpa$taxid = NA
  
  # For each MPA line, extract the taxid(s) matching the report entries
  for(i in 2:nrow(mpa)){
    t_names = mpa[i,1] %>% as.character() %>%
      strsplit('\\|') %>% unlist() %>%
      str_remove('.*__') %>%
      str_replace_all('[^[:alnum:]]', '_')
    mpa$taxid[i] = paste0('*',
                          paste(kr$V7[match(t_names, kr$V8)], collapse = '*'),
                          '*')
  }
  
  write.table(mpa, file = paste0(out_prefix, '.mpa.prism.txt'))
  
  # Identify microbiome entries (bacteria/fungi/viruses) at species level
  n = str_which(mpa$V1, '__Bacteria|__Fungi|__Viruses')
  
  # Minimum reads threshold based on host counts
  min_reads = ceiling(sum(kr$V2[kr$V7 %in% c(2,4751,10239)])/min_read_per)
  min_reads = ifelse(min_reads > 9, min_reads, 10)
  
  # Compute unique‐fraction for each species
  kr$uf = kr$V5/kr$V2
  kr2 = kr  # backup full set
  kr = kr[n, ] %>% subset((V2 > min_reads | (uf > min_uniq_frac & V2 > 9)))
  sp_taxids = kr$V7[str_detect(kr$V6, 'S')]
  
  #---- filter the Kraken output for only microbial reads ----
  mic_out <- paste0(out_prefix, ".kraken.microbiome.output.txt")
  taxids_to_remove = c(0,1, setdiff(kr2$V7, kr2$V7[n])) # remove nonmicrobial Kraken reads
  
  # Use awk to drop any line whose taxid is in the “remove” list
  awk_cmd <- sprintf(
    "awk -F '\t' '$3 !~ /(%s)/' %s >> %s",
    paste0('\\(taxid ', taxids_to_remove, '\\)', collapse='|'),
    kraken_out, mic_out
  )
  code <- system(awk_cmd)
  if (code != 0L) stop("awk filtering failed (exit code ", code, ")")
  
  # delete full kraken output file
  system(paste0('rm ', out_path, sample, '.kraken.output.txt'))
  
  #---- extract read IDs for downstream FASTA conversion ----
  id_file <- paste0(out_prefix, '.kraken.microbiome.ids.txt')
  cut_cmd <- paste(
    "cut -d '\t' -f2", mic_out, ">", id_file
  )
  code <- system(cut_cmd)
  if (code != 0L) stop("cut IDs failed (exit code ", code, ")")
  
  #---- pull microbial reads out of FASTQ and convert to FASTA ----
  
  # Determine if FASTQ headers contain /1 or /2 (only read 1 line)
  fq_check <- if (paired) paste0(out_prefix, "_1.fq") else paste0(out_prefix, ".fq")
  fq_header <- system(paste("head -n 1", shQuote(fq_check)), intern = TRUE)
  has_suffix <- grepl("/[12]\\b", fq_header)
  
  # If suffix is present, modify id_file in place to include /1 and /2
  if (has_suffix) {
    tmp_ids <- paste0(id_file, ".tmp")
    system(paste0("awk '{print $0\"/1\"; print $0\"/2\"}' ", shQuote(id_file), " > ", shQuote(tmp_ids)))
    file.rename(tmp_ids, id_file)
  }
  
  if (paired) {
    f1 <- paste0(out_prefix, "_1.fq"); f2 <- paste0(out_prefix, "_2.fq")
    tmp1 <- paste0(out_prefix, "_1.fq-temp")
    tmp2 <- paste0(out_prefix, "_2.fq-temp")
    
    cmd1 <- paste(seqkit_path, "grep --pattern-file", shQuote(id_file), shQuote(f1), ">", shQuote(tmp1))
    system(cmd1)
    system(paste(seqkit_path, "fq2fa", shQuote(tmp1), ">", shQuote(paste0(out_prefix, "_1.fa"))))
    file.remove(tmp1, f1)
    
    cmd2 <- paste(seqkit_path, "grep --pattern-file", shQuote(id_file), shQuote(f2), ">", shQuote(tmp2))
    system(cmd2)
    system(paste(seqkit_path, "fq2fa", shQuote(tmp2), ">", shQuote(paste0(out_prefix, "_2.fa"))))
    file.remove(tmp2, f2)
    
  } else {
    f1 <- paste0(out_prefix, ".fq")
    tmp1 <- paste0(out_prefix, "_1.fq-temp")
    
    cmd1 <- paste(seqkit_path, "grep --pattern-file", shQuote(id_file), shQuote(f1), ">", shQuote(tmp1))
    system(cmd1)
    system(paste(seqkit_path, "fq2fa", shQuote(tmp1), ">", shQuote(paste0(out_prefix, "_1.fa"))))
    file.remove(tmp1, f1)
  }
  
  
  # delete temporary microbiome ids txt file
  system(paste0('rm ', out_path, sample, '.kraken.microbiome.ids.txt'))
  
  # Return the MPA table and the filtered Kraken report
  list(mpa = mpa, kr_report = kr2, taxid = sp_taxids)
}


#' Run minimap2 on Kraken-extracted FASTA(s) and produce a SAM file
#'
#' After you've run `run_kraken_extract()`, this function locates the
#' `${sample}_1.fa` (and optional `${sample}_2.fa`) in `out_path` and aligns
#' them to a provided `.mmi` index using **minimap2**. It writes a SAM file
#' (`${sample}-minimap.sam`) and then parses it to identify unmapped reads,
#' rewriting the input FASTA file(s) to contain only those unmapped reads.
#'
#' @param sample String. Sample name (same prefix used in `run_kraken_extract()`).
#' @param out_path String. Directory where `${sample}_1.fa` (and `${sample}_2.fa`, if paired) reside and where outputs are written.
#' @param minimap2_path String. Full path to the `minimap2` executable.
#' @param minimap2_index String. Full path to the `.mmi` index file.
#' @param threads Integer. Number of threads for minimap2 (default: 8).
#' @param paired Logical. `TRUE` if a `${sample}_2.fa` mate file is expected (default: `FALSE`).

prism_minimap = function(sample,
                         out_path,
                         minimap2_path,
                         minimap2_index,
                         threads,
                         paired) {
  
  # build paths to inputs & output
  fasta1  <- file.path(out_path, paste0(sample, "_1.fa"))
  fasta2  <- if (paired == T) file.path(out_path, paste0(sample, "_2.fa")) else NULL
  out_sam <- file.path(out_path, paste0(sample, "-minimap.sam"))
  
  #---- 1) PARAMETER VALIDATION ----
  if (!file.exists(minimap2_path) || !file.access(minimap2_path, 1) == 0) {
    stop("minimap2 executable not found or not executable: ", minimap2_path)
  }
  if (!file.exists(minimap2_index)) {
    stop("minimap2 index not found: ", minimap2_index)
  }
  if (!file.exists(fasta1)) {
    stop("Extracted FASTA #1 not found: ", fasta1)
  }
  if (paired && !file.exists(fasta2)) {
    stop("Extracted FASTA #2 not found (paired=TRUE): ", fasta2)
  }
  
  #---- 2) BUILD & RUN COMMAND ----
  cmd <- if (paired) {
    sprintf('%s -ax sr -t %d %s %s %s > %s',
            shQuote(minimap2_path), threads,
            shQuote(minimap2_index),
            shQuote(fasta1), shQuote(fasta2),
            shQuote(out_sam))
  } else {
    sprintf('%s -ax sr -t %d %s %s --no-pairing > %s',
            shQuote(minimap2_path), threads,
            shQuote(minimap2_index),
            shQuote(fasta1),
            shQuote(out_sam))
  }
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] running Minimap2: ", cmd)
  status <- system(cmd)
  
  #---- 3) CHECK EXIT STATUS ----
  if (status != 0L) {
    stop("minimap2 failed with exit code ", status)
  }
  
  #---- 4) VERIFY SAM ----
  if (!file.exists(out_sam)) {
    stop("SAM file was not created: ", out_sam)
  }
  info <- file.info(out_sam)
  if (is.na(info$size) || info$size == 0) {
    stop("SAM file is empty: ", out_sam)
  }
  
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] SAM written to: ", out_sam)
  
  
  #---- 5) extract unmapped reads and rewrite FASTA(s) ----
  # 1) read the SAM and grab query IDs where the 0x4 (unmapped) bit is set
  sam_df <- tryCatch(
    suppressWarnings(
      read.delim(out_sam, comment.char='@', header=FALSE, stringsAsFactors = F) %>% mutate(V2 = as.numeric(V2))
    ),
    error = function(e) stop("Failed to read SAM: ", out_sam, "\n", e$message)
  )
  unmapped_ids <- sam_df$V1[ bitwAnd(sam_df$V2, 4) != 0 ]
  
  if (length(unmapped_ids)==0) {
    message("No non-human reads found in ", out_sam)
  } else {
    # 2) load the input FASTA(s) you aligned
    fa1 <- ShortRead::readFasta(fasta1)
    fa1_ids <- ShortRead::id(fa1) %>% as.character() %>% str_remove('\\s.*')
    fa1  <- fa1[fa1_ids %in% unmapped_ids]
    out1   <- file.path(out_path, paste0(sample, "_1.fa"))
    writeFasta(fa1, file=out1)
    message("Wrote ", length(fa1), " unmapped reads to ", out1)
    
    # 3) if paired, do the same for read2
    if (!is.null(fasta2)) {
      fa2 <- ShortRead::readFasta(fasta2)
      fa2_ids <- ShortRead::id(fa2) %>% as.character() %>% str_remove('\\s.*')
      fa2   <- fa2[fa2_ids %in% unmapped_ids]
      out2    <- file.path(out_path, paste0(sample, "_2.fa"))
      writeFasta(fa2, file=out2)
      message("Wrote ", length(fa2), " unmapped reads to ", out2)
    } else {fa2 = NULL}
  }
  
  list(fa1 = fa1, fa2 = fa2)
}


#' Run STAR to remove host reads from FASTA, keeping only unmapped reads
#'
#' This function aligns either a single FASTA file (`_1.fa`) or paired FASTA files
#' (`_1.fa` + `_2.fa`) against a specified STAR genome index, then filters out
#' all reads that map to the host, retaining only unmapped reads.
#'
#' @param sample String. Sample name used to construct input and output file paths.
#' @param out_path String. Directory where `<sample>_1.fa` (and optionally `<sample>_2.fa`) are located.
#' @param star_path String. Full path to the STAR executable.
#' @param star_genome_dir String. Full path to the STAR genome index directory.
#' @param paired Logical. `TRUE` if paired-end FASTA files (`_1.fa` and `_2.fa`) are present.
#' @param id_df Data frame. Contains read IDs in column `id`, used to annotate with a `star` flag.

prism_star = function(sample,
                      out_path,
                      star_path,
                      star_genome_dir,
                      paired,
                      id_df,
                      fa1,
                      fa2) {
  
  #— basic sanity checks —
  if (!file.exists(star_path) || !file.access(star_path, 1) == 0) {
    stop("STAR not found or not executable: ", star_path)
  }
  if (!dir.exists(star_genome_dir)) {
    stop("STAR genomeDir not found: ", star_genome_dir)
  }
  # input FASTA paths
  fa1_path <- file.path(out_path, paste0(sample, "_1.fa"))
  if (!file.exists(fa1_path)) stop("FASTA #1 not found: ", fa1_path)
  if (paired) {
    fa2_path <- file.path(out_path, paste0(sample, "_2.fa"))
    if (!file.exists(fa2_path)) stop("FASTA #2 not found: ", fa2_path)
  }
  
  #— make STAR output directory —
  star_dir <- file.path(out_path, "star")
  if (!dir.exists(star_dir)) dir.create(star_dir, recursive = TRUE)
  
  #— build & run STAR —
  prefix <- file.path(star_dir, paste0(sample, "_"))
  cmd <- if (paired == T) {
    paste(
      star_path,
      "--genomeDir", shQuote(star_genome_dir),
      "--outFileNamePrefix", shQuote(prefix),
      "--readFilesIn",
      shQuote(fa1_path), shQuote(fa2_path)
    )
  } else {
    paste(
      star_path,
      "--genomeDir", shQuote(star_genome_dir),
      "--outFileNamePrefix", shQuote(prefix),
      "--readFilesIn", shQuote(fa1_path)
    )
  }
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] running STAR: ", cmd)
  status <- system(cmd)
  if (status != 0) stop("STAR failed (exit code ", status, ")")
  
  #— load the resulting SAM, but bail out if empty —
  sam_file <- file.path(out_path, "star", paste0(sample, "_Aligned.out.sam"))
  if (!file.exists(sam_file) || file.info(sam_file)$size == 0) {
    message("No alignments found (empty SAM); skipping STAR‐based filtering.")
    # return FASTAs untouched
    fa1_keep <- ShortRead::readFasta(file.path(out_path, paste0(sample, "_1.fa")))
    if (paired) {
      fa2_keep <- ShortRead::readFasta(file.path(out_path, paste0(sample, "_2.fa")))
      return(list(fa1 = fa1_keep, fa2 = fa2_keep))
    }
    return(list(fa1 = fa1_keep))
  }
  
  sam_df <- tryCatch(
    read.delim(sam_file, comment.char='@', header=FALSE, stringsAsFactors=FALSE),
    error = function(e) {
      message("Could not parse SAM (perhaps no alignments); skipping filtering.")
      return(NULL)
    }
  )
  if (is.null(sam_df) || nrow(sam_df) == 0) {
    message("No human mapped reads in SAM; skipping STAR‐based filtering.")
    fa1_keep <- ShortRead::readFasta(file.path(out_path, paste0(sample, "_1.fa")))
    if (paired) {
      fa2_keep <- ShortRead::readFasta(file.path(out_path, paste0(sample, "_2.fa")))
      return(list(fa1 = fa1_keep, fa2 = fa2_keep))
    }
    return(list(fa1 = fa1_keep))
  }
  
  #— get FASTA header ids —
  fa1_ids <- as.character(ShortRead::id(fa1)) %>% sub("\\s.*", "", .)
  if (paired == T) {
    fa2_ids <- as.character(ShortRead::id(fa2)) %>% sub("\\s.*", "", .)
  }
  
  #— identify mapped reads (FLAG &4 means unmapped → keep; so mapped = bitwAnd(flag,4)==0) —
  mapped_ids <- sam_df$V1[bitwAnd(sam_df$V2, 4) == 0]
  # id_df$star <- ifelse(id_df$id %in% mapped_ids, 1, 0)
  
  #— subset to only unmapped reads and write back FASTA —
  fa1_keep <- fa1[ !(fa1_ids %in% mapped_ids) ]
  ShortRead::writeFasta(fa1_keep, file = fa1_path)
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] wrote filtered FASTA: ", fa1_path)
  
  if (paired == T) {
    fa2_keep <- fa2[ !(fa2_ids %in% mapped_ids) ]
    ShortRead::writeFasta(fa2_keep, file = fa2_path)
    message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] wrote filtered FASTA: ", fa2_path)
    return(list(fa1 = fa1_keep, fa2 = fa2_keep))
  } else {
    return(list(fa1 = fa1_keep))
  }
}

#' Subsample microbial reads by taxon
#'
#' This function takes Kraken-extracted FASTA files (`_1.fa` and optionally `_2.fa`)
#' and subsamples microbial reads by taxon ID. For each taxon in `taxid`, it selects
#' up to `max_sample` reads, prioritizing species-level reads first and then including
#' higher or lower taxonomic levels (from `mpa`) as needed to reach the desired number.
#'
#' The resulting subsampled FASTA(s) are written as `<sample>_sub_fa1` (and
#' `<sample>_sub_fa2` if paired).
#'
#' @param sample String. Sample name used to construct input and output FASTA paths.
#' @param out_path String. Directory containing the input FASTA(s) and where output files will be written.
#' @param fa1 ShortRead object. FASTA reads from `${sample}_1.fa`.
#' @param fa2 ShortRead object (optional). FASTA reads from `${sample}_2.fa` if paired-end.
#' @param taxid Numeric or character vector. Taxon IDs to target for subsampling.
#' @param mpa Data frame. Mapping of taxonomic lineages, typically parsed from a MetaPhlAn or similar file, containing a `taxid` column.
#' @param max_sample Integer. Maximum number of reads to retain per taxon.
#' @param paired Logical. `TRUE` if paired-end FASTA files (`_1.fa` and `_2.fa`) are present.

prism_subsample = function(sample,
                           out_path,
                           fa1,  
                           fa2,
                           taxid,
                           mpa,
                           max_sample,
                           paired) {
  
  #--- construct input FASTA paths ---
  fa1_file <- file.path(out_path, paste0(sample, "_1.fa"))
  if (!file.exists(fa1_file)) {
    stop("Input FASTA #1 not found: ", fa1_file)
  }
  if (paired) {
    fa2_file <- file.path(out_path, paste0(sample, "_2.fa"))
    if (!file.exists(fa2_file)) {
      stop("Input FASTA #2 not found: ", fa2_file)
    }
  }
  
  
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Started subsampling for ", sample)
  
  #--- read in FASTA #1 and extract headers ---
  headers = ShortRead::id(fa1)
  
  #--- parse taxids from headers (in 1k‐read chunks if large) ---
  if(length(headers) > 1000){
    header_taxid = list()
    i1 = seq(0, length(headers), by = 1000)[-1] 
    i2 = i1 - 999  
    if(max(i1) < length(fa1)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa1))}
    for(i in 1:length(i2)){
      # print(i)
      header_taxid[[i]] = str_remove(headers[i2[i]:i1[i]] %>% as.character(), '.*taxid\\|')
    }
    header_taxid = unlist(header_taxid)
  } else {header_taxid = str_remove(headers %>% as.character(), '.*taxid\\|')}
  
  #--- initial sampling: find indices of each taxid in header_taxid ---
  #--- limit to max_sample per taxid ---
  n = which(header_taxid %in% taxid)
  micro_headers = headers[n]
  micro_header_taxid = header_taxid[n]
  n = sapply(taxid %>% as.character(), function(x) which(micro_header_taxid %in% x))
  for(i in 1:length(n)){if(length(n[[i]]) > max_sample){n[[i]] = n[[i]][1:max_sample]}}
  if(class(n) != 'list'){n=list(as.vector(n)); names(n)=taxid}
  micro_headers = lapply(n, function(x) micro_headers[x])
  
  # adjust subsampling using taxon taxonomy for undersampled reads if possible 
  for(i in 1:length(micro_headers)){
    # print(i)
    ## include higher resolution reads (e.g. subspecies, strain, etc.) if not enough species level reads are sampled
    if(length(micro_headers[[i]]) < max_sample){
      lin = str_subset(mpa$taxid, paste0('\\*', names(n)[i], '\\*')) %>% 
        str_extract(paste0('\\*', names(n)[i], '\\*.*')) %>%
        str_remove('^\\*') %>% 
        str_remove('\\*$') %>% 
        str_split('\\*') %>% 
        unlist() %>%
        str_subset('NA', negate = T) %>% 
        as.numeric() %>%
        unique()
      micro_headers[[i]] = headers[which(header_taxid %in% lin)[1:max_sample] %>% na.omit()]
    }
    
    ## include lower resolution reads (e.g. genus/family/class etc.) if not enough species level reads are sampled
    if(length(micro_headers[[i]]) < max_sample){
      lin = str_subset(mpa$taxid, paste0('\\*', names(n)[i], '\\*')) %>% 
        str_remove('^\\*') %>% 
        str_remove('\\*$') %>% 
        str_split('\\*') %>% 
        unlist() %>%
        str_subset('NA', negate = T) %>% 
        as.numeric() %>%
        unique()
      idx = which(lin %in% c(2, 4751, 10239))
      lin = lin[idx:length(lin)]
      counter = which(lin == names(micro_headers)[i])
      while(length(micro_headers[[i]]) < max_sample & counter > 0){
        counter = counter - 1
        tx = lin[counter]
        micro_headers[[i]] = c(micro_headers[[i]], headers[which(header_taxid == tx)[1:max_sample] %>% na.omit()])
        if(length(micro_headers[[i]]) > max_sample){micro_headers[[i]] = micro_headers[[i]][1:max_sample]}
      }
    }
  }
  
  # extract subsampled fastas
  micro_headers = micro_headers %>% lapply(as.character) %>% unlist() %>% unname() %>% unique()
  micro_ids = micro_headers %>% str_remove('\\s.*')
  fa1 = subset(fa1, headers %in% micro_headers)
  if(file.exists(paste0(out_path, sample, '_sub_fa1'))){system(paste('rm', paste0(out_path, sample, '_sub_fa1')))}
  writeFasta(fa1, file = paste0(out_path, sample, '_sub_fa1'))
  
  #--- build ID‐to‐taxid table ---
  ids = ShortRead::id(fa1) %>% as.character() 
  tax = ids %>% str_remove('.*taxid\\|') %>% as.numeric()
  id_df = data.frame(taxid = tax, id = ids %>% str_remove('\\s.*'))
  
  #--- process FASTA #2 similarly ---
  if(paired == T){
    # fa2
    headers2 = ShortRead::id(fa2)
    
    if(length(headers2) > 1000){
      ids2 = list()
      i1 = seq(0, length(headers2), by = 1000)[-1] 
      i2 = i1 - 999  
      if(max(i1) < length(fa2)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa2))}
      for(i in 1:length(i2)){
        # print(i)
        ids2[[i]] = str_remove(headers2[i2[i]:i1[i]] %>% as.character(), '\\s.*')
      }
      ids2 = unlist(ids2)
    } else {ids2 = str_remove(headers2 %>% as.character(), '\\s.*')}
    
    fa2 = subset(fa2, str_remove(ids2,'/1$|/2$') %in% str_remove(micro_ids, '/1$|/2$')) # remove read mate identifiers if present
    if(file.exists(paste0(out_path, sample, '_sub_fa2'))){system(paste('rm', paste0(out_path, sample, '_sub_fa2')))}
    writeFasta(fa2, file = paste0(out_path, sample, '_sub_fa2'))
  } 
  
  id_df
}


#' Build a custom BLAST database from selected taxids
#'
#' This function creates a **custom BLAST nucleotide database** containing only
#' sequences corresponding to a specified set of taxonomic IDs. It extracts those
#' sequences from an existing BLAST database, maps accessions to taxids, and builds
#' a new, smaller BLAST database suitable for targeted searches.
#'
#' @param blast_path String. Path to the directory containing BLAST+ executables
#' (e.g., `blastdbcmd`, `makeblastdb`).
#' @param blast_db_path String. Full path to the source BLAST database (without file extensions).
#' @param uit Numeric or character vector. List of taxonomic IDs to include in the custom database.
#' @param custom_db_path String. Output directory where the custom BLAST database and intermediate files will be written.
#' @param out_path String. Directory for any additional downstream outputs (reserved for pipeline consistency).
#' @param full_accession_map_path String. Path to a full accession-to-taxid mapping file used to construct the taxid map for the custom database.

prism_custom_blastdb <- function(blast_path, blast_db_path, uit, custom_db_path, out_path, full_accession_map_path) {
  
  #---- Step 1: Validate Inputs ----
  if (length(uit) == 0) stop("No taxids provided to generate custom BLAST database.")
  if (!dir.exists(custom_db_path)) {
    dir.create(custom_db_path, recursive = TRUE)
  }
  
  #---- Step 2: Define Paths ----
  taxid_file      <- file.path(custom_db_path, "uit_taxids.txt")
  custom_fasta    <- file.path(custom_db_path, "customdb.fasta")
  taxid_map_file  <- file.path(custom_db_path, "taxid_map.txt")
  blast_db_prefix <- file.path(custom_db_path, "custom_taxa_db")
  headers_file    <- file.path(custom_db_path, "headers.txt")
  
  #---- Step 3: Write taxid list ----
  writeLines(as.character(uit), taxid_file)
  
  #---- Step 4: FASTA extraction via taxidlist (faster) ----
  cmd_extract <- sprintf(
    "%s/blastdbcmd -db %s -dbtype nucl -outfmt %%f -taxidlist %s -out %s",
    shQuote(blast_path),
    shQuote(blast_db_path),
    shQuote(taxid_file),
    shQuote(custom_fasta)
  )
  system(cmd_extract)
  
  #---- Step 5: Build taxid_map from master accession-taxid map ----
  system(paste("grep '^>'", shQuote(custom_fasta), ">", shQuote(headers_file)))
  headers <- readLines(headers_file)
  accessions <- unique(str_match(headers, ">([^\\t\\s]+)\\s*")[,2])
  
  # Sort accessions in R and write directly to the sorted file
  sorted_accession_file <- file.path(custom_db_path, "accessions_sorted.txt")
  writeLines(sort(accessions), sorted_accession_file)
  
  # Perform the join: accession (col 1) to accession_map (col 2)
  join_cmd <- sprintf(
    "join -1 1 -2 2 -t $'\\t' %s %s > %s",
    shQuote(sorted_accession_file),
    shQuote(full_accession_map_path),
    shQuote(taxid_map_file)
  )
  system(join_cmd)
  
  #---- Step 6: Build custom BLAST DB ----
  cmd_make <- sprintf(
    "%s/makeblastdb -in %s -dbtype nucl -parse_seqids -blastdb_version 5 -taxid_map %s -out %s",
    shQuote(blast_path),
    shQuote(custom_fasta),
    shQuote(taxid_map_file),
    shQuote(blast_db_prefix)
  )
  system(cmd_make)
  
  # Clean up intermediate files
  file.remove(taxid_file, taxid_map_file, custom_fasta)
}


#' Run BLAST on FASTA(s)
#'
#' Aligns subsampled FASTA files against a BLAST database that includes human, common model organism, and pre-identified microbial species taxids
#'
#' @param sample              Sample name (string)
#' @param out_path            Directory where `_sub_fa1` (and `_sub_fa2`) live and where outputs will go
#' @param blast_path          Path to the `blast` executables that contains `blastn`, `blastdbcmd`, `makeblastdb`
#' @param blast_db_path       Path to the BLAST database
#' @param model_org_taxids    Path to a text file of model-organism taxids (one per line)
#' @param id_df               Data.frame with at least columns `taxid` and `star`
#' @param mode                Must be either "subsample" or "full"
#' @param uit                 uniquely identified species taxids
#' @param kr2                 Data.frame from Kraken report with columns `V6` (rank) and `V7` (taxid)
#' @param barcode_only        Logical; if TRUE, skip read1 BLAST
#' @param paired              Logical; if TRUE, also BLAST `_sub_fa2`
#' @param make_custom_db      Logica; if TRUE, create a temporary subsetted BLAST database with only the taxids of interest. Speeds up code but may require up 5-100GB of disk space
#' @param custom_db_path      Output directory for subsetted BLAST database if it is created
#' @param fasta_size_threshold_custom_db minimun fasta query file size in megabytes to create a custom blast database
#'
#' @return Invisibly returns a character vector of the BLAST output file paths.
#' @export

prism_blast = function(sample,
                       out_path,
                       blast_path,
                       blast_db_path,
                       model_org_taxids,
                       id_df,
                       mode,
                       uit = NULL,
                       kr_report,
                       barcode_only,
                       paired,
                       threads, 
                       use_custom_db = T,
                       custom_db_path = NULL,
                       fasta_size_threshold_custom_db = 5,
                       remove_custom_db = F,
                       prism_path) {
  
  #---- check inputs ----
  if (!file.exists(model_org_taxids)) {
    stop("Model‐organism taxid file not found: ", model_org_taxids)
  }
  if (!dir.exists(blast_path) || file.access(blast_path, 4) != 0) {
    stop("blastn directory not found or not accessible: ", blast_path)
  }
  if (!dir.exists(out_path)) {
    stop("Output directory not found: ", out_path)
  }
  
  if(mode == 'subsample'){
    message("Started BLAST on subsampled FASTA1 at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  } else {
    message("Started BLAST on full FASTA1 at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  }
  
  
  #---- read model‐organism taxids ----
  model_orgs <- read.delim(model_org_taxids, header = FALSE, stringsAsFactors = FALSE)[[1]] %>% as.character()
  
  # ── sanity checks ─────────────────────────────────────────────────────────────
  
  # 1) did we actually read anything?
  if (length(model_orgs) == 0L) {
    stop("No taxids found in ", shQuote(model_org_taxids), " – file is empty.")
  }
  
  # 2) are they all integer strings?
  invalid <- model_orgs[!grepl("^[0-9]+$", model_orgs)]
  if (length(invalid) > 0L) {
    stop("Invalid taxid(s) in ", shQuote(model_org_taxids), ": ",
         paste(invalid, collapse = ", "))
  }
  
  #---- determine the taxids to BLAST ----
  if(mode == 'subsample'){
    #---- select taxids with >=10 non‐human reads and species rank in kr_report ----
    tx_table <- table(id_df$taxid)
    tx_keep  <- names(tx_table)[tx_table > 9]
    species_taxids <- kr_report$V7[str_detect(kr_report$V6, "S")]
    tx <- intersect(tx_keep, species_taxids)
    
    FA1 = shQuote(file.path(out_path, paste0(sample, "_sub_fa1")))
    FA2 = shQuote(file.path(out_path, paste0(sample, "_sub_fa2")))
    out1 = file.path(out_path, paste0(sample, "_sub_fa1-blast.csv"))
    out2 = file.path(out_path, paste0(sample, "_sub_fa2-blast.csv"))
    
  } else if(mode == 'full') {
    
    # make sure species taxids are provide in uit
    if (length(uit) == 0L) {stop("No taxids provided")}
    tx = uit
    FA1 = shQuote(file.path(out_path, paste0(sample, "_1.fa")))
    FA2 = shQuote(file.path(out_path, paste0(sample, "_2.fa")))
    out1 = file.path(out_path, paste0(sample, "-final-blast1.csv"))
    out2 = file.path(out_path, paste0(sample, "-final-blast2.csv"))
    
  }
  
  # create custom BLAST database only with taxa of interest -- increases code speed
  if(use_custom_db == T & file.size(file.path(paste0(out_path,sample, '_1.fa'))) > fasta_size_threshold_custom_db*10^6){
    if(is.null(custom_db_path)){custom_db_path = paste0(out_path, '/customdb')} # if no custom database path is provided, make one
    if(!dir.exists(custom_db_path)){dir.create(custom_db_path, recursive = T)}
    if(length(list.files(custom_db_path)) == 0){ # if the custom database folder is empty, make the database 
      prism_custom_blastdb(blast_path, blast_db_path, uit = c(9606, model_orgs, tx), custom_db_path, out_path, full_accession_map_path = file.path(prism_path, 'sorted_accession_map.txt'))
    }
  }
  
  # if a custom subsetted blast database path is supplied, use it instead of the full blast database
  if(!is.null(custom_db_path)){blast_db_path = file.path(custom_db_path, "custom_taxa_db")}
  
  # if a custom database was already made in the default location (custom_db_path = NULL) use it
  if(dir.exists(paste0(out_path, '/customdb')) & length(list.files(paste0(out_path, '/customdb'))) > 0){
    custom_db_path = paste0(out_path, '/customdb')
    blast_db_path = file.path(custom_db_path, "custom_taxa_db")
  }
  
  setwd(str_remove(blast_db_path, basename(blast_db_path)))
  
  #---- BLAST read1 ----
  if (barcode_only == F) {
    cmd1 <- paste0(
      blast_path, "/blastn ",
      "-max_hsps 1 ",
      "-culling_limit 10 ",
      "-max_target_seqs 10 ",
      "-task megablast ",
      "-num_threads ", threads, " ", 
      "-query ",  FA1, " ",
      "-db ",     shQuote(blast_db_path), " ",
      '-outfmt "6 delim=, qseqid sacc staxids sstart qcovs pident bitscore" ',
      "-out ",    shQuote(out1), " ",
      "-taxids ", paste0(c(9606, model_orgs, tx), collapse = ",")
    )
    code1 <- system(cmd1, ignore.stderr = TRUE)
    if (code1 != 0L){stop("BLAST read1 failed (exit code ", code1, "):\n", cmd1)}
  }
  
  #---- BLAST read2 if paired ----
  if (paired == T) {
    
    if(mode == 'subsample'){
      message("Started BLAST on subsampled FASTA2 at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    } else {
      message("Started BLAST on full FASTA2 at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    }
    
    cmd2 <- paste0(
      blast_path, "/blastn ",
      "-max_hsps 1 ",
      "-culling_limit 10 ",
      "-max_target_seqs 10 ",
      "-task megablast ",
      "-num_threads ", threads, " ", 
      "-query ",  FA2, " ",
      "-db ",     shQuote(blast_db_path), " ",
      '-outfmt "6 delim=, qseqid sacc staxids sstart qcovs pident bitscore" ',
      "-out ",    shQuote(out2), " ",
      "-taxids ", paste0(c(9606, model_orgs, tx), collapse = ",")
    )
    code2 <- system(cmd2, ignore.stderr = TRUE)
    if (code2 != 0L){stop("BLAST read2 failed (exit code ", code2, "):\n", cmd2)}
  }
  
  # remove custom subsetted blast database 
  if(remove_custom_db == T & !is.null(custom_db_path)){system(paste('rm -r', custom_db_path))}
  
}


#' Filter and classify BLAST-aligned reads to retain high-quality microbial hits
#'
#' This function filters BLAST alignment results from a set of CSV files to retain only
#' high-quality, non-human, microbial reads. It removes human, low-quality, and
#' non-microbial hits, assigns taxonomy to each retained read at the most specific
#' (typically species) rank, and returns both the filtered microbial reads and
#' unmapped/filtered-out reads.
#'
#' @param out_path String. Directory containing the BLAST result CSV files (e.g., `*-blast.csv`).
#' @param sample String. Sample name, used for consistency in downstream processing.
#' @param kr_report Data frame or data.table. Kraken report table containing columns
#' for taxonomic hierarchy, taxonomic rank, name, and taxid.
#' @param mpa Data frame. MetaPhlAn-like taxonomic profile, used to identify microbial taxa.
#' @param paired Logical. `TRUE` if paired-end data were processed; used for bookkeeping only.
#' @param blast_file_pattern String. File name pattern (e.g., `"blast.csv"`) used to locate BLAST output files in `out_path`.
#' @param min_qcovs Numeric. Minimum query coverage threshold used to filter BLAST alignments.

prism_filter <- function(out_path, sample, kr_report, mpa, paired, blast_file_pattern, min_qcovs) {
  
  # ————————————————
  # STEP 1: Load BLAST files
  # ————————————————
  
  # List files matching the given pattern in the output path
  blast_files <- list.files(out_path, pattern = blast_file_pattern, full.names = TRUE)
  
  # Error check: ensure BLAST files exist
  if (length(blast_files) == 0) {
    stop("No '*-blast.csv' files found in ", out_path)
  }
  
  # Load all BLAST files and combine them
  # Assign column names and ensure taxid field is valid numeric string
  blast_list <- lapply(seq_along(blast_files), function(i) {
    b <- data.table::fread(blast_files[i], header = FALSE)
    setnames(b, c('id','sacc','staxids','pos','qcovs','pident','bitscore'))  # 7 columns
    b[, read := i]  # Now add 8th column
    return(b)
  })
  blast <- rbindlist(blast_list)
  setcolorder(blast, c("id", "read", "sacc", "staxids", "pos", "qcovs", "pident", "bitscore"))
  
  blast = subset(blast, qcovs > min_qcovs)
  
  # ————————————————————————————————————
  # STEP 2: Remove reads mapping to human (taxid 9606)
  # ————————————————————————————————————
  
  # Make sure staxids is a character only once
  blast[, staxids := as.character(staxids)]
  
  # Flag human-mapped reads (taxid 9606) using efficient group logic
  blast[, is_human := any(staxids == "9606"), by = id]
  
  # Extract human-mapped reads
  unmapped <- blast[is_human == TRUE][, `:=`(
    type = "human"
  )]
  
  # Keep only non-human reads
  blast <- blast[is_human == FALSE]
  
  # Remove temporary flag
  blast[, is_human := NULL]
  
  # Clean up taxid column (remove semicolon and everything after)
  blast[, staxids := sub(";.*", "", staxids)]
  
  # ————————————————————————————————————————————
  # STEP 3: Filter low-quality and non-microbial reads
  # ————————————————————————————————————————————
  
  # Filter microbial taxonomic profiles
  kr_report <- as.data.table(kr_report)
  
  micro_taxids <- as.character(
    kr_report[grep("__Bacteria|__Fungi|__Viruses", mpa$V1), V7]
  )
  
  # Flag low-quality reads
  blast[, lowqual := as.integer(qcovs > min_qcovs)]
  
  # Compute quantiles once per ID
  # Step 1: Precompute quantiles per ID
  qstats <- blast[, .(
    q90_pident = quantile(pident, 0.9, na.rm = TRUE),
    q90_bitscore = quantile(bitscore, 0.9, na.rm = TRUE),
    q90_qcovs = quantile(qcovs, 0.9, na.rm = TRUE)
  ), by = id]
  
  # Step 2: Join back to main blast table
  blast <- merge(blast, qstats, by = "id", all.x = TRUE)
  
  # Keep only reads that meet high-quality thresholds
  blast <- blast[
    pident >= q90_pident | bitscore >= q90_bitscore | qcovs >= q90_qcovs
  ]
  
  # Flag reads as "other" if they are non-microbial and best hit
  blast[, other := as.integer(
    !(staxids %in% micro_taxids) &
      (pident == max(pident, na.rm = TRUE) | bitscore == max(bitscore, na.rm = TRUE))
  ), by = id]
  
  # Append "other" and "lowqual" reads to unmapped
  unmapped <- rbindlist(list(
    unmapped,
    blast[other == 1 & !(staxids %in% micro_taxids)][, .(id, read, staxids, sacc, pos, qcovs, pident, bitscore, type = "other")],
    blast[lowqual == 0][, .(id, read, staxids, sacc, pos, qcovs, pident, bitscore, type = "lowqual")]
  ), use.names = TRUE, fill = TRUE)
  
  # Keep only high-quality microbial reads
  blast <- blast[
    other == 0 & staxids %in% micro_taxids & lowqual == 1
  ]
  
  # Drop helper columns
  blast[, c("lowqual", "other", "q90_pident", "q90_bitscore", "q90_qcovs") := NULL]
  
  # ——————————————————————————————————————————————
  # STEP 4: Assign species-level taxonomic rank per read
  # ——————————————————————————————————————————————
  
  blast <- as.data.table(blast)
  
  # Step 4.1: Clean taxonomy names in kr_report for matching
  kr_tax_lookup <- kr_report[, .(
    taxid = as.character(V7),
    name_clean = str_squish(str_replace_all(V8, "[^[:alnum:]]", " ")),
    rank = str_to_lower(str_remove_all(V6, "[0-9]+"))
  )]
  
  # Step 4.2: Build a lookup table from taxid → name + rank + order
  rank_order <- data.table(rank = c('r','d','k','p','c','o','f','g','s'), order = 1:9)
  kr_tax_lookup <- merge(kr_tax_lookup, rank_order, by = "rank", all.x = TRUE)
  
  # Step 4.3: For each read ID, identify the best taxonomic hit based on:
  # - taxid frequency per read
  # - if multiple hits, choose one with highest rank specificity (highest 'order')
  
  # Count how many times each (read, taxid) pair appears
  tax_counts <- blast[, .N, by = .(id, staxids)]
  
  # Join with kraken taxonomy info
  tax_counts <- merge(tax_counts, kr_tax_lookup, by.x = "staxids", by.y = "taxid", all.x = TRUE)
  
  # Remove unknown or unmatched taxids
  tax_counts <- tax_counts[!is.na(order)]
  
  # Choose the most specific (highest order) name per read
  best_tax <- tax_counts[order(-order, -N), .SD[1], by = id]
  
  # Rename columns for clarity
  setnames(best_tax, c("rank", "name_clean"), c("rank", "tax_name"))
  
  # Step 4.4: Join best taxonomic label back to blast table
  blast <- merge(blast, best_tax[, .(id, rank, tax_name)], by = "id", all.x = TRUE)
  
  # ————————————————————————————————————
  # STEP 5: Return final list
  # ————————————————————————————————————
  
  list(
    micro = blast %>% dplyr::select(id, read, rank, staxids, tax_name, sacc, pos, qcovs, pident, bitscore),   # Filtered reads (microbial, good quality, species-level)
    unmapped = unmapped %>%  dplyr::select(id, read, staxids, type, sacc, pos, qcovs, pident, bitscore)   # Human, low-quality, or non-microbial reads  
  )
}


#' Identify taxa with uniquely mappable reads from BLAST alignments
#'
#' This function analyzes filtered BLAST alignment results to detect taxa (by taxid)
#' that have a sufficient number of **uniquely mapped reads**, i.e., reads mapping
#' exclusively to that taxon. It uses the results of \code{\link{prism_filter}} to
#' exclude non-microbial, low-quality, and human reads before performing the analysis.
#'
#' @param out_path String. Directory containing BLAST result CSV files (e.g., `*-blast.csv`).
#' @param sample String. Sample name, used for consistent file handling.
#' @param kr_report Data frame or data.table. Kraken report providing taxonomic hierarchy and IDs.
#' @param mpa Data frame. MetaPhlAn-like taxonomy table used to identify microbial taxa.
#' @param paired Logical. `TRUE` if paired-end reads were processed (informational only).
#' @param blast_file_pattern String. Pattern used to match BLAST CSV files in `out_path`.
#' @param nreads Integer. Minimum number of uniquely mappable reads required per taxon (default: 10).
#' @param min_qcovs Numeric. Minimum query coverage threshold passed to \code{prism_filter()}.

prism_multimapping = function(out_path, sample, kr_report, mpa, paired, blast_file_pattern, nreads=10, min_qcovs){
  
  #— 1) Filter out non-microbial reads from Blast results (human and model organism) ————————————————
  out = prism_filter(out_path, sample, kr_report, mpa, paired, blast_file_pattern, min_qcovs)
  blast = out$micro
  
  #— 2) Identify taxids with >=10 uniquely mappable reads ————————————————
  ta = table(blast$id, blast$staxids)
  ta[ta>0] = 1
  rs = rowSums(ta)
  ta = ta[rs == 1,,drop = F]
  tx = colSums(ta)
  
  uit = colnames(ta)[tx>=nreads]
  
  if (length(uit) == 0L && is.character(uit)) {
    stop("ERROR: no uniquely identifiable species found; aborting.")
  }
  
  uit
}

#' prism_genbank
#'
#' This function integrates filtered BLAST results with species-level taxonomic resolution
#' and maps them to GenBank gene/protein annotations based on read positions.
#'
#' @param blast input BLAST dataframe with mapping accession and position
#'
#' @return A data frame of BLAST-aligned reads annotated with species name, gene name, and product
#'         information from GenBank, along with taxonomic classification details.
#'         
prism_genbank = function(blast){
  
  genbank_dict = readRDS(paste0(prism_path, 'genbank/genbank_indices.RDS'))
  uacc = unique(blast$sacc)
  fn = genbank_dict$gbfile[unique(genbank_dict$list_indices[genbank_dict$accessions_flat %in% uacc])]
  
  message("Reading GenBank files ...")
  
  genbank = list()
  for(i in seq_along(fn)){
    # message(sprintf("📄 Reading GenBank file %d of %d ...", i, length(fn)))
    genbank[[i]] = readRDS(
      paste0(prism_path, 'genbank/rds/', fn[i], '.RDS')
    ) %>% subset(sacc %in% uacc)
  }
  
  genbank <- rbindlist(lapply(genbank, as.data.table), use.names = TRUE)
  
  message("Merging with GenBank gene/protein annotations ...")
  
  blast <- as.data.table(blast)
  blast[, `:=`(pos1 = pos, pos2 = pos)]
  genbank <- as.data.table(genbank)
  setkey(genbank, sacc, start, end)
  setkey(blast,    sacc, pos1,   pos2)
  
  res <- foverlaps(
    blast, genbank,
    by.x    = c("sacc", "pos1", "pos2"),
    by.y    = c("sacc", "start", "end"),
    nomatch = NA
  )
  
  res[, interval_width := fifelse(!is.na(start) & !is.na(end), end - start, Inf)]
  res_most_specific <- res[, .SD[which.min(interval_width)], by = .(id)]
  
  final <- res_most_specific[ , .(
    id, read, sacc, staxids, pos, qcovs, pident, bitscore,
    gene    = gene,
    product = product
  ) ] %>%
    droplevels() %>%
    arrange(id, pos) %>%
    left_join(
      kr_report %>%
        mutate(V7 = as.character(V7)) %>%
        select(V6, V7, V8) %>%
        rename(rank = V6, staxids = V7, tax_name = V8),
      by = 'staxids'
    ) %>%
    select(id, read, staxids, tax_name, rank, sacc, gene, product, everything())
  
  return(final %>% droplevels())
}


#' Compute multi-mapping ratios per taxon from BLAST alignments
#'
#' Given a BLAST alignment table, this function summarizes per-taxon multi-mapping
#' behavior. It identifies reads that map to a single taxid vs. multiple taxids,
#' computes the fraction of multi-mapped reads per taxid, and derives accession-level
#' overlap ratios between uniquely mapped and multi-mapped reads.
#'
#' @param blast Data frame or data.table. BLAST results with at least the columns:
#' \code{id} (read ID), \code{staxids} (subject taxid), and \code{sacc} (subject accession).

prism_multimapping_ratios = function(blast){
  
  # Create a contingency table of counts between `id` and `staxids`, convert to data frame,
  # filter to only rows with Freq > 0, group by `Var1` (i.e., `id`), and add a column `n` for the number of unique `staxids` per `id`
  t = table(blast$id, blast$staxids) %>% 
    data.frame() %>% 
    subset(Freq > 0) %>% 
    group_by(Var1) %>% 
    mutate(n = n())
  
  # Separate into two data frames:
  #  - `uniq`: reads (`id`) with only 1 associated `staxid`
  #  - `multi`: reads with multiple associated `staxid`s
  uniq = subset(t, n == 1)
  multi = subset(t, n > 1)
  
  # Create a summary data frame `d1`:
  #  - Count total number of occurrences per `staxid`
  #  - Count how many times each `staxid` appears in reads that match to multiple `staxid`s
  #  - Join these two tables on `staxid` (`Var1`)
  d1 = left_join(table(blast$staxids) %>% data.frame(), 
                 table(t$Var2) %>% data.frame(), 
                 by = 'Var1') 
  
  # Rename columns to more meaningful names
  colnames(d1) = c('staxids', 'total', 'multi')
  
  # Calculate the ratio of reads with multiple hits for each `staxid`
  # Add two new columns (`sacc_ratio4`, `sacc_ratio5`) initialized as NA
  d1 = d1 %>% 
    mutate(ratio = multi / total) %>% 
    dplyr::select(staxids, ratio) %>% 
    mutate(sacc_ratio4 = NA, sacc_ratio5 = NA)
  
  # Subset `blast` to get rows where `id` is uniquely mapped
  xuniq = blast %>% subset(id %in% uniq$Var1)
  
  # Subset `blast` to get rows where `id` maps to multiple `staxid`s
  xmulti = blast %>% subset(!(id %in% uniq$Var1))  # equivalent to id %in% uniq$Var1 == F
  
  # For each unique `staxid`, calculate:
  #  - sacc_ratio4: how many saccs are unique to unique matches (`xuniq`) vs shared with multi (`xmulti`)
  #  - sacc_ratio5: how many saccs are unique to multi matches (`xmulti`) vs shared with unique
  for (j in 1:nrow(d1)) {
    # Get unique accession numbers (sacc) for the current `staxid` in both unique and multi mappings
    s1 = unique(xuniq$sacc[xuniq$staxids == d1$staxids[j]])
    s2 = unique(xmulti$sacc[xmulti$staxids == d1$staxids[j]])
    
    # Calculate ratios of non-overlapping to overlapping saccs
    d1$sacc_ratio4[j] = length(setdiff(s1, s2)) / length(intersect(s1, s2))
    d1$sacc_ratio5[j] = length(setdiff(s2, s1)) / length(intersect(s1, s2))
  }
  
  d1 %>% mutate(staxids = as.character(staxids))
}

#' Summarize GenBank annotation diversity and abundance per taxon
#'
#' This function computes per-taxon summary statistics from GenBank-derived
#' annotation data, including counts and diversity indices for gene and product
#' annotations. It measures both abundance (number of entries) and diversity
#' (Shannon index) of annotated features per taxid.
#'
#' @param prod Data frame or tibble containing at least the columns:
#' \code{staxids} (taxonomic ID), \code{gene}, and \code{product}.

prism_genbank_stats = function(prod){
  prod %>% group_by(staxids) %>% 
    summarize(fprod = n(), fugene = length(unique(gene)), fuprod = length(unique(product)),
              prod_div = product %>% table() %>% vegan::diversity() %>% mean(),
              gene_div = gene %>% table() %>% vegan::diversity() %>% mean(),
              .groups = 'drop') %>% 
    mutate(fprod = fprod/sum(fprod), fugene = fugene/sum(fugene), fuprod = fuprod/sum(fuprod)) 
}


#' Quantify k-mer taxonomy signals and potential misclassification per taxon
#'
#' For a set of taxa present in BLAST results, this function parses an external
#' k-mer assignment report (`output.txt` in `out_path`), extracts lines for those
#' taxids, and aggregates k-mer counts per read across taxonomic ranks
#' (k, p, c, o, f, g, s) as well as special buckets (`unclassified`, `human`,
#' `unrelated`). It also summarizes lineage-level counts from a Kraken report to
#' derive simple misclassification/lineage proportion metrics. Outputs are
#' normalized to per-taxon proportions.
#'
#' @param out_path String. Directory containing the k-mer assignment `output.txt`
#'   file and where temporary files will be written.
#' @param sample String. Sample name used to name temporary files.
#' @param kr_report Data frame or data.table. Kraken report with taxonomy columns
#'   including rank (`V6`), taxid (`V7`), counts/reads (`V3`, `V5`), and names (`V8`).
#' @param blast Data frame or data.table. BLAST alignments containing at least a
#'   \code{staxids} column; unique values of \code{staxids} define the target taxa.

prism_misclass_kmertax = function(out_path, sample, kr_report, blast){
  
  # Get unique staxids from the `prod` object
  tx = unique(blast$staxids)
  
  # Find the output file that contains kmer info
  out = list.files(out_path, full.names = T) %>% str_subset('output.txt')
  
  # Construct awk command to extract lines matching the relevant taxids from the output file
  str = paste0("awk -F '\t' '$3 ~ /", paste0('\\(taxid ', tx, '\\)', collapse = '|'), "/' ", out, ' >> ', out_path, sample, '.tempout.txt')
  system(str)  # Run the awk command
  
  # Read the filtered output file into R
  out = read.delim(paste0(out_path, sample, '.tempout.txt'), header = F)
  
  # Clean and reformat output for further analysis
  microbiome_output_file = out %>% 
    dplyr::select(-V1) %>% 
    separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
    mutate(taxid = str_remove(taxid, '\\)') %>% trimws() %>% as.numeric(), name = trimws(name)) %>% 
    dplyr::rename(id = V2) %>%
    tibble()
  
  # Initialize empty lists to store outputs
  li = list()
  misclass = list()
  
  # Loop through each taxid in the list
  for(taxa in tx){
    
    # Extract lineage information for the current taxid from the `mpa` object
    lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
      str_extract(paste0('\\*', taxa, '\\*.*')) %>%
      str_remove('^\\*') %>% 
      str_remove('\\*$') %>% 
      str_split('\\*') %>% 
      unlist() %>%
      str_subset('NA', negate = T) %>% 
      as.numeric() %>%
      unique()
    
    # Get full lineage that ends in the current taxid
    full.lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*$')) %>% 
      str_remove('^\\*') %>% 
      str_remove('\\*$') %>% 
      str_split('\\*') %>% 
      unlist() %>% 
      as.numeric() 
    
    # Combine full lineage with any additional entries from `lin`
    full.lin = c(full.lin, lin) %>% unique()
    
    # Map lineage taxids to taxonomic ranks from `kr_report`
    ranks = kr_report$V6[kr_report$V7 %in% full.lin] %>% as.character() %>% str_to_lower() %>% str_replace('d', 'k') %>% str_remove('[0-9]')
    
    # Build a data frame of taxids and their ranks
    d = data.frame(taxid = full.lin %>% as.character(), rank = ranks %>% factor(level = unique(ranks)), stringsAsFactors = F)
    
    if('k' %in% d$rank == F){d = d %>% mutate(rank = str_replace(rank, 'r', 'k')) %>% distinct(rank, .keep_all = T)}
    
    # Count number of entries, total and unique reads from `kr_report`
    d2 = apply(d, 1, function(x){
      i = str_which(mpa$taxid, x[1])
      data.frame(n = length(i), r = sum(kr_report$V3[i]), u = sum(kr_report$V5[i]))
    }) %>% bind_rows()
    
    # Combine and filter by key taxonomic ranks
    d3 = cbind(taxid = taxa, d, d2) %>% 
      subset(rank %in% c('k','p','c','o','f','g','s')) %>% 
      group_by(taxid, rank) %>% 
      summarize_if(.predicate = is.numeric, .funs = sum)
    
    # Save lineage stats to `misclass` list
    misclass[[taxa]] = d3
    
    # Prepare dataframe for comparing against microbiome output
    d = data.frame(taxid = as.character(full.lin), rank = ranks, stringsAsFactors = F)
    
    # Subset microbiome output to current lineage
    out = subset(microbiome_output_file, taxid %in% lin) %>% 
      separate(V5, into = c('r1', 'r2'), sep = '\\|\\:\\|') 
    
    # Limit rows for performance reasons
    if(nrow(out) > 1000){out = out[1:1000, ]}
    
    # Initialize temporary list to store sample-level summaries
    df.list = list()
    
    # Loop through each read
    for(i in 1:nrow(out)){
      
      # Parse and summarize kmer hits from both r1 and r2 regions
      r = data.frame(pos = out[['r1']][i] %>% strsplit('\\s') %>% unlist(), stringsAsFactors = F) %>% 
        separate(pos, into = c('taxid', 'nkmer'), sep = ':', convert = T, fill = 'right')  %>% 
        group_by(taxid) %>% 
        summarize(n = sum(nkmer), .groups = 'drop') %>% 
        na.omit() %>% 
        rbind(data.frame(pos = out[['r2']][i] %>% strsplit('\\s') %>% unlist(), stringsAsFactors = F) %>% 
                separate(pos, into = c('taxid', 'nkmer'), sep = ':', convert = T, fill = 'right')  %>% 
                group_by(taxid) %>% 
                summarize(n = sum(nkmer), .groups = 'drop') %>% 
                na.omit()
        ) %>% 
        group_by(taxid) %>% 
        summarize(n = sum(n), .groups = 'drop') %>%
        mutate(taxid = as.character(taxid)) %>% 
        drop_na() %>% 
        left_join(d, by = 'taxid') 
      
      # Annotate special cases: unclassified, human, unrelated
      if(0 %in% r$taxid){r$rank[r$taxid == 0] = 'unclassified'}
      if(any(9606 %in% r$taxid)){r$rank[r$taxid %in% 9606] = 'human'}
      if(any(is.na(r$rank))){r$rank[which(is.na(r$rank))] = 'unrelated'}
      if(length(which(r$rank == 'unrelated')) > 1){
        r$n[r$rank == 'unrelated'] = sum(r$n[r$rank == 'unrelated'])
        r = r[-which(duplicated(r$rank)), ]
      }
      
      # Compile final summary for the read
      df.list[[i]] = data.frame(
        id = out$id[i],
        taxid = taxa,
        unclassified = ifelse('unclassified' %in% r$rank, r$n[r$rank == 'unclassified'], 0),
        host = ifelse('human' %in% r$rank, r$n[r$rank == 'human'], 0),
        unrelated = ifelse('unrelated' %in% r$rank, r$n[r$rank == 'unrelated'], 0),
        k = ifelse('k' %in% r$rank, r$n[r$rank == 'k'], 0),
        p = ifelse('p' %in% r$rank, r$n[r$rank == 'p'], 0),
        c = ifelse('c' %in% r$rank, r$n[r$rank == 'c'], 0),
        o = ifelse('o' %in% r$rank, r$n[r$rank == 'o'], 0),
        f = ifelse('f' %in% r$rank, r$n[r$rank == 'f'], 0),
        g = ifelse('g' %in% r$rank, r$n[r$rank == 'g'], 0),
        s = ifelse('s' %in% r$rank, r$n[r$rank == 's'], 0)
      )
    }
    
    # Aggregate all read-level summaries for the taxid
    li[[as.character(taxa)]] = bind_rows(df.list) %>% 
      group_by(taxid) %>% 
      summarize_if(.predicate = is.numeric, sum) %>% 
      bind_rows()
  }
  
  # Clean up temporary file
  system(paste0('rm ',  out_path, sample, '.tempout.txt'))
  
  # Combine all summaries into final outputs
  li = bind_rows(li) %>% arrange(taxid) %>% mutate(taxid = as.character(taxid)) %>% rename(staxids = taxid)
  li = li %>% pivot_longer(-c(staxids)) %>% group_by(staxids) %>% mutate(value = value/sum(value)) %>% pivot_wider(id_cols = c(staxids), names_from = name, values_from = value) %>% ungroup()      
  
  misclass = bind_rows(misclass) %>% pivot_longer(-c(taxid, rank)) %>% group_by(taxid, name) %>% 
    mutate(value = value/value[rank == 'k']) %>% subset(rank != 'k') %>% 
    pivot_wider(id_cols = c(taxid), names_from = c(rank, name), values_from = value) %>% ungroup() %>% 
    rename(staxids = taxid)
  
  # Return a list with two named elements: the taxonomy summary and misclassification summary
  list(ktax = li, misclass = misclass)
}


#' Extract species-level Kraken counts and fractions for BLAST-supported taxa
#'
#' Filters a Kraken report to species rank and returns per-taxon counts and
#' fractions only for taxa present in the BLAST results. Produces raw counts,
#' unique-read counts, and two normalized fractions: proportion of microbial
#' reads and proportion relative to a contamination panel.
#'
#' @param kr_report Data frame or data.table. Kraken report with at least columns:
#'   \code{V2} (count), \code{V3} (assigned reads), \code{V5} (unique reads),
#'   \code{V6} (rank), \code{V7} (taxid), and \code{V8} (name).
#' @param blast Data frame or data.table containing a \code{staxids} column; used
#'   to restrict output to taxa present in BLAST results.

prism_krcounts = function(kr_report, blast){
  kr_report %>% 
    mutate(
      V8 = str_squish(V8),
      species = V8,
      n2 = V2,
      n3 = V3,
      rank = V6,
      staxids = V7,
      uniq = V5,
      fmicro = V2 / sum(V2[V7 %in% c(2, 4751, 10239)]),
      fcontam = V2 / sum(V2[V7 %in% c(2100, 1747, 1282, 75775, 562, 1270, 1423, 573, 34062, 470)])
    ) %>% 
    subset(str_detect(rank, 'S')) %>% 
    dplyr::select(-rank) %>% 
    subset(staxids %in% blast$staxids) %>% 
    select(staxids, n2, n3, uniq, fmicro, fcontam) %>% 
    mutate(staxids = as.character(staxids))
}


#' Build per-taxon feature matrix and contamination scores, then summarize results
#'
#' This function aggregates multiple downstream features per taxon (`staxids`)—
#' GenBank annotation diversity, BLAST multi-mapping ratios, k-mer taxonomy
#' proportions, Kraken lineage proportions, and Kraken species counts—into a
#' single feature matrix. It then loads an XGBoost contamination model to score
#' each taxon and returns summarized outputs for reporting.
#'
#' @param out_path String. Directory for intermediate files used by k-mer parsing.
#' @param sample String. Sample identifier used for temporary filenames.
#' @param kr_report Data frame/data.table. Kraken report with at least columns
#'   \code{V2} (count), \code{V3} (assigned), \code{V5} (unique), \code{V6} (rank),
#'   \code{V7} (taxid), \code{V8} (name).
#' @param blast Data frame/data.table. Filtered BLAST alignments containing
#'   \code{id}, \code{staxids}, \code{sacc}, \code{pident}, \code{bitscore}, \code{qcovs}.
#' @param prod Data frame/tibble. GenBank-derived annotations with per-read
#'   \code{staxids}, \code{gene}, \code{product}, and (optionally) \code{tax_name}, \code{rank}.

prism_profile = function(out_path, sample, kr_report, blast, prod){
  # (1) Gene/product stats
  div = prism_genbank_stats(prod)
  
  # (2) Multimapping analysis 
  multimap = prism_multimapping_ratios(blast)
  
  # (3) kmer taxonomy 
  out = prism_misclass_kmertax(out_path, sample, kr_report, blast)
  ktax = out$ktax
  
  # (4) Kraken misclassifications
  misclass = out$misclass
  
  # (5) Kraken counts 
  krcounts = prism_krcounts(kr_report, blast)
  
  # Join them all by 'staxids' using full_join
  feature_mat = reduce(list(div, multimap, ktax, misclass, krcounts), full_join, by = 'staxids')
  
  # load xgboost contamination model and predict
  xg = readRDS(file = paste0(prism_path, 'prismxg.RDS'))
  pred = suppressMessages(predict(xg, feature_mat %>% dplyr::select(xg$feature_names) %>% mutate_all(~ replace(., is.infinite(.), NA)) %>% as.matrix()))
  
  # make species counts file 
  counts = prod %>% 
    subset(str_detect(rank, 'S')) %>% 
    group_by(tax_name, staxids) %>% 
    summarize(n = n(), .groups = 'drop') %>% 
    left_join(data.frame(staxids = feature_mat$staxids, pred = pred), by = 'staxids') %>% 
    arrange(-n) %>% 
    mutate(tax_name = tax_name %>% str_replace_all('_', ' '))
  
  # merge predictions with product data
  prod = prod %>% left_join(data.frame(staxids = feature_mat$staxids, pred = pred), by = 'staxids')
  
  list(counts = counts, prod = prod, xgmat = feature_mat)
}



#' Write PRISM-annotated FASTA(s) with BLAST-derived headers
#'
#' Creates PRISM microbiome FASTA file(s) by subsetting the input FASTA(s) to
#' reads present in the BLAST table and rewriting their headers to include
#' PRISM annotations (taxid and accession). Outputs final files as
#' `<sample>_1.fa` (and `<sample>_2.fa` if paired) in `out_path_final`.
#'
#' @param fa1 ShortRead object. Reads from R1.
#' @param fa2 ShortRead object (optional). Reads from R2 if paired-end.
#' @param blast Data frame/tibble with at least columns \code{id}, \code{staxids},
#'   \code{sacc}, and (optionally) \code{pos}. Used to select reads and build
#'   annotated headers.
#' @param out_path String. Working directory for intermediate files.
#' @param out_path_final String. Destination directory for final FASTA outputs.
#' @param sample String. Sample identifier used in output filenames.
#' @param paired Logical. `TRUE` if paired-end reads are provided and should be processed.

prism_FASTA = function(fa1, fa2, blast, out_path, out_path_final, sample, paired){
  
  # make PRISM microbiome fasta file(s)
  
  #--- read in FASTA #1 and extract headers ---
  headers = ShortRead::id(fa1) %>% as.character()
  
  #--- parse taxids from headers (in 1k‐read chunks if large) ---
  if(length(headers) > 1000){
    header_ids = list()
    i1 = seq(0, length(headers), by = 1000)[-1] 
    i2 = i1 - 999  
    if(max(i1) < length(fa1)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa1))}
    for(i in 1:length(i2)){
      # print(i)
      header_ids[[i]] = str_remove(headers[i2[i]:i1[i]], '\\s.*')
    }
    header_ids = unlist(header_ids)
  } else {header_ids = str_remove(headers, '\\s.*')}
  
  fa1_new = subset(fa1, header_ids %in% blast$id) # subset for reads that were blasted
  ShortRead::writeFasta(fa1_new, file = paste0(out_path, sample, '-new_1.fa'))
  headers = subset(headers, header_ids %in% blast$id)
  header_ids = subset(header_ids, header_ids %in% blast$id)
  df = data.frame(header = headers, id = header_ids) %>% left_join(blast %>% select(id, staxids, sacc, pos), by = 'id')
  new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc)
  write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
  str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
               out_path, sample, '-new_1.fa > ', out_path_final, sample, '_1.fa')
  system(str)
  str = paste0('rm ', out_path, sample, '_1.fa ', out_path, sample, '-new_1.fa ', out_path_final, 'new_headers.txt')
  system(str)
  
  if(paired == T){
    headers2 = ShortRead::id(fa2) %>% as.character()
    #--- parse taxids from headers (in 1k‐read chunks if large) ---
    if(length(headers2) > 1000){
      header_ids2 = list()
      i1 = seq(0, length(headers2), by = 1000)[-1] 
      i2 = i1 - 999  
      if(max(i1) < length(fa2)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa2))}
      for(i in 1:length(i2)){
        # print(i)
        header_ids2[[i]] = str_remove(headers2[i2[i]:i1[i]], '\\s.*')
      }
      header_ids2 = unlist(header_ids2)
    } else {header_ids2 = str_remove(headers2, '\\s.*')}
    
    fa2_new = subset(fa2, header_ids2 %in% blast$id) # subset for reads that were blasted
    ShortRead::writeFasta(fa2_new, file = paste0(out_path, sample, '-new_2.fa'))
    df = data.frame(header = headers2, id = header_ids2) %>% left_join(blast %>% select(id, staxids, sacc, pos), by = 'id')
    new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc)
    write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
    str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
                 out_path, sample, '-new_2.fa > ', out_path_final, sample, '_2.fa')
    system(str)
    str = paste0('rm ', out_path, sample, '_2.fa ', out_path, sample, '-new_2.fa ', out_path_final, 'new_headers.txt')
    system(str)
  }
  
}

