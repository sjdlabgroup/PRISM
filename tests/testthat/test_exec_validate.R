library(testthat)
source(file.path(dirname(testthat::test_path()), "../../functions.R"))


test_that("prism_kraken validates missing executables/db and inputs", {
  expect_error(
    prism_kraken(
      sample="S1", data_path=tempdir(),
      kraken_path=tempfile(),           # nonexistent
      kraken_db_path=tempfile(),        # nonexistent
      seqkit_path=tempfile(),           # nonexistent
      out_path=tempdir(),
      min_read_per=1e4, min_uniq_frac=5,
      paired=FALSE, fq1_end="_1.fastq", fq2_end="_2.fastq", threads=1
    ),
    regexp="kraken2 executable|Kraken2 DB|SeqKit executable"
  )
})

test_that("prism_minimap validates missing executable/index/inputs", {
  expect_error(
    prism_minimap(
      sample="S1", out_path=tempdir(),
      minimap2_path=tempfile(), minimap2_index=tempfile(),
      threads=1, paired=FALSE
    ),
    regexp="minimap2 executable|index|Extracted FASTA #1"
  )
})

test_that("prism_star validates missing STAR and genomeDir", {
  out <- tempdir()
  f1  <- file.path(out,"S1_1.fa"); writeLines(c(">r1","ACGT"), f1)
  expect_error(
    prism_star("S1", out, star_path=tempfile(), star_genome_dir=tempfile(),
               paired=FALSE, id_df=data.frame(id="r1"),
               fa1=ShortRead::readFasta(f1), fa2=NULL),
    regexp="STAR not found|genomeDir not found"
  )
})

