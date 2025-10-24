# End-to-end(ish) test: BLAST CSV (synthetic) -> prism_filter -> prism_multimapping
# No external tools are invoked; we only exercise in-R logic on realistic tables.

library(testthat)
source(file.path(dirname(testthat::test_path()), "../../functions.R"))

skip_if_not_installed("data.table")
skip_if_not_installed("dplyr")
skip_if_not_installed("stringr")

library(data.table)
library(dplyr)
library(stringr)

# ---- Minimal mocks for taxonomy inputs ------------------------------------------------
mock_kr <- function() {
  # Include bacterial (562), fungal-ish (4751), human (9606), and “other” (2157) taxids.
  data.frame(
    V2 = c(200, 150, 60, 10),         # counts
    V3 = c(180, 140, 55,  5),         # clade counts
    V5 = c(100,  70, 30,  3),         # unique k-mers
    V6 = c("S", "S", "S", "S"),       # all species rank for simplicity
    V7 = c(562, 4751, 9606, 2157),    # taxids
    V8 = c("Escherichia_coli", "Cutibacterium_acnes", "Homo_sapiens", "Archaea_sp"),
    stringsAsFactors = FALSE
  )
}

mock_mpa <- function() {
  data.frame(
    V1 = c(
      "k__Bacteria|p__Proteobacteria|g__Escherichia|s__Escherichia_coli",
      "k__Bacteria|p__Actinobacteria|g__Cutibacterium|s__Cutibacterium_acnes",
      "k__Other|x__NonMicro|g__Foo|s__Bar"
    ),
    stringsAsFactors = FALSE
  )
}

# ---- Synthetic BLAST CSV that yields all major groups (micro/human/other) -------------
write_blast_csv_pair <- function(out_dir, prefix = "S") {
  f1 <- file.path(out_dir, sprintf("%s_sub_fa1-blast.csv", prefix))
  f2 <- file.path(out_dir, sprintf("%s_sub_fa2-blast.csv", prefix))

  # Paired files are similar but not identical; together they cover categories.
  b1 <- data.table(
    id      = c(
      # uniquely mapped microbial 562 (≥3 reads)
      "u562_1","u562_1","u562_2","u562_2","u562_3","u562_3",
      # extra microbial (still 562)
      "u562_4",
      # human
      "h1",
      # other (non-microbial best hit)
      "o1","o1"
    ),
    sacc    = c("A1","A2","A3","A4","A5","A6","A7","H1","O1","O2"),
    staxids = c("562","562","562","562","562","562","562","9606","2157","2157"),
    pos     = c( 10,  20,  35,  50,  65,  80,  95, 400, 300, 305),
    qcovs   = c( 96,  95,  97,  95,  96,  95,  90,  98,  95,  95),
    pident  = c( 99,  98,  99,  98,  99,  98,  95,  99,  97,  96),
    bitscore= c(220, 210, 215, 205, 210, 205, 190, 230, 190, 185)
  )

  b2 <- data.table(
    id      = c("u562_1","u562_2","u562_3","h2","o2"),
    sacc    = c("A1x","A3x","A5x","H2","O3"),
    staxids = c("562", "562", "562", "9606", "2157"),
    pos     = c( 12,    36,    66,    420,   320),
    qcovs   = c( 90,    90,    90,     98,     95),
    pident  = c( 95,    95,    95,     99,     97),
    bitscore= c(195,   195,   195,     230,    185)
  )

  fwrite(b1, f1, col.names = FALSE)
  fwrite(b2, f2, col.names = FALSE)
  c(f1, f2)
}

# ---- Integration test ---------------------------------------------------------------
test_that("Integration: prism_filter -> prism_multimapping yields expected groups and UIT", {
  out_dir <- tempdir()
  files <- write_blast_csv_pair(out_dir, "Int")

  # Use a typical min_qcovs; prism_filter removes qcovs <= min_qcovs up front by design.
  res <- suppressWarnings(prism_filter(
    out_path = out_dir,
    sample   = "Int",
    kr_report= mock_kr(),
    mpa      = mock_mpa(),
    paired   = FALSE,
    blast_file_pattern = "blast.csv",
    min_qcovs = 50
  ))

  # Sanity on schemas
  expect_true(is.data.frame(res$micro))
  expect_true(is.data.frame(res$unmapped))
  expect_true(all(c("id","read","staxids","sacc","pos","qcovs","pident","bitscore") %in% names(res$micro)))
  expect_true(all(c("id","read","staxids","sacc","pos","qcovs","pident","bitscore") %in% names(res$unmapped)))

  # Ensure major groups have representation (by design of synthetic input)
  expect_true(any(res$micro$staxids   == "562"))     # microbial present
  expect_true(any(res$unmapped$staxids== "9606"))    # human present
  expect_true(any(res$unmapped$staxids== "2157"))    # other present

  # Feed into multimapping with nreads=3; should promote 562
  uit <- suppressWarnings(prism_multimapping(
    out_path = out_dir,
    sample   = "Int",
    kr_report= mock_kr(),
    mpa      = mock_mpa(),
    paired   = FALSE,
    blast_file_pattern = "blast.csv",
    nreads   = 3,
    min_qcovs= 50
  ))
  expect_true(is.character(uit) || is.factor(uit))
  expect_true("562" %in% as.character(uit))
})
