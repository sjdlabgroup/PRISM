library(testthat)
source(file.path(dirname(testthat::test_path()), "../../functions.R"))
skip_if_not_installed("data.table"); skip_if_not_installed("dplyr"); skip_if_not_installed("stringr")
library(data.table); library(dplyr); library(stringr)

mock_kr <- function() data.frame(
  V2=c(120,80,50,10), V3=c(110,70,45,5), V5=c(60,30,25,4),
  V6=c("S","S","S","S"), V7=c(562,4751,9606,2157),
  V8=c("Escherichia_coli","Cutibacterium_acnes","Homo_sapiens","Archaea_sp"),
  stringsAsFactors=FALSE
)

mock_mpa <- function() data.frame(
  V1=c(
    "k__Bacteria|p__Proteo|g__Escherichia|s__Escherichia_coli",
    "k__Bacteria|p__Actino|g__Cutibacterium|s__Cutibacterium_acnes",
    "k__Other|x__NonMicro|g__Foo|s__Bar"
  ),
  stringsAsFactors=FALSE
)

write_blast_csv <- function(path){
  b <- data.table(
    id=c("r_micro1","r_micro1","r_micro2","r_micro2","r_h1","r_other1","r_other1","r_low1","r_low2"),
    sacc=c("A1","A2","A3","A4","H1","O1","O2","L1","L2"),
    staxids=c("562","562","562","562","9606","2157","2157","562","4751"),
    pos=c(100,120,80,95,999,300,305,90,110),
    qcovs=c(95,70,96,60,98,95,92,20,25),  # keep 'other' high across thresholds; ensure some low-qual
    pident=c(98,85,99,84,99,98,97,60,55),
    bitscore=c(200,150,210,140,220,200,195,50,40)
  )
  fwrite(b, path, col.names=FALSE)
}

test_that("prism_filter: schema and group presence are sane with representative inputs", {
  out_dir <- tempdir()
  f1 <- file.path(out_dir, "S1_sub_fa1-blast.csv"); f2 <- file.path(out_dir, "S1_sub_fa2-blast.csv")
  write_blast_csv(f1); write_blast_csv(f2)

  res <- suppressWarnings(prism_filter(
    out_path=out_dir, sample="S1", kr_report=mock_kr(), mpa=mock_mpa(),
    paired=FALSE, blast_file_pattern="blast.csv", min_qcovs=50
  ))

  expect_true(is.data.frame(res$micro)); expect_true(is.data.frame(res$unmapped))
  expect_true(all(c("id","read","staxids","sacc","pos","qcovs","pident","bitscore") %in% names(res$micro)))
  expect_true(all(c("id","read","staxids","sacc","pos","qcovs","pident","bitscore") %in% names(res$unmapped)))
  expect_true(any(res$micro$staxids=="562"))
  expect_true(any(res$unmapped$staxids=="9606"))
  expect_true(any(res$unmapped$staxids=="2157"))
})

test_that("prism_filter: threshold sensitivity yields fewer micro reads at stricter cutoff", {
  out_dir <- tempdir()
  f1 <- file.path(out_dir, "S2_sub_fa1-blast.csv"); write_blast_csv(f1)

  res_strict <- suppressWarnings(prism_filter(out_dir,"S2",mock_kr(),mock_mpa(),FALSE,"blast.csv",90))
  res_lenient<- suppressWarnings(prism_filter(out_dir,"S2",mock_kr(),mock_mpa(),FALSE,"blast.csv",50))
  expect_true(nrow(res_lenient$micro) >= nrow(res_strict$micro))
})

test_that("prism_filter: groups exist for human and other with this synthetic design", {
  out_dir <- tempdir()
  f1 <- file.path(out_dir, "S3_sub_fa1-blast.csv"); write_blast_csv(f1)

  res <- suppressWarnings(prism_filter(out_dir,"S3",mock_kr(),mock_mpa(),FALSE,"blast.csv",50))
  expect_true(any(res$unmapped$staxids=="9606"))
  expect_true(any(res$unmapped$staxids=="2157"))
  expect_true(any(res$micro$staxids=="562"))
})
