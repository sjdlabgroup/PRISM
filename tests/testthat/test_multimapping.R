library(testthat)
source(file.path(dirname(testthat::test_path()), "../../functions.R"))
skip_if_not_installed("data.table"); skip_if_not_installed("dplyr"); skip_if_not_installed("stringr")
library(data.table); library(dplyr); library(stringr)

mock_kr <- function() data.frame(
  V2=c(150,120,60,10), V3=c(140,110,55,5), V5=c(70,50,30,3),
  V6=c("S","S","S","S"), V7=c(562,1270,9606,2157),
  V8=c("Escherichia_coli","Propionibacterium_acnes","Homo_sapiens","Archaea_sp"),
  stringsAsFactors=FALSE
)

mock_mpa <- function() data.frame(
  V1=c(
    "k__Bacteria|g__Escherichia|s__Escherichia_coli",
    "k__Bacteria|g__Propionibacterium|s__Propionibacterium_acnes",
    "k__Other|x__NonMicro|g__Foo|s__Bar"
  ),
  stringsAsFactors=FALSE
)

write_blast_csv_multi <- function(p1,p2){
  b1 <- data.table(
    id=c("u562_1","u562_2","u562_3","u1270_1","u1270_2","m1","m1","m2","m2","h1","o1","low1"),
    sacc=c("A","B","C","D","E","M1a","M1b","M2a","M2b","H","O","L"),
    staxids=c("562","562","562","1270","1270","562","1270","562","1270","9606","2157","562"),
    pos=c(10,20,30,10,20,10,15,10,15,999,300,200),
    qcovs=c(99,98,97,95,95,96,96,96,96,99,95,25),
    pident=c(99,99,99,98,98,98,98,98,98,99,96,55),
    bitscore=c(220,215,210,205,205,200,200,200,200,230,180,50)
  )
  b2 <- rbind(
    b1[1:5,][, `:=`(sacc=paste0(sacc,"x"), qcovs=qcovs-10, pident=pident-10, bitscore=bitscore-50)],
    data.table(id="h2", sacc="H2", staxids="9606", pos=500, qcovs=98, pident=99, bitscore=230)
  )
  fwrite(b1,p1,col.names=FALSE); fwrite(b2,p2,col.names=FALSE)
}

test_that("prism_multimapping: identifies ≥N uniquely mappable taxa", {
  out_dir <- tempdir()
  f1 <- file.path(out_dir, "Smm_1-blast.csv")
  f2 <- file.path(out_dir, "Smm_2-blast.csv")
  write_blast_csv_multi(f1,f2)

  uit <- suppressWarnings(prism_multimapping(
    out_path=out_dir, sample="Smm", kr_report=mock_kr(), mpa=mock_mpa(),
    paired=FALSE, blast_file_pattern="blast.csv", nreads=3, min_qcovs=50
  ))
  expect_true(is.character(uit) || is.factor(uit))
  expect_identical(as.character(uit), "562")

  uit2 <- suppressWarnings(prism_multimapping(
    out_path=out_dir, sample="Smm", kr_report=mock_kr(), mpa=mock_mpa(),
    paired=FALSE, blast_file_pattern="blast.csv", nreads=2, min_qcovs=50
  ))
  expect_setequal(as.character(uit2), c("562","1270"))
})

test_that("prism_multimapping: multi-mapped reads are excluded from unique counts", {
  out_dir <- tempdir()
  f1 <- file.path(out_dir, "Smm2_1-blast.csv")
  f2 <- file.path(out_dir, "Smm2_2-blast.csv")
  write_blast_csv_multi(f1,f2)

  uit <- suppressWarnings(prism_multimapping(
    out_path=out_dir, sample="Smm2", kr_report=mock_kr(), mpa=mock_mpa(),
    paired=FALSE, blast_file_pattern="blast.csv", nreads=3, min_qcovs=50
  ))
  expect_true("562" %in% as.character(uit))
})
