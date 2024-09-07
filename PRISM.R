suppressMessages(library(optparse))
suppressMessages(library(ShortRead))
suppressMessages(library(tidyverse))
suppressMessages(library(furrr))
suppressMessages(library(data.table))


option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--data_path"), action="store", help = "path to fastq files"),
  make_option(c("--kraken_path"), action="store", help = "path to Kraken2 executable"),
  make_option(c("--kraken_db_path"), action="store", help = "path to Kraken2 reference database"),
  make_option(c("--seqkit_path"), action="store", help = "path to SeqKit"),
  make_option(c("--star_path"), action="store", help = "path to STAR"),
  make_option(c("--star_genome_dir"), action="store", help = "path to STAR genome index (for --genomeDir parameter)"),
  make_option(c("--model_org_taxids"), action="store", default = 'NA', help = "path to .txt file with model organism taxids"),
  make_option(c("--blast_path"), action="store", help = "path to blastn"),
  make_option(c("--blast_db_path"), action="store", help = "path to blast reference database"),
  make_option(c("--prism_path"), action="store", help = "path to PRISM directory"),
  make_option(c("--min_read_per"), action="store", default = 10^4, help = "minimum reads per X to analyze. Default is 1 per 10,000 (min_read_per=10^4)"),
  make_option(c("--min_uniq_frac"), action="store", default = 5, help = "minimum ratio of unique k-mer to Kraken reads. Default=5"),
  make_option(c("--paired"), action="store", default = T, help = "paired end (T) or single end (F) reads"),
  make_option(c("--max_sample"), action="store", default = 1000, help = "max reads to sample"),
  make_option(c("--barcode_only"), action="store", default = F, help = "does read 1 only contain a barcode (T) or not (F)"),
  make_option(c("--fq1_end"), action="store", default = "_1.fastq", help = "fastq file name end format after sample name (e.g. sample_1.fastq, or sample_R1.fastq)"),
  make_option(c("--fq2_end"), action="store", default = "_2.fastq", help = "fastq file name end format after sample name (e.g. sample_1.fastq, or sample_R1.fastq)")
)
opt = parse_args(OptionParser(option_list = option_list))

sample = opt$sample
data_path = opt$data_path
kraken_path = opt$kraken_path
kraken_db_path = opt$kraken_db_path
seqkit_path = opt$seqkit_path
star_path = opt$star_path
star_genome_dir = opt$star_genome_dir
blast_path = opt$blast_path
blast_db_path = opt$blast_db_path
prism_path = ifelse(str_detect(opt$prism_path, '/$'), opt$prism_path, paste0(opt$prism_path, '/'))
model_org_taxids = ifelse(opt$model_org_taxids == 'NA', paste0(prism_path, opt$model_org_taxids), opt$model_org_taxids)
ranks = 'S'
min_read_per = opt$min_read_per
min_uniq_frac = opt$min_uniq_frac
paired = opt$paired %>% as.logical()
barcode_only = opt$barcode_only %>% as.logical()
max_sample = opt$max_sample
out_path = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/data/')
out_path_final = paste0(ifelse(str_detect(data_path, '/$'), data_path, paste0(data_path, '/')), sample, '_prism/')
fq1_end = opt$fq1_end
fq2_end = opt$fq2_end

if(dir.exists(out_path_final)){system(paste('rm -r', out_path_final))}
if(!dir.exists(out_path)){dir.create(out_path, recursive = T)}


cat(paste('\nSample name:', sample, '\n'))

cat(paste('Started Kraken2 at', Sys.time(), '\n'))
if(paired == T){
  # run Kraken
  str = paste0(kraken_path, ' \\\n',
               '--db ', kraken_db_path, ' \\\n',
               '--threads 24 \\\n',
               '--paired \\\n',
               '--use-names \\\n',
               '--report-minimizer-data \\\n',
               '--classified-out ', paste0(out_path, sample), '#.fq', ' \\\n',
               '--output ', paste0(out_path, sample, '.kraken.output.txt'), ' \\\n',
               '--report ', paste0(out_path, sample, '.kraken.report.txt'), ' \\\n',
               paste0(data_path, sample, fq1_end, ' \\\n'),
               paste0(data_path, sample, fq2_end, ' \n\n')
  )
} else {
  # run Kraken
  str = paste0(kraken_path, ' \\\n',
               '--db ', kraken_db_path, ' \\\n',
               '--threads 24 \\\n',
               '--use-names \\\n',
               '--report-minimizer-data \\\n',
               '--classified-out ', paste0(out_path, sample), '.fq', ' \\\n',
               '--output ', paste0(out_path, sample, '.kraken.output.txt'), ' \\\n',
               '--report ', paste0(out_path, sample, '.kraken.report.txt'), ' \\\n',
               paste0(data_path, sample, fq1_end, ' \n\n')
  )
}
# create MPA style report and standard Kraken report
str = paste0(str, 
             'cut -f1-3,6-8 ', 
             paste0(out_path, sample, '.kraken.report.txt'),
             ' > ',
             paste0(out_path, sample, '.kraken.report.std.txt'),
             '\n\n',
             '/projects/sd948/bassel/Kraken2Uniq/kraken2-master/kreport2mpa.py \\\n',
             '-r ', paste0(out_path, sample, '.kraken.report.std.txt'), ' \\\n',
             '-o ', paste0(out_path, sample, '.kraken.report.mpa.txt'), ' \\\n',
             '--intermediate-ranks',
             '\n\n'
) 

system(str)

###### extract kraken assigned microbiome reads and output 
cat(paste('Started extracting kraken microbial reads at', Sys.time(), '\n'))

# get reports and taxids
kr = paste0(out_path, sample, '.kraken.report.txt')
mpa = paste0(out_path, sample, '.kraken.report.mpa.txt')
kr = read.delim(kr, header = F) %>% mutate(V8 = trimws(V8))
kr = kr[-which(kr$V8 %in% c('root', 'unclassified')), ] %>% mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
mpa = read.delim(mpa, header = F)
mpa$taxid = NA
for(i in 2:nrow(mpa)){
  t_names = mpa[i,1] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    str_remove('.*__') %>% 
    str_replace_all('[^[:alnum:]]', '_') 
  mpa$taxid[i] = paste0('*', paste(kr$V7[match(t_names, kr$V8)], collapse = '*'), '*')
}
n = str_which(mpa$V1, 'k__Bacteria|k__Fungi|k__Viruses')
min_reads = ceiling(sum(kr$V2[kr$V7 %in% c(2,4751,10239)])/min_read_per)
min_reads = ifelse(min_reads > 9, min_reads, 10)
kr$uf = kr$V5/kr$V2
kr2 = kr # to be used later
kr = kr[n, ] %>% subset(str_remove_all(V6, '[0-9]') %in% ranks & (V2 > min_reads | (uf > min_uniq_frac & V2 > 9)))
taxid = kr$V7

# extract output
taxids_to_remove = c(0,1, setdiff(kr2$V7, kr2$V7[n]))
str = paste0("awk -F '\t' '$3 !~ /", paste0('\\(taxid ', taxids_to_remove, '\\)', collapse = '|'), "/' ", out_path, sample, '.kraken.output.txt >> ', out_path, sample, '.kraken.microbiome.output.txt')
system(str)
system(paste0('rm ', out_path, sample, '.kraken.output.txt'))
# extract ids only into another file
str = paste0("cut -d '\t' -f2 ", out_path, sample, '.kraken.microbiome.output.txt > ', out_path, sample, '.kraken.microbiome.ids.txt')
system(str)

if(paired == T){
  # read 1
  str = paste0(seqkit_path, ' grep --pattern-file ', out_path, sample, '.kraken.microbiome.ids.txt ', out_path, sample, '_1.fq', ' > ', out_path, sample, '_1.fq-temp')
  system(str)
  str = paste0(seqkit_path, ' fq2fa ', out_path, sample, '_1.fq-temp', ' > ', out_path, sample, '_1.fa')
  system(str)
  system(paste0('rm ', out_path, sample, '_1.fq-temp'))
  system(paste0('rm ', out_path, sample, '_1.fq'))
  # read 2
  str = paste0(seqkit_path, ' grep --pattern-file ', out_path, sample, '.kraken.microbiome.ids.txt ', out_path, sample, '_2.fq', ' > ', out_path, sample, '_2.fq-temp')
  system(str)
  str = paste0(seqkit_path, ' fq2fa ', out_path, sample, '_2.fq-temp', ' > ', out_path, sample, '_2.fa')
  system(str)
  system(paste0('rm ', out_path, sample, '_2.fq-temp'))
  system(paste0('rm ', out_path, sample, '_2.fq'))
} else {
  str = paste0(seqkit_path, ' grep --pattern-file ', out_path, sample, '.kraken.microbiome.ids.txt ', out_path, sample, '.fq', ' > ', out_path, sample, '_1.fq-temp')
  system(str)
  str = paste0(seqkit_path, ' fq2fa ', out_path, sample, '_1.fq-temp', ' > ', out_path, sample, '_1.fa')
  system(str)
  system(paste0('rm ', out_path, sample, '_1.fq-temp'))
  system(paste0('rm ', out_path, sample, '.fq'))
}
system(paste0('rm ', out_path, sample, '.kraken.microbiome.ids.txt'))

### subsample
cat(paste('Started subsampling at', Sys.time(), '\n'))
fq1 = readFasta(paste0(out_path, sample, '_1.fa'))
headers = ShortRead::id(fq1)

# get read taxids
if(length(headers) > 1000){
  header_taxid = list()
  i1 = seq(0, length(headers), by = 1000)[-1] 
  i2 = i1 - 999  
  if(max(i1) < length(fq1)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fq1))}
  for(i in 1:length(i2)){
    # print(i)
    header_taxid[[i]] = str_remove(headers[i2[i]:i1[i]], '.*taxid\\|')
  }
  header_taxid = unlist(header_taxid)
} else {header_taxid = str_remove(headers, '.*taxid\\|')}

# start subsampling
n = which(header_taxid %in% taxid)
micro_headers = headers[n]
micro_header_taxid = header_taxid[n]
n = sapply(taxid %>% as.character(), function(x) which(micro_header_taxid %in% x))
for(i in 1:length(n)){if(length(n[[i]]) > max_sample){n[[i]] = n[[i]][1:max_sample]}}
micro_headers = lapply(n, function(x) micro_headers[x])

# adjust subsampling using taxon phylogeny for undersampled reads if possible 
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
fa1 = subset(fq1, headers %in% micro_headers)
if(file.exists(paste0(out_path, sample, '_sub_fa1'))){system(paste('rm', paste0(out_path, sample, '_sub_fa1')))}
writeFasta(fa1, file = paste0(out_path, sample, '_sub_fa1'))

ids = ShortRead::id(fa1) %>% as.character() 
tax = ids %>% str_remove('.*taxid\\|') %>% as.numeric()
id_df = data.frame(taxid = tax, id = ids %>% str_remove('\\s.*'))

if(paired == T){
  # fq2
  fq2 = readFasta(paste0(out_path, sample, '_2.fa'))
  headers2 = ShortRead::id(fq2)
  
  if(length(headers2) > 1000){
    ids2 = list()
    i1 = seq(0, length(headers2), by = 1000)[-1] 
    i2 = i1 - 999  
    if(max(i1) < length(fq2)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fq2))}
    for(i in 1:length(i2)){
      # print(i)
      ids2[[i]] = str_remove(headers2[i2[i]:i1[i]], '\\s.*')
    }
    ids2 = unlist(ids2)
  } else {ids2 = str_remove(headers2, '\\s.*')}
  
  fa2 = subset(fq2, ids2 %in% micro_ids)
  if(file.exists(paste0(out_path, sample, '_sub_fa2'))){system(paste('rm', paste0(out_path, sample, '_sub_fa2')))}
  writeFasta(fa2, file = paste0(out_path, sample, '_sub_fa2'))
} 

##### star - remove human reads from subsampled fasta
cat(paste('Started STAR at', Sys.time(), '\n'))
if(paired == T){
  if(!dir.exists(paste0(out_path, 'star/'))){dir.create(paste0(out_path, 'star/'), recursive = T)}
  
  str = paste0(star_path, ' \\\n', 
               '--genomeDir ', star_genome_dir, ' \\\n',
               '--outFileNamePrefix ', out_path, 'star/', sample, '_sub_ \\\n', 
               '--readFilesIn ',  out_path, sample, '_sub_fa1 ', out_path, sample, '_sub_fa2 ',
               '\n\n'
  )
  system(str)
  
  sam = tryCatch(read.delim(paste0(out_path, 'star/', sample, '_sub_Aligned.out.sam'), comment.char = '@', header = F), error = function(e){NA})
  if(any(!is.na(sam))){
    fa_ids = ShortRead::id(fa1) %>% as.character() %>% str_remove('\\s.*')
    fa1 = subset(fa1, fa_ids %in% sam$V1 == F)
    fa_ids = ShortRead::id(fa2) %>% as.character() %>% str_remove('\\s.*')
    fa2 = subset(fa2, fa_ids %in% sam$V1 == F)
    writeFasta(fa1, file = paste0(out_path, sample, '_sub_fa1'))
    writeFasta(fa2, file = paste0(out_path, sample, '_sub_fa2'))  
    id_df$star = ifelse(id_df$id %in% sam$V1, 1, 0)
    write.table(id_df, file = paste0(out_path, sample, '_sub_ids.txt'))
  }
  
} else {
  
  if(!dir.exists(paste0(out_path, 'star/'))){dir.create(paste0(out_path, 'star/'), recursive = T)}
  str = paste0(star_path, ' \\\n', 
               '--genomeDir ', star_genome_dir, ' \\\n',
               '--outFileNamePrefix ', out_path, 'star/', sample, '_sub_ \\\n', 
               '--readFilesIn ',  out_path, sample, '_sub_fa1 ',
               '\n\n'
  )
  system(str)
  
  sam = tryCatch(read.delim(paste0(out_path, 'star/', sample, '_sub_Aligned.out.sam'), comment.char = '@', header = F), error = function(e){NA})
  if(!is.na(sam)){
    fa_ids = ShortRead::id(fa1) %>% as.character() %>% str_remove('\\s.*')
    fa1 = subset(fa1, fa_ids %in% sam$V1 == F)
    writeFasta(fa1, file = paste0(out_path, sample, '_sub_fa1'))
    id_df$star = ifelse(id_df$id %in% sam$V1, 1, 0)
    write.table(id_df, file = paste0(out_path, sample, '_sub_ids.txt'))
  }
}
if('star' %in% colnames(id_df) == F){id_df$star = 0; write.table(id_df, file = paste0(out_path, sample, '_sub_ids.txt'))}


# blast
cat(paste('Started species Blast at', Sys.time(), '\n'))

t = read.delim(model_org_taxids, header = F) %>% unlist() %>% as.character()
tx = table(id_df$taxid[id_df$star == 0])
tx = names(tx)[tx > 9]
tx = intersect(tx, kr2$V7[str_detect(kr2$V6, 'S')])

# read 1
if(barcode_only == F){
  str = paste0(
    'export BLASTDB=', blast_db_path, '; ', 
    blast_path, ' ', 
    '-max_hsps 1 ', 
    '-culling_limit 10 ', 
    '-max_target_seqs 10 ',
    '-query ',  out_path, sample, '_sub_fa1 ',
    '-db nt ',
    '-outfmt "6 delim=, qseqid sacc staxids sstart pident bitscore sstrand" ',
    '-out  ', paste0(out_path, sample, '_sub_fa1-blast.csv '), 
    '-taxids ', paste0(c(9606, t, tx), collapse = ',') 
  )
  system(str, ignore.stderr = T)    
}

# read 2
if(paired == T){
  str = paste0(
    'export BLASTDB=', blast_db_path, '; ', 
    blast_path, ' ', 
    '-max_hsps 1 ', 
    '-culling_limit 10 ', 
    '-max_target_seqs 10 ',
    '-query ',  out_path, sample, '_sub_fa2 ',
    '-db nt ',
    '-outfmt "6 delim=, qseqid sacc staxids sstart pident bitscore sstrand" ',
    '-out  ', paste0(out_path, sample, '_sub_fa2-blast.csv '), 
    '-taxids ', paste0(c(9606, t, tx), collapse = ',') 
  )
  system(str, ignore.stderr = T)  
}

####### Assess multimapping 
cat(paste('Started multimapping analysis', Sys.time(), '\n'))

# load blast data and remove human mappable reads
f = list.files(paste0(out_path), full.names = T) %>% str_subset('-blast.csv')
blast = c(); for(i in f){blast = rbind(blast, read.csv(i, header = F, stringsAsFactors = F))}
colnames(blast) =  c('id', 'sacc', 'staxids', 'pos', 'ppos', 'bitscore', 'strand')

blast = blast %>% 
  tibble() %>%
  group_by(id) %>% 
  mutate(human = ifelse(any(staxids == 9606), 1, 0)) %>% 
  subset(human == 0) %>% 
  dplyr::select(-human) %>% 
  mutate(staxids = str_remove(staxids, ';.*'))

id_df$blast = ifelse(id_df$id %in% blast$id , 0, 1)
id_df$blast = ifelse(id_df$star == 1, 0, id_df$blast)
write.table(id_df, file = paste0(out_path, sample, '_sub_ids.txt'))

# filter multimapping assignments - keep only "good" mapping percentage and bitscore
# remove reads whose best match is a model organism 
blast = blast %>% 
  group_by(id) %>% 
  filter(ppos >= quantile(ppos, 0.9) | bitscore >= quantile(bitscore, 0.9)) %>% 
  group_by(id) %>% 
  mutate(other = ifelse(staxids %in% kr2$V7 == F & (ppos == max(ppos) | bitscore == max(bitscore)), 1, 0)) %>% 
  subset(other == 0 & staxids %in% kr2$V7) %>%
  dplyr::select(-other) %>%
  distinct(id, staxids, .keep_all = T)

# get most specific taxonomic assignment for each read
future::plan(multisession)
rank_order = data.frame(rank = c('r','d','k','p','c','o','f','g','s'), order = 1:9)
x = table(blast$id, blast$staxids)
df = future_map(1:nrow(x), .progress = T, function(i){
  mpa[which(kr2$V7 %in% colnames(x)[which(x[i,]> 0)]), 'V1'] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    table() %>% 
    as.data.frame() %>%
    arrange(-Freq) %>% 
    filter(Freq == max(Freq)) %>% 
    mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>% 
    dplyr::select(-.) %>% 
    dplyr::select(-rank) %>% 
    left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = c('name')) %>% 
    mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>% 
    left_join(rank_order, by = 'rank') %>% 
    subset(order == max(order, na.rm = T)) %>% 
    tail(n=1) %>% 
    dplyr::select(rank, taxid, name) %>% 
    mutate(id = rownames(x)[i], name = str_replace_all(name, '_', ' '))
}) %>% bind_rows() %>% tibble()

# remove taxids with <10 counts or mean bitscore < 85
species_taxids = table(df$name) %>% 
  data.frame() %>% 
  dplyr::rename(name = Var1) %>% 
  subset(Freq > 9) %>% 
  left_join(df %>% dplyr::select(name, taxid, rank) %>% distinct(), by = 'name') %>% 
  subset(rank == 's') %>% 
  dplyr::select(taxid) %>% 
  unlist() %>% 
  unname()

multiblast_stats = blast %>% 
  subset(id %in% df$id & staxids %in% species_taxids) %>% 
  mutate(staxids = as.integer(staxids)) %>% 
  group_by(staxids) %>% 
  summarize(p = mean(ppos), b = mean(bitscore), n = n()) %>% 
  dplyr::rename(taxid = staxids) %>% 
  left_join(kr2 %>% 
              mutate(name = V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish(), 
                     taxid = V7 %>% as.integer()) %>% 
              dplyr::select(name, taxid), by = 'taxid') %>% 
  arrange(-n)

multiblast_data = list(read_data = df, initial_stats = multiblast_stats)

read_length = nchar(fa1@sread[1])
if(read_length < 51){min_bit = 60} else {min_bit = 85}

while(any(multiblast_stats$b < min_bit | multiblast_stats$n < 9)){
  blast = blast %>% 
    subset(staxids %in% multiblast_stats$taxid[multiblast_stats$b > min_bit]) %>% 
    subset(staxids %in% multiblast_stats$taxid[multiblast_stats$n > 9]) %>% 
    group_by(id) %>% 
    subset(staxids %in% kr2$V7)  
  
  x = table(blast$id, blast$staxids)
  
  df = future_map(1:nrow(x), .progress = T, function(i){
    mpa[which(kr2$V7 %in% colnames(x)[which(x[i,]> 0)]), 'V1'] %>% as.character() %>% 
      strsplit('\\|') %>% 
      unlist() %>% 
      table() %>% 
      as.data.frame() %>%
      arrange(-Freq) %>% 
      filter(Freq == max(Freq)) %>% 
      mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>% 
      dplyr::select(-.) %>% 
      dplyr::select(-rank) %>% 
      left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = 'name') %>% 
      mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>% 
      left_join(rank_order, by = 'rank') %>% 
      subset(order == max(order, na.rm = T)) %>% 
      tail(n=1) %>% 
      dplyr::select(rank, taxid, name) %>% 
      mutate(id = rownames(x)[i], name = str_replace_all(name, '_', ' '))
  }) %>% bind_rows() 
  
  species_taxids = table(df$name) %>% 
    data.frame() %>% 
    rename(name = Var1) %>% 
    subset(Freq > 9) %>% 
    left_join(df %>% dplyr::select(name, taxid, rank) %>% distinct(), by = 'name') %>% 
    subset(rank == 's') %>% 
    dplyr::select(taxid) %>% 
    unlist() %>% 
    unname()
  
  multiblast_stats = blast %>% 
    subset(id %in% df$id & staxids %in% species_taxids) %>% 
    mutate(staxids = as.integer(staxids)) %>% 
    group_by(staxids) %>% 
    summarize(p = mean(ppos), b = mean(bitscore), n = n()) %>% 
    rename(taxid = staxids) %>% 
    left_join(kr2 %>% 
                mutate(name = V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish(), 
                       taxid = V7 %>% as.integer()) %>% 
                dplyr::select(name, taxid), by = 'taxid') %>% 
    arrange(-n)
}

multiblast_data$final_stats = multiblast_stats

saveRDS(multiblast_data, file = paste0(out_path, sample, '_multiblast.RDS'))


# filter multimapping assignments - keep only "good" mapping percentage and bitscore
# remove reads whose best match is a model organism 
# keep only species identified above as having unique reads
f = list.files(paste0(out_path), full.names = T) %>% str_subset('-blast.csv')
blast = list()
for(i in 1:length(f)){
  x = read.csv(f[i], header = F, stringsAsFactors = F)
  colnames(x) = c('id', 'sacc', 'staxids', 'pos', 'ppos', 'bitscore', 'strand')
  x$staxids = as.character(x$staxids)
  blast[[i]] = tibble(read = i, x)
}
blast = bind_rows(blast)

blast = blast %>% 
  tibble() %>% 
  group_by(id) %>% 
  mutate(human = ifelse(any(staxids == 9606), 1, 0)) %>% 
  subset(human == 0) %>% 
  dplyr::select(-human) %>% 
  mutate(staxids = str_remove(staxids, ';.*'))

blast = blast %>% 
  subset(staxids %in% multiblast_data$final_stats$taxid) %>% 
  group_by(id) %>% 
  filter(ppos >= quantile(ppos, 0.9) | bitscore >= quantile(bitscore, 0.9)) %>% 
  group_by(id) %>% 
  mutate(other = ifelse(staxids %in% kr2$V7 == F & (ppos == max(ppos) | bitscore == max(bitscore)), 1, 0)) %>% 
  subset(other == 0 & staxids %in% kr2$V7) %>%
  dplyr::select(-other) %>%
  distinct(id, staxids, .keep_all = T)

#### get gene/product data 
genbank_dict = readRDS(paste0(prism_path, 'genbank/genbank_dict2.RDS')); genbank_dict$accession = as.character(genbank_dict$accession)
total_accessions = genbank_dict$accession %>% unlist() %>% paste() %>% strsplit(',') %>% unlist()


# reassess any remaining multimapping reads 
x = table(blast$id, blast$staxids)
df2 = future_map(1:nrow(x), .progress = T, function(i){
  mpa[which(kr2$V7 %in% colnames(x)[which(x[i,]> 0)]), 'V1'] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    table() %>% 
    as.data.frame() %>%
    arrange(-Freq) %>% 
    filter(Freq == max(Freq)) %>% 
    mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>% 
    dplyr::select(-.) %>% 
    dplyr::select(-rank) %>% 
    left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = 'name') %>% 
    mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>% 
    left_join(rank_order, by = 'rank') %>% 
    subset(order == max(order, na.rm = T)) %>% 
    tail(n=1) %>% 
    dplyr::select(rank, taxid, name) %>% 
    mutate(id = rownames(x)[i], name = str_replace_all(name, '_', ' '))
}) %>% bind_rows() 

multimapped = table(df2$name) %>% 
  data.frame() %>% 
  dplyr::rename(name = Var1) %>% 
  subset(Freq > 9) %>%
  left_join(df2 %>% dplyr::select(name, taxid, rank) %>% distinct(), by = 'name')

ranks_old = multimapped$rank
ranks_new = 1
multimapped_new = multimapped

while(any(multimapped_new$rank != 's') & all(ranks_old != ranks_new)){
  for(i in 1:nrow(multimapped)){
    if(multimapped_new$rank[i] == 's'){next}
    lin = str_subset(mpa$taxid, paste0('\\*', multimapped$taxid[i], '\\*')) %>% 
      str_extract(paste0('\\*', multimapped$taxid[i], '\\*.*')) %>%
      str_remove('^\\*') %>% 
      str_remove('\\*$') %>% 
      str_split('\\*') %>% 
      unlist() %>%
      str_subset('NA', negate = T) %>% 
      as.numeric() %>%
      unique() 
    lin = lin[-1]
    
    if(length(which(unique(multimapped_new$taxid) %in% lin)) == 1){
      idx = which(multimapped_new$taxid %in% lin)[1]
      multimapped_new$taxid[i] = multimapped_new$taxid[idx]
      multimapped_new$rank[i] = multimapped_new$rank[idx]
      multimapped_new$name[i] = multimapped_new$name[idx]
    } 
  }
  ranks_old = ranks_new
  ranks_new = multimapped_new$rank
}

multimap_df = df2 %>% left_join(left_join(multimapped, multimapped_new, by = 'Freq', suffix = c('', '_new')) %>% dplyr::select(-Freq), by = c('name', 'rank', 'taxid'))

blast = blast %>% 
  left_join(multimap_df, by = c('id')) %>% 
  subset(staxids == taxid_new)


## get products
cat(paste('Started microbial gene and product mapping at', Sys.time(), '\n'))

# future::plan(multisession)
uacc = unique(blast$sacc)
ii = future_map(1:length(uacc), .progress = T, function(x){str_which(genbank_dict$accession, uacc[x])}) %>% unlist() %>% unique()
genbank = list()
for(i in 1:length(ii)){
  # print(i)
  genbank[[i]] = readRDS(paste0(prism_path, 'genbank/rds/', genbank_dict$file[ii[i]]))
}
genbank = data.table::rbindlist(genbank) %>% tibble()
genbank = mutate_all(genbank, .funs = as.character)
genbank$start = as.numeric(genbank$start)
genbank$end = as.numeric(genbank$end)

read_length = nchar(fa1@sread[1])

dat = list()
sacc = unique(blast$sacc)
for(i in 1:length(sacc)){
  g2 = subset(genbank, accession == sacc[i])
  b2 = blast[which(blast$sacc == sacc[i]), ]
  dat2 = list()
  for(j in 1:nrow(b2)){
    # print(paste(i,j))
    g = subset(g2, b2$pos[j] >= start & b2$pos[j] <= end)
    if(nrow(g) == 0){g = subset(g2, b2$pos[j]+read_length >= start & b2$pos[j]+read_length <= end)}
    if(nrow(g) == 0){g[1,]=NA}
    g$id = b2$id[j]
    g$accession = b2$sacc[j]
    g = left_join(b2[j,], g %>% dplyr::rename(sacc = accession) %>% dplyr::select(-taxon), by = c('id', 'sacc'))
    dat2[[j]] = g
  }
  dat[[i]] = bind_rows(dat2) %>% tibble()
}

dat = bind_rows(dat) %>% droplevels()

dat$nbitscore = dat$bitscore/read_length

saveRDS(dat, file = paste0(out_path, sample, '-products.RDS'))
prod = dat; rm(dat)


#### KMER PHYLOGENY ANALYSIS and PHYLOGENY MISCLASSIFICATION ANALYSIS
cat(paste('Started kmer phylogeny analysis at', Sys.time(), '\n'))

tx = unique(prod$taxid)
out = list.files(out_path, full.names = T) %>% str_subset('output.txt')
str = paste0("awk -F '\t' '$3 ~ /", paste0('\\(taxid ', tx, '\\)', collapse = '|'), "/' ", out, ' >> ', out_path, sample, '.tempout.txt')
system(str)

out = read.delim(paste0(out_path, sample, '.tempout.txt'), header = F)
microbiome_output_file = out %>% 
  dplyr::select(-V1) %>% 
  separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
  mutate(taxid = str_remove(taxid, '\\)') %>% trimws() %>% as.numeric(), name = trimws(name)) %>% 
  dplyr::rename(id = V2) %>%
  tibble() 

li = list()
misclass = list()
for(taxa in tx){
  lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_extract(paste0('\\*', taxa, '\\*.*')) %>%
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>%
    str_subset('NA', negate = T) %>% 
    as.numeric() %>%
    unique()
  
  full.lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*$')) %>% 
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() 
  
  full.lin = c(full.lin, lin) %>% unique()
  
  ranks = kr2$V6[kr2$V7 %in% full.lin] %>% as.character() %>% str_to_lower() %>% str_replace('d', 'k') %>% str_remove('[0-9]')
  
  d = data.frame(taxid = full.lin %>% as.character(), rank = ranks %>% factor(level = unique(ranks)), stringsAsFactors = F)
  d2 = apply(d, 1, function(x){i = str_which(mpa$taxid, x[1]); data.frame(n = length(i), r = sum(kr2$V3[i]), u = sum(kr2$V5[i]))}) %>% bind_rows() 
  d3 = cbind(taxid = taxa, d, d2) %>% 
    subset(rank %in% c('k','p','c','o','f','g','s')) %>% 
    group_by(taxid, rank) %>% 
    summarize_if(.predicate = is.numeric, .funs = sum) 
  
  misclass[[taxa]] = d3
  
  d = data.frame(taxid = as.character(full.lin), rank = ranks, stringsAsFactors = F)
  
  out = subset(microbiome_output_file, taxid %in% lin) %>% separate(V5, into = c('r1', 'r2'), sep = '\\|\\:\\|') 
  
  if(nrow(out) > 1000){out = out[1:1000, ]}
  
  df.list = list()
  for(i in 1:nrow(out)){
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
    
    if(0 %in% r$taxid){r$rank[r$taxid == 0] = 'unclassified'}
    if(any(9606 %in% r$taxid)){r$rank[r$taxid %in% 9606] = 'human'}
    if(any(is.na(r$rank))){r$rank[which(is.na(r$rank))] = 'unrelated'}
    if(length(which(r$rank == 'unrelated')) > 1){r$n[r$rank == 'unrelated'] = sum(r$n[r$rank == 'unrelated']); r = r[-which(duplicated(r$rank)), ]}
    
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
  li[[as.character(taxa)]] = bind_rows(df.list) %>% group_by(taxid) %>% summarize_if(.predicate = is.numeric, sum) %>% bind_rows()
}

system(paste0('rm ',  out_path, sample, '.tempout.txt'))
li = bind_rows(li) %>% arrange(taxid) %>% mutate(taxid = as.character(taxid))

write.table(li, file = paste0(out_path, sample, '-kmerphylo.txt'))

misclass = bind_rows(misclass)

write.table(misclass, file = paste0(out_path, sample, '-misclass.txt'))


#### HUMAN SIMILARITY

cat(paste('Started human similarity analysis at', Sys.time(), '\n'))
system(paste0('mkdir ', out_path, 'human_similarity'))

if(paired == T){
  str = paste0(star_path, ' \\\n', 
               '--genomeDir ', star_genome_dir, ' \\\n',
               '--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \\\n',
               '--outFileNamePrefix ', out_path, 'human_similarity/', sample, '_ \\\n', 
               '--readFilesIn ',  out_path, sample, '_sub_fa1 ', out_path, sample, '_sub_fa2 '
  )
  system(str)
} else {
  str = paste0(star_path, ' \\\n', 
               '--genomeDir ', star_genome_dir, ' \\\n',
               '--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \\\n',
               '--outFileNamePrefix ', out_path, 'human_similarity/', sample, '_ \\\n', 
               '--readFilesIn ',  out_path, sample, '_sub_fa1 '
  )
  system(str)
}


ids = read.table(paste0(out_path, sample, '_sub_ids.txt'))
sam = read.delim(paste0(out_path, 'human_similarity/', sample, '_Aligned.out.sam'), comment.char = '@', header = F)[, c('V1','V3','V4','V5')]
colnames(sam) = c('id','seq','pos','qual')
sam = left_join(sam, ids %>% dplyr::select(id,taxid), by = 'id')
ids = subset(ids, star ==0 & blast==0)

d = left_join(table(sam$taxid) %>% data.frame(), table(ids$taxid) %>% data.frame(), by = 'Var1') %>% 
  left_join(table(distinct(sam, id, .keep_all = T)$taxid) %>% data.frame(), by = 'Var1') %>% 
  dplyr::rename(taxid = Var1) %>% 
  group_by(taxid) %>% 
  mutate(taxid = as.character(taxid)) %>% 
  summarize(hs_ratio1 = Freq.x/Freq.y, 
            hs_ratio2 = Freq/Freq.y,
            .groups = 'drop')

sam = sam %>% 
  subset(taxid %in% unique(prod$taxid)) %>% 
  mutate(taxid = as.character(taxid)) %>% 
  group_by(taxid) %>% 
  summarize(hs_n = n(), 
            hs_uniq = length(unique(seq)), 
            hs_div = ifelse(length(seq)>0, sapply(1:10, function(x) sample(seq, 50, replace = T) %>% table() %>% vegan::diversity()) %>% mean(), NA),
            hs_qual = mean(qual),
            .groups = 'drop'
  ) %>% 
  left_join(d, by = 'taxid') %>% 
  arrange(taxid %>% as.numeric())

write.table(sam, file = paste0(out_path, sample, '-humansimilarity.txt'))


## MULTIBLAST ANALYSIS

cat(paste('Started multiblast analysis at', Sys.time(), '\n'))

f = list.files(out_path, full.names = T) %>% str_subset('-blast.csv')

dat = list()
for(k in 1:length(f)){
  dat[[k]] = read.csv(f[k], header = F)
  colnames(dat[[k]]) = c('id', 'sacc', 'staxids', 'pos', 'ppos', 'bitscore', 'strand')
  dat[[k]]$staxids = as.character(dat[[k]]$staxids)
}
dat = bind_rows(dat)

x = dat
t = table(x$id, x$staxids) %>% data.frame() %>% subset(Freq>0) %>% group_by(Var1) %>% mutate(n=n())

uniq = subset(t, n == 1)
multi = subset(t, n > 1)

d1 = left_join(table(x$staxids) %>% data.frame(), table(t$Var2) %>% data.frame(), by = 'Var1') 
colnames(d1) = c('staxids', 'total','multi')
d1 = d1 %>% mutate(ratio = multi/total) %>% dplyr::select(staxids, ratio) %>% 
  mutate(sacc_ratio1 = NA, sacc_ratio2= NA, sacc_ratio3 = NA,sacc_ratio4 = NA,sacc_ratio5 = NA)

xuniq = x %>% subset(id %in% uniq$Var1)
xmulti = x %>% subset(id %in% uniq$Var1 == F) 

for(j in 1:nrow(d1)){
  s1 = unique(xuniq$sacc[xuniq$staxids == d1$staxids[j]])
  s2 = unique(xmulti$sacc[xmulti$staxids == d1$staxids[j]])
  d1$sacc_ratio1[j] = length(setdiff(s1,s2))/length(unique(c(s1,s2)))
  d1$sacc_ratio2[j] = length(setdiff(s2,s1))/length(unique(c(s1,s2)))
  d1$sacc_ratio3[j] = length(intersect(s1,s2))/length(unique(c(s1,s2)))
  d1$sacc_ratio4[j] = length(setdiff(s1,s2))/length(intersect(s1,s2))
  d1$sacc_ratio5[j] = length(setdiff(s2,s1))/length(intersect(s1,s2))
}

xuniq = xuniq %>% 
  group_by(staxids) %>% 
  summarize(mb_strand = abs(0.5-length(which(strand == 'plus'))/n()),
            mb_div = sapply(1:10, function(x) sample(sacc, 50, replace = T) %>% table() %>% vegan::diversity()) %>% mean(),
            mb_uniq = length(unique(sacc)),
            mb_ppos = mean(ppos),
            mb_bit = mean(bitscore), 
            type = 'uniq',
            .groups = 'drop') %>% 
  pivot_longer(-c(staxids, type))

xmulti = xmulti %>% 
  group_by(staxids) %>% 
  summarize(mb_strand = abs(0.5-length(which(strand == 'plus'))/n()),
            mb_div = ifelse(length(sacc)>0, sapply(1:10, function(x) sample(sacc, 50, replace = T) %>% table() %>% vegan::diversity()) %>% mean(), NA),
            mb_uniq = length(unique(sacc)),
            mb_ppos = mean(ppos),
            mb_bit = mean(bitscore), 
            type = 'multi',
            .groups = 'drop') %>% 
  pivot_longer(-c(staxids, type))

x2 = rbind(xuniq, xmulti) %>% 
  group_by(staxids, name) %>% 
  summarize(value = ifelse(any(type == 'multi'), value[type == 'uniq']/value[type == 'multi'], NA), .groups = 'drop') %>% 
  pivot_wider(id_cols = staxids, names_from = name, values_from = value)

x2 = left_join(d1, x2, by = c('staxids')) %>% dplyr::rename(taxid=staxids) %>% subset(taxid %in% unique(prod$taxid)) %>% 
  mutate_all(~ replace(., is.infinite(.), NA))

write.table(x2, file = paste0(out_path, sample, '-multiblast.txt'))


##################### Subset the microbiome fastq file and re blast
cat(paste('Started final Blast at', Sys.time(), '\n'))

tx = distinct(blast %>% dplyr::select(-taxid) %>% left_join(id_df %>% dplyr::select(id, taxid), by = 'id'), id, staxids, taxid) %>% 
  subset(staxids %in% multiblast_data$final_stats$taxid) %>% ungroup() %>% dplyr::select(taxid) %>% unique() %>% unlist() %>% as.character()

fa1_new = subset(fq1, header_taxid %in% tx)
headers = ShortRead::id(fa1_new)
# get read taxids
if(length(headers) > 1000){
  header_ids = list()
  i1 = seq(0, length(headers), by = 1000)[-1] 
  i2 = i1 - 999  
  if(max(i1) < length(fa1_new)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa1_new))}
  for(i in 1:length(i2)){
    # print(i)
    header_ids[[i]] = str_remove(headers[i2[i]:i1[i]], '\\s.*')
  }
  header_ids = unlist(header_ids)
} else {header_ids = str_remove(headers, '\\s.*')}
ShortRead::writeFasta(fa1_new, file = paste0(out_path, sample, '-new_1.fa'))

if(paired == T){    
  fa2_new = subset(fq2, header_taxid %in% tx)
  headers2 = ShortRead::id(fa2_new)
  # get read taxids
  if(length(headers2) > 1000){
    header_ids2 = list()
    i1 = seq(0, length(headers2), by = 1000)[-1] 
    i2 = i1 - 999  
    if(max(i1) < length(fa2_new)){i2 = c(i2, max(i1)+1); i1=c(i1, length(fa2_new))}
    for(i in 1:length(i2)){
      # print(i)
      header_ids2[[i]] = str_remove(headers2[i2[i]:i1[i]], '\\s.*')
    }
    header_ids2 = unlist(header_ids2)
  } else {header_ids2 = str_remove(headers2, '\\s.*')}
  
  ShortRead::writeFasta(fa2_new, file = paste0(out_path, sample, '-new_2.fa'))
}

# read 1
if(barcode_only == F){
  str = paste0(
    'export BLASTDB=', blast_db_path, '; ', 
    blast_path, ' ', 
    '-max_hsps 1 ', 
    '-culling_limit 5 ', 
    '-max_target_seqs 5 ',
    '-query ',  out_path, sample, '-new_1.fa ',
    '-db nt ',
    '-num_threads 12 ',
    '-outfmt "6 delim=, qseqid sacc staxids sstart pident bitscore sstrand qcovs" ',
    '-out  ', paste0(out_path, sample, '-final-blast.csv '), 
    '-taxids ', paste0(c(9606, multiblast_data$final_stats$taxid), collapse = ',') 
  )
  system(str, ignore.stderr = T)    
}

# read 2
if(paired == T){
  str = paste0(
    'export BLASTDB=', blast_db_path, '; ', 
    blast_path, ' ', 
    '-max_hsps 1 ', 
    '-culling_limit 5 ', 
    '-max_target_seqs 5 ',
    '-query ',  out_path, sample, '-new_2.fa ',
    '-db nt ',
    '-num_threads 12 ',
    '-outfmt "6 delim=, qseqid sacc staxids sstart pident bitscore sstrand qcovs" ',
    '-out  ', paste0(out_path, sample, '-final-blast2.csv '), 
    '-taxids ', paste0(c(9606, multiblast_data$final_stats$taxid), collapse = ',') 
  )
  system(str, ignore.stderr = T)  
}

############## multimapping and gene/protein mapping 

cat(paste('Started final multimapping and gene/protein mapping', Sys.time(), '\n'))

f = list.files(paste0(out_path), full.names = T) %>% str_subset('-final-blast')
blast2 = list()
for(i in 1:length(f)){
  x = read.csv(f[i], header = F, stringsAsFactors = F)
  colnames(x) = c('id', 'sacc', 'staxids', 'pos', 'ppos', 'bitscore', 'strand', 'qcovs')
  x$staxids = as.character(x$staxids)
  blast2[[i]] = tibble(read = i, x)
}
blast2 = bind_rows(blast2)

blast2 = blast2 %>% 
  tibble() %>% 
  group_by(id) %>% 
  mutate(human = ifelse(any(staxids == 9606), 1, 0)) %>% 
  subset(human == 0) %>% 
  dplyr::select(-human) %>% 
  mutate(staxids = str_remove(staxids, ';.*'))

rank_order = data.frame(rank = c('r','d','k','p','c','o','f','g','s'), order = 1:9)

x = table(blast2$id, blast2$staxids) %>% as.data.frame.matrix() 
x = x %>% mutate(duplicate_group = group_indices(x, across(everything())))
z = distinct(x)
future::plan(multisession)
df = future_map(1:nrow(z), .progress = T, function(i){
  mpa[which(kr2$V7 %in% colnames(z)[which(z[i,]> 0)]), 'V1'] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    table() %>% 
    as.data.frame() %>%
    arrange(-Freq) %>% 
    filter(Freq == max(Freq)) %>% 
    mutate(rank = str_remove(., '_.*'), name = str_remove(., '.*__') %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish()) %>% 
    dplyr::select(-.) %>% 
    dplyr::select(-rank) %>% 
    left_join(data.frame(name = kr2$V8 %>% str_replace_all('[^[:alnum:]]', ' ') %>% str_squish() , taxid = kr2$V7, rank = kr2$V6), by = c('name')) %>% 
    mutate(rank = str_to_lower(rank) %>% str_remove_all('[0-9]+')) %>% 
    left_join(rank_order, by = 'rank') %>% 
    subset(order == max(order, na.rm = T)) %>% 
    tail(n=1) %>% 
    dplyr::select(rank, taxid, name) %>% 
    mutate(duplicate_group = z$duplicate_group[i], name = str_replace_all(name, '_', ' '))
}) %>% bind_rows() %>% tibble()

multidf = x %>% rownames_to_column('id') %>% left_join(df, by = 'duplicate_group') %>% dplyr::select(id,rank,taxid,name)

write.table(multidf, file = paste0(out_path, sample, '-multidf.txt'))

# gene/protein mapping
future::plan(multisession)
uacc = unique(blast2$sacc) %>% as.character()
ii = future_map(1:length(uacc), .progress = T, function(x){str_which(genbank_dict$accession, uacc[x])}) %>% unlist() %>% unique()
genbank = list()
for(i in 1:length(ii)){
  # print(i)
  genbank[[i]] = readRDS(paste0(prism_path, 'genbank/rds/', genbank_dict$file[ii[i]]))
}
genbank = data.table::rbindlist(genbank) %>% tibble()
genbank = mutate_all(genbank, .funs = as.character)
genbank$start = as.numeric(genbank$start)
genbank$end = as.numeric(genbank$end)


x = blast2 %>% data.table::data.table()
x = x[, .(Group_List = .(list(.SD))), by = sacc]
s = x[[1]]; b = x[[2]]; names(b) = x$sacc
dt <- data.table::as.data.table(genbank)
dt_list <- dt[, .(Group_List = .(list(.SD))), by = accession]

bl = list()
for(i in 1:length(b)){
  # print(i)
  new_b = b[[i]] %>% data.frame()
  if(names(b)[i] %in% dt_list$accession == F){bl[[i]] = new_b; next}
  g = dt_list[[2]][which(dt_list$accession == names(b)[i])][[1]] %>% 
    data.frame() %>% distinct() %>% subset(end > start) %>% rownames_to_column('rn') %>% arrange(start)
  if(nrow(g) == 0){bl[[i]] = new_b; next}
  xx = g %>% 
    dplyr::select(-taxon,-locus) %>% 
    pivot_longer(-c(rn,definition,version,organism,gene,product,protein)) %>% 
    subset(!is.na(value))
  re = c(); for(j in 1:(nrow(xx)-1)){re[j] = ifelse(xx$value[j] > xx$value[j+1], 1, 0)}
  if(any(re == 1)){xx = xx[-which(re == 1), ]}
  int = findInterval(new_b$pos, xx$value, rightmost.closed = T)
  if(any(int == 0)){int[int == 0] = NA}
  if(any(names(table(xx$name[int])) == 'end')){int[xx$name[int] == 'end'] = NA}
  bl[[i]] = new_b %>% mutate(sacc = names(b)[i]) %>% cbind(g[xx$rn[int], ] %>% dplyr::select(-c(rn,taxon))) %>% tibble()
}


def = read.delim(paste0(prism_path, 'cog-20.def.tab'), header = F)
fun = read.delim(paste0(prism_path, 'fun-20.tab'), header = F)

bfinal = bind_rows(bl) %>% arrange(id, product, gene) %>% distinct(id, .keep_all = T) %>% 
  mutate(cog = def$V3[match(str_to_lower(gene), str_to_lower(def$V4))], 
         cat = def$V2[match(str_to_lower(gene), str_to_lower(def$V4))]) %>% 
  mutate(cat = fun$V3[match(cat, fun$V1)]) %>% 
  left_join(multidf %>% dplyr::select(id, rank, name) %>% dplyr::rename(tax_name = name), by = 'id') %>% 
  dplyr::select(id, staxids, rank, tax_name, everything()) %>% 
  dplyr::select(-organism, -locus)


## get taxonomy for each read
rank_order2 = data.frame(rank = c('R','R1','R2','R3','R4','R5','R6','R7','R8','R9','R10',
                                  'D','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10',
                                  'K','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10',
                                  'P','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10',
                                  'C','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',
                                  'O','O1','O2','O3','O4','O5','O6','O7','O8','O9','O10',
                                  'F','F1','F2','F3','F4','F5','F6','F7','F8','F9','F10',
                                  'G','G1','G2','G3','G4','G5','G6','G7','G8','G9','G10',
                                  'S','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'
),
order = 1:99)

tax = unique(bfinal$staxids)

tax_df1 = list()
for(i in 1:length(tax)){
  tax_df1[[i]] = data.frame(name = mpa[which(kr2$V7 == tax[i]),'V1'] %>% str_split('\\|') %>% unlist() %>% str_remove('.*__') %>% str_replace_all('_', ' '), 
                            taxid = mpa[which(kr2$V7 == tax[i]), 'taxid'] %>% str_remove('\\*') %>% strsplit('\\*') %>% unlist() %>% as.character()) %>% 
    left_join(kr2 %>% dplyr::select(V6,V7) %>% dplyr::rename(rank = V6, taxid = V7) %>% mutate(taxid = as.character(taxid)), by = 'taxid') %>% 
    pivot_longer(-rank) %>%
    mutate(name = paste0(rank,'_',name), staxids = tax[i]) 
}

tax_df1 = tax_df1 %>% bind_rows() %>% mutate(name = factor(name, levels = c(paste0(rank_order2$rank, '_taxid'), paste0(rank_order2$rank, '_name')))) %>% 
  distinct(name, value, staxids, .keep_all = T) %>% arrange(name) %>% 
  pivot_wider(id_cols = staxids, names_from = name, values_from = value)

tax_df2 = list()
for(i in 1:length(tax)){
  tax_df2[[i]] = data.frame(name = mpa[which(kr2$V7 == tax[i]),'V1'] %>% str_split('\\|') %>% unlist() %>% str_remove('.*__') %>% str_replace_all('_', ' '), 
                            taxid = mpa[which(kr2$V7 == tax[i]), 'taxid'] %>% str_remove('\\*') %>% strsplit('\\*') %>% unlist() %>% as.character()) %>% 
    left_join(kr2 %>% dplyr::select(V6,V7) %>% dplyr::rename(rank = V6, taxid = V7) %>% mutate(taxid = as.character(taxid)), by = 'taxid') %>% 
    mutate(staxids = tax[i])
}
tax_df2 = bind_rows(tax_df2)

saveRDS(list(tax_df1 = tax_df1, tax_df2 = tax_df2), file = paste0(out_path, sample, '-taxdf.RDS'))


### counts table and prism predictions 
ct = bfinal %>% group_by(staxids) %>% summarize(n = n()) %>% 
  left_join(tax_df2, by = 'staxids') %>% group_by(name, taxid, rank) %>% summarize(n = sum(n), .groups = 'drop') %>% 
  mutate(rank = factor(rank, levels = rank_order2$rank)) %>% arrange(rank)

ct = rbind(data.frame(name = 'Total', taxid = '00', rank = 'T', n = sum(kr$V2[kr$V7 %in% c(0,1)])),
           data.frame(name = 'Homo sapiens', taxid = '9606', rank = 'S', n = sum(kr$V2[kr$V7 %in% c(9606)])),
           ct
)

prod = readRDS(paste0(out_path, sample, '-products.RDS')) %>%   group_by(taxid, name) %>% 
  summarize(n = n(), fprod = n(), fugene = length(unique(gene)), fuprod = length(unique(product)),
            prod_div = sapply(1:10, function(x) sample(product, 100, replace = T) %>% table() %>% vegan::diversity()) %>% mean(),
            gene_div = sapply(1:10, function(x) sample(gene, 100, replace = T) %>% table() %>% vegan::diversity()) %>% mean(),
            .groups = 'drop') %>% 
  mutate(fprod = fprod/sum(fprod), fugene = fugene/sum(fugene), fuprod = fuprod/sum(fuprod))

counts = read.delim(paste0(out_path, sample, '.kraken.report.txt'), header =F) %>% tibble() %>% mutate(V8 = str_squish(V8)) %>% summarize(species = V8, n2=V2, n3=V3, rank = V6, taxid = V7, uniq = V5, fmicro = V2/sum(V2[V7 %in% c(2,4751,10239)]), fcontam = V2/sum(V2[V7 %in% c(2100, 1747, 1282, 75775, 562, 1270, 1423,573,34062,470)])) %>%subset(str_detect(rank, 'S')) %>% dplyr::select(-rank)
kphylo = read.table(paste0(out_path, sample, '-kmerphylo.txt')) %>% pivot_longer(-c(taxid)) %>% group_by(taxid) %>% mutate(value = value/sum(value)) %>% pivot_wider(id_cols = c(taxid), names_from = name, values_from = value) %>% ungroup()              
hs = read.table(paste0(out_path, sample, '-humansimilarity.txt'))               
mc = read.table(paste0(out_path, sample, '-misclass.txt')) %>% pivot_longer(-c(taxid, rank)) %>% group_by(taxid, name) %>% mutate(value = value/value[rank == 'k']) %>%subset(rank != 'k') %>% pivot_wider(id_cols = c(taxid), names_from = c(rank, name), values_from = value) %>% ungroup()
mb = read.table(paste0(out_path, sample, '-multiblast.txt'))  

test = mb %>% 
  left_join(counts, by = c('taxid')) %>% 
  left_join(kphylo, by = c('taxid')) %>% 
  left_join(hs, by = c('taxid')) %>% 
  left_join(prod, by = c('taxid')) %>% 
  left_join(mc, by = c('taxid')) %>% 
  dplyr::select(species, taxid, everything())

write.table(test, file = paste0(out_path, sample, '-xgmat.txt'))

xg = readRDS(file = paste0(prism_path, 'prismxg.RDS'))

pred = predict(xg, test %>% dplyr::select(xg$feature_names) %>% as.matrix())

ct = left_join(ct, data.frame(taxid = test$taxid %>% as.character(), pred = pred), by = 'taxid')

write.csv(ct, file = paste0(out_path_final, sample, '-counts.csv'), row.names = F, quote = F)

# add contamination score to results file
bfinal = left_join(bfinal %>% mutate(staxids = as.character(staxids)), data.frame(staxids = test$taxid %>% as.character(), pred = pred), by = 'staxids')
data.table::fwrite(bfinal, file = paste0(out_path_final, sample, '-results.csv'))

# make PRISM microbiome fasta file
headers <- headers %>% str_remove('kraken.*') %>% str_trim()
fa1_new = subset(fa1_new, header_ids %in% bfinal$id) # subset for reads that were blasted
ShortRead::writeFasta(fa1_new, file = paste0(out_path, sample, '-new_1.fa'))
headers = subset(headers, header_ids %in% bfinal$id)
header_ids = subset(header_ids, header_ids %in% bfinal$id)
df = data.frame(header = headers, id = header_ids) %>% left_join(bfinal %>% select(id, staxids, sacc, pos), by = 'id')
new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc, ' pos:', df$pos)
write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
             out_path, sample, '-new_1.fa > ', out_path_final, sample, '_1.fa')
system(str)
str = paste0('rm ', out_path, sample, '_1.fa ', out_path, sample, '-new_1.fa ', out_path_final, 'new_headers.txt')
system(str)

if(paired == T){
  headers2 <- headers2 %>% str_remove('kraken.*') %>% str_trim()
  fa2_new = subset(fa2_new, header_ids2 %in% bfinal$id) # subset for reads that were blasted
  ShortRead::writeFasta(fa2_new, file = paste0(out_path, sample, '-new_2.fa'))
  headers2 = subset(headers2, header_ids2 %in% bfinal$id)
  header_ids2 = subset(header_ids2, header_ids2 %in% bfinal$id)
  df = data.frame(header = headers2, id = header_ids2) %>% left_join(bfinal %>% select(id, staxids, sacc, pos), by = 'id')
  new_headers = paste0(">", df$header, ' | PRISM | staxids:', df$staxids, ' sacc:', df$sacc, ' pos:', df$pos)
  write(new_headers, file = paste0(out_path_final, 'new_headers.txt'))
  str = paste0("awk 'NR==FNR { h[++i] = $0; next } /^>/ { print h[++j]; next } { print }' ", out_path_final, 'new_headers.txt ', 
               out_path, sample, '-new_2.fa > ', out_path_final, sample, '_2.fa')
  system(str)
  str = paste0('rm ', out_path, sample, '_2.fa ', out_path, sample, '-new_2.fa ', out_path_final, 'new_headers.txt')
  system(str)
}


cat(paste('Finished at', Sys.time(), '\n\n'))





