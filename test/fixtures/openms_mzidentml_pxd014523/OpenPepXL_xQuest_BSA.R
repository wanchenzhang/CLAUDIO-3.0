# Labeled Benchmark, xQuest vs. OpenPepXL

library("plyr")
library("dplyr")
library("Biostrings")

correct_protein_positions <- function(tab, fasta_db)
{
  # read fasta database
  fasta_database <- readAAStringSet(fasta_db)
  
  AbsPos1 <- character(length(tab$ID))
  AbsPos2 <- character(length(tab$ID))
  for (i in 1:length(tab$ID))
  {
    seq1 <- gsub("(Oxidation)", "", tab$sequence[i], fixed = T)
    seq1 <- gsub("(Carbamidomethyl)", "", seq1, fixed = T)
    match1 <- vmatchPattern(seq1, fasta_database)
    pos1 <- as.numeric(start(match1)[which(sapply(start(match1), length) > 0)])[1] + as.numeric(tab$xl_pos[i])
    AbsPos1[i] <- pos1
    
    # use the unmodified sequence from
    tab$sequence[i] <- seq1
    
    if (nchar(tab$beta_sequence[i]) > 0)
    {
      seq2 <- gsub("(Oxidation)", "", tab$beta_sequence[i], fixed = T)
      seq2 <- gsub("(Carbamidomethyl)", "", seq2, fixed = T)
      match2 <- vmatchPattern(seq2, fasta_database)
      pos2 <- as.numeric(start(match2)[which(sapply(start(match2), length) > 0)])[1] + as.numeric(tab$xl_pos2[i])
      AbsPos2[i] <- pos2
      tab$beta_sequence[i] <- seq2
    }
    else
    {
      AbsPos2[i] < "-"
    }
  }
  #tab <- cbind(tab, AbsPos1, AbsPos2)
  tab$AbsPos1 <- AbsPos1
  tab$AbsPos2 <- AbsPos2
  return(tab)
}

preprocess_raw_table <- function(raw, fasta_db, correct_positions) 
{
  protein1 <- raw$accessions
  protein2 <- raw$accessions_beta
  # set protein 2 for loop links
  protein2[which(raw$xl_type == "loop-link")] <- raw$accessions[which(raw$xl_type == "loop-link")]
  if (!is.null(raw$beta_sequence))
  {
    beta_sequence <- raw$beta_sequence
  }
  else
  {
    beta_sequence <- rep("-", times=length(raw$xl_type))
    for (i in 1:length(raw$xl_type))
    {
      if (raw$xl_type[i] == "cross-link")
      {
        # extract beta from ID: split by "-", take the only result vector [[1]], the second string from that
        beta_sequence[i] <- strsplit(raw$OpenXQuest.id[i], split="-", fixed=T)[[1]][2]
      }
    }
  }
  beta_sequence[which(raw$xl_type == "loop-link")] <- raw$sequence[which(raw$xl_type == "loop-link")]
  id <- paste(protein1, protein2, raw$XL_Protein_position_alpha, raw$XL_Protein_position_beta, sep="-")
  
  # make Chimera XlinkAnalyzer table
  XLA_table <- data.frame(ID = id, Protein1 = protein1, Protein2 = protein2, 
                          AbsPos1 = raw$XL_Protein_position_alpha, AbsPos2 = raw$XL_Protein_position_beta, 
                          score = raw$score, 
                          xl_type = raw$xl_type, sequence = raw$sequence, beta_sequence = beta_sequence, 
                          xl_pos = raw$xl_pos, xl_pos2 = raw$xl_pos2,
                          stringsAsFactors=F)
  
  if (correct_positions)
  {
    XLA_table <- correct_protein_positions(XLA_table, fasta_db)
  }
  
  # sort by position to identify unique links (switch all alpha and beta info, if pos2 < pos1)
  for (i in 1:length(XLA_table$ID))
  {
    if (as.numeric(XLA_table$AbsPos2[i]) < as.numeric(XLA_table$AbsPos1[i]))
    {
      AbsPos1 <- XLA_table$AbsPos2[i]
      AbsPos2 <- XLA_table$AbsPos1[i]
      Protein1 <- XLA_table$Protein2[i]
      Protein2 <- XLA_table$Protein1[i]
      
      XLA_table$AbsPos1[i] <- AbsPos1
      XLA_table$AbsPos2[i] <- AbsPos2
      XLA_table$Protein1[i] <- Protein1
      XLA_table$Protein2[i] <- Protein2
      
      pep2 <- XLA_table$sequence[i]
      pep1 <- XLA_table$beta_sequence[i]
      pos2 <- XLA_table$xl_pos[i]
      pos1 <- XLA_table$xl_pos2[i]
      
      XLA_table$sequence[i] <- pep1
      XLA_table$beta_sequence[i] <- pep2
      XLA_table$xl_pos[i] <- pos1
      XLA_table$xl_pos2[i] <- pos2
      
      XLA_table$ID[i] <- paste(Protein1, Protein2, AbsPos1, AbsPos2, sep="-")
    }
  }
  # renew ID column
  XLA_table$ID <- paste(XLA_table$Protein1, XLA_table$Protein2, XLA_table$AbsPos1, XLA_table$AbsPos2, sep="-")
  
  return(XLA_table)
}

add_topolink_lines <- function(XLA_table)
{
  # define conversion table between one- and three-letter codes for amino acids
  oneLetter <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  threeLetter <- c("ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR")

  oneLetterCodesAlpha <- substr(XLA_table$sequence, XLA_table$xl_pos+1, XLA_table$xl_pos+1)
  threeLetterCodesAlpha <- mapvalues(oneLetterCodesAlpha, oneLetter, threeLetter)
  
  oneLetterCodesBeta <- substr(XLA_table$beta_sequence, XLA_table$xl_pos2+1, XLA_table$xl_pos2+1)
  threeLetterCodesBeta <- mapvalues(oneLetterCodesBeta, oneLetter, threeLetter)
  
  XLA_table$OneLetterCodeAlpha <- oneLetterCodesAlpha
  XLA_table$ThreeCodeAlpha <- threeLetterCodesAlpha
  XLA_table$OneLetterCodeBeta <- oneLetterCodesBeta
  XLA_table$ThreeCodeBeta <- threeLetterCodesBeta
  
  topolink_lines <- paste0("  observed ", threeLetterCodesAlpha, " A ", XLA_table$AbsPos1, " ", threeLetterCodesBeta, " A ", XLA_table$AbsPos2)
  XLA_table$topolink <- topolink_lines
  return(XLA_table)
}


#input_dir <- "/home/eugen/MSData/2017_PaperBenchmark/Labeled_Bench_2019_03/xQuest/"
input_dir <- "/home/eugen/MSData/2017_PaperBenchmark/Labeled_Bench_2019_03/OpenPepXL/"
input_files <- list.files(input_dir)

#output_dir <- "/home/eugen/MSData/2017_PaperBenchmark/Labeled_Bench_2019_03/Structure/xQuest/"
output_dir <- "/home/eugen/MSData/2017_PaperBenchmark/Labeled_Bench_2019_03/Structure/OpenPepXL/"

fasta_db <- "/home/eugen/MSData/2017_PaperBenchmark/Labeled_Bench_2019_03/Structure/BSA_4F5S.fasta"
correct_positions <- TRUE

file <- input_files[2]

for (file in input_files)
{
  raw_file_path <- paste0(input_dir, file)
  raw_file <- read.table(raw_file_path,
                         header = T, sep = "\t", stringsAsFactors=F, comment.char = "@")
  raw_file <- raw_file[raw_file$xl_rank == 1,]
  XLA_table <- preprocess_raw_table(raw_file, fasta_db, correct_positions)
  XLA_table <- XLA_table[order(XLA_table$score, decreasing = T),]
  
  # filter out some unnecessary cases
  min_link_dist <- 3
  same_pep <- XLA_table$sequence == XLA_table$beta_sequence
  small_dist <- abs(as.numeric(XLA_table$AbsPos1) - as.numeric(XLA_table$AbsPos2)) <= min_link_dist
  throw_out <- same_pep & small_dist
  XLA_table <- XLA_table[!throw_out,]
  
  # count occurences of each unique link
  XLA_table$spectral_count <- sapply(XLA_table$ID, function(x) {sum(x == XLA_table$ID)})
  
  # only use links found in at least two spectra
  XLA_table <- XLA_table[which(XLA_table$spectral_count > 1),]
  
  # make hits unique
  XLA_table <- distinct(XLA_table, ID, .keep_all = T)
  
  if (XLA_table$score[1] < 10)
  {
    XLA_table$score <- XLA_table$score * 20
  }

  write.csv(XLA_table, file=paste0(output_dir, file))

  XLA_table <- add_topolink_lines(XLA_table)
  write.table(XLA_table$topolink, file=paste0(output_dir, file, "_topolink"), sep=" ", row.names = F, col.names = F, quote = F)
}



