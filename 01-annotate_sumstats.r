#!/usr/bin/env Rscript

# Load relevant libraries
library(data.table)

# Load HRC file------
hrc <- fread("#/HRC_GRCh37.vcf", 
             stringsAsFactors = F, skip = 49, select = c("#CHROM", "POS", "ID", "REF", "ALT"))
setnames(hrc, c("#CHROM", "POS", "ID", "REF", "ALT"), c("CHR", "BP", "ID", "HRC_REF", "HRC_ALT"))

# Process each WBC file------
wbc_names <- c("eos", "bas", "neu", 
               "lym", "mon", "wbc_all")

# Create vector which will have all SNPs from all immune cells
snp_vec <- c()

for(wbc_name in wbc_names) {
  in_file <- paste0("#", wbc_name, ".txt")
  mydat <- fread(in_file, stringsAsFactors=F, 
               select = c("rs_number", "reference_allele", "other_allele", 
                          "eaf", "beta", "se", "p-value", "q_statistic", 
                          "q_p-value", "i2", "n_studies", "n_samples", "effects"))
  mydat <- mydat[(mydat$'i2' > 0.40 & mydat$'q_p-value' > 0.05) | (mydat$'i2' <= 0.40) | is.nan(mydat$'i2'),]
  mydat[, c("CHR", "BP") := tstrsplit(rs_number, ":", fixed=TRUE)]
  mydat[, c("BP", "minor_a", "major_a") := tstrsplit(BP, "_", fixed=TRUE)]
  mydat$BP <- as.integer(mydat$BP)
  mydat <- merge(mydat, hrc, by = c("CHR", "BP"))
  mydat[, ID := ifelse(ID == ".", rs_number, ID)]
  exp_file <- paste0("#", wbc_name, ".txt")
  out_file <- paste0("#", wbc_name, ".txt")
  mydat_exp <- mydat[mydat$'p-value' < 5e-8,]
  snp_vec <- c(snp_vec, mydat_exp$ID)
  fwrite(mydat_exp, exp_file, row.names = F, quote = F, sep = "\t", na = NA)
  fwrite(mydat, out_file, row.names = F, quote = F, sep = "\t", na = NA)
}

snp_vec <- unique(snp_vec)

for(wbc_name in wbc_names) {
  in_file <- paste0("#", wbc_name, ".txt")
  mydat <- fread(in_file, stringsAsFactors=F, 
               select = c("rs_number", "reference_allele", "other_allele", 
                          "eaf", "beta", "se", "p-value", "q_statistic", 
                          "q_p-value", "i2", "n_studies", "n_samples", "effects"))
  mydat <- mydat[(mydat$'i2' > 0.40 & mydat$'q_p-value' > 0.05) | (mydat$'i2' <= 0.40) | is.nan(mydat$'i2'),]
  mydat[, c("CHR", "BP") := tstrsplit(rs_number, ":", fixed=TRUE)]
  mydat[, c("BP", "minor_a", "major_a") := tstrsplit(BP, "_", fixed=TRUE)]
  mydat$BP <- as.integer(mydat$BP)
  mydat <- merge(mydat, hrc, by = c("CHR", "BP"))
  mydat[, ID := ifelse(ID == ".", rs_number, ID)]
  exp_file <- paste0("#", wbc_name, ".txt")
  out_file <- paste0("#", wbc_name, ".txt")
  mydat_exp <- mydat[mydat$ID %in% snp_vec,]
  mydat_exp[,phenotype := wbc_name]
  fwrite(mydat_exp, exp_file, row.names = F, quote = F, sep = "\t", na = NA)
  fwrite(mydat, out_file, row.names = F, quote = F, sep = "\t", na = NA)
}

# Process each CRC file------
crc_names <- c("colon", "distal", "female", "left", 
               "male", "overall", "proximal", "rectal")

# Create vector which will have all SNPs from all immune cells
snp_vec <- c()

for(crc_name in crc_names) {
  in_file <- paste0("#", crc_name, ".txt")
  mydat <- fread(in_file, stringsAsFactors=F, 
               select = c("Chromosome", "Position", "Allele1", "SNP",
                          "Allele2", "Freq1", "Effect", "StdErr", "P.value", 
                          "Direction", "HetISq", "HetPVal", "Neff", "REF", "ALT"))
  mydat <- mydat[(mydat$'HetISq' > 0.40 & mydat$'HetPVal' > 0.05) | (mydat$'HetISq' <= 0.40) | is.nan(mydat$'HetISq'),]
  mydat$Position <- as.integer(mydat$Position)
  mydat[, MarkerName := paste0(Chromosome, ":", Position, "_", Allele1, "_", Allele2)]
  exp_file <- paste0("#", crc_name, ".txt")
  out_file <- paste0("#", crc_name, ".txt")
  mydat[, SNP := ifelse(SNP == ".", MarkerName, SNP)]
  mydat_exp <- mydat[mydat$'P.value' < 5e-8,]
  snp_vec <- c(snp_vec, mydat_exp$SNP)
  fwrite(mydat_exp, exp_file, row.names = F, quote = F, sep = "\t", na = NA)
  fwrite(mydat, out_file, row.names = F, quote = F, sep = "\t", na = NA)
}

snp_vec <- unique(snp_vec)

for(crc_name in crc_names) {
  in_file <- paste0("#", crc_name, ".txt")
  mydat <- fread(in_file, stringsAsFactors=F, 
               select = c("Chromosome", "Position", "Allele1", "SNP",
                          "Allele2", "Freq1", "Effect", "StdErr", "P.value", 
                          "Direction", "HetISq", "HetPVal", "Neff", "REF", "ALT"))
  mydat <- mydat[(mydat$'HetISq' > 0.40 & mydat$'HetPVal' > 0.05) | (mydat$'HetISq' <= 0.40) | is.nan(mydat$'HetISq'),]
  mydat$Position <- as.integer(mydat$Position)
  mydat[, MarkerName := paste0(Chromosome, ":", Position, "_", Allele1, "_", Allele2)]
  exp_file <- paste0("#", crc_name, ".txt")
  out_file <- paste0("#", crc_name, ".txt")
  mydat[, SNP := ifelse(SNP == ".", MarkerName, SNP)]
  mydat_exp <- mydat[mydat$SNP %in% snp_vec,]
  mydat_exp[, phenotype := crc_name]
  fwrite(mydat_exp, exp_file, row.names = F, quote = F, sep = "\t", na = NA)
  fwrite(mydat, out_file, row.names = F, quote = F, sep = "\t", na = NA)
}

