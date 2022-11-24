library(optparse)
library(tidyverse)
library(dplyr)
library(stringr)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


source("RIO.R")

input <- function(infile) {
	 pfix = prefix()
   parameters <<- readParameters(infile)
   sample_name <<- parameters['sample_name', 2]
   fq <<- paste(pfix, parameters['fq', 2], sep="/")
   kraken_report <<- paste(pfix, parameters['kraken_report', 2], sep="/")
   mpa_report <<- paste(pfix, parameters['mpa_report', 2], sep="/")
   keep_original <<- parameters['keep_original', 2]
   ntaxid <<- as.numeric(parameters['ntaxid', 2])
}

run <- function() {}

output <- function(out_path) {

kr = read.delim(kraken_report, header = F)
kr = kr[-c(1:2), ]
mpa = read.delim(mpa_report, header = F)
n = str_which(mpa$V1, 'k__Bacteria|k__Fungi|k__Viruses')
taxid = kr$V7[n]
taxid.list = split(taxid, ceiling(seq_along(taxid)/ntaxid))

if(file.exists(paste0(out_path, sample_name, '_line_numbers.txt'))){
  system(paste0('rm ', out_path, sample_name, '_line_numbers.txt'))
}

for(i in 1:length(taxid.list)){
  print(paste('Finding reads', i, '/', length(taxid.list)))
  
  taxid = paste0("taxid|", taxid.list[[i]], collapse = "\\|") 
  taxid = paste0("'", taxid, "'")
  str = paste0("grep -wn ", taxid, " ", fq, " | grep -Eo '^[^:]+' >> ", out_path, sample_name, "_line_numbers.txt")
  system(str)
}

h = read.delim(paste0(out_path, sample_name, '_line_numbers.txt'), header = F)
r = h$V1+1
d = data.frame(h=h$V1,r=r) %>% rownames_to_column('n') %>% pivot_longer(-n)
write.table(d$value, file = paste0(out_path, sample_name, '_line_numbers.txt'), row.names = F, col.names = F)
print(paste0(out_path, sample_name, '_line_numbers.txt'))
print('Extracting reads')
str = paste0("awk 'NR==FNR{ a[$1]; next }FNR in a' ", out_path, sample_name, '_line_numbers.txt ', fq, " > ", out_path, sample_name, ".fa")
system(str)

str = paste0("sed -i 's/@/>/g' ", out_path, sample_name, '.fa')
system(str)

str = paste0(out_path, sample_name, '_line_numbers.txt')
system(paste('rm', str))

if(keep_original == F){
  system(paste('rm', fq))
}

print('Done')

}
