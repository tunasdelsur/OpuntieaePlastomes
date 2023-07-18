### install genbankr package ####
##
## if (!requireNamespace("BiocManager", quietly = TRUE))
##  install.packages("BiocManager")
##
## BiocManager::install("genbankr")
#######################################

library(genbankr)


#setwd

### Parse files

f0 = "/path_to_gb_file"
readLines(f0) -> gb0

### split
grep("FEATURES     ", gb0) -> s1
grep("ORIGIN     ", gb0) -> s2

gb0[1:s1] -> gb.begin
gb0[s2:length(gb0)] -> gb.end
gb0[(s1+1):(s2-1)] -> annot

### Find all genes

annot[grep("     gene      ", annot, fixed=T)] -> genes.pos
annot[grep("     gene      ", annot)+1] -> genes.names

### split positions
unlist(lapply(strsplit(genes.pos, " "), FUN=function(x)(x[[length(x)]]))) -> genes.pos
sub("complement(", "", genes.pos, fixed=T) -> genes.pos
gsub(")", "", genes.pos, fixed=T) -> genes.pos
sub("join(", "", genes.pos, fixed=T) -> genes.pos
gsub(",", "", genes.pos, fixed=T) -> genes.pos

strsplit(genes.pos, "..", fixed=T) -> genes.pos
lapply(genes.pos, as.numeric) -> genes.pos
lapply(genes.pos, FUN=function(x)(x[c(1,length(x))])) -> genes.pos

### get names

strsplit(genes.names, "\"", fixed=T) -> genes.names
unlist(lapply(genes.names, FUN=function(x)(x[[length(x)]]))) -> genes.names
lapply(strsplit(genes.names, " "), "[", 1) -> genes.names

names(genes.pos) <- genes.names
do.call(rbind, genes.pos) -> genes.pos
as.data.frame(genes.pos) -> genes.pos
colnames(genes.pos) <- c("start", "end")
genes.pos[order(genes.pos$start),] -> genes.pos

### create misc annotations

# model
#n1 <- "     misc_feature    1..200"
#n2 <- "                     /note=\"testing (zzz)\""
#n3 <- "                     /info=\"annotated by R\""
#n4 <- "                     /annotator=\"R\""

### Size of intergenic to annotate

min.annot.bp = 1

new.annot <- vector()

for (i in 2:nrow(genes.pos)) {
  
  genes.pos[(i-1):i,] -> pos0
  pos0$end[1]+1 -> start0
  pos0$start[2]-1 -> end0
  end0 - start0 -> l0
  if (l0 > min.annot.bp) {
    paste(start0, "..", end0, sep="") -> pos1
    paste(paste(rownames(pos0), collapse="-"), "spacer") -> name0
    
    n1 <- paste("     misc_feature    ", pos1, sep="")
    n2 <- paste("                     /note=\"", name0, "\"", sep = "")
    n3 <- "                     /info=\"annotated by R\""
    n4 <- "                     /annotator=\"R\""
    
    c(n1,n2,n3,n4) -> annot0
    c(new.annot,annot0) -> new.annot

  }
}

### export

c(gb.begin, annot, new.annot, gb.end) -> gb.new
writeLines(gb.new, con=paste(sub(".gb", "", f0, fixed=T), "_annotated_IntergenicSpacer.gb", sep=""))


### rps12 must be manually adjustes
