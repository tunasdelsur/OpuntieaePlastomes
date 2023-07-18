library(ape)
library(phangorn)
library(phyloch)
library(phyloU)
library(skimmingLociR)

setwd("F:/......./7_testing_topologies")
getwd() -> wd

read.csv("Gene_markers_orders.csv") -> gene.order

######################################
### Data
######################################

# t1 = full plastome
# t2 = genes

### Trees

setwd("F:/......../7_testing_topologies/Trees/TesteDelta1")

read.tree("T1_BT_SisterOpunti_rooted.tre") -> t1
read.tree("T2_MSA_SiterOpuntia_rooted.tre") -> t2

### Aligns

setwd("F:/......./7_testing_topologies/Alignments/1.Raw")

list.files(pattern=".fasta") -> files

sapply(files, read.dna, "fasta", simplify=F) -> aligns.raw
sub(".fasta", "", names(aligns.raw)) -> l0
sub("_re", "", l0) -> l0
paste(l0, "_raw", sep="") -> names(aligns.raw)

setwd("F:/......./7_testing_topologies/Alignments/2.Gblock-NoGaps")

list.files(pattern=".fasta") -> files

sapply(files, read.dna, "fasta", simplify=F) -> aligns.g
sub(".fasta", "", names(aligns.g)) -> l0
sub(".Gblock_NoGap", "", l0, fixed=T) -> l0
paste(l0, "_gblock", sep="") -> names(aligns.g)

setwd(wd)

c(aligns.raw, aligns.g) -> aligns
length(aligns) == length(aligns.raw)++length(aligns.g)

######################################
### sh.test
######################################

## SH-test
## H0 = 'The two trees being compared are congruent' 
## Delta ln > 0 support t1
## Delta ln < 0 suppport t2

names(aligns) -> labs
do.call(rbind, strsplit(labs, "_")) -> dat
colnames(dat) <- c("region", "type")
data.frame(dat, delta=NA, pvalue=NA) -> dat


for (i in 1:length(aligns)) {
  
  aligns[[i]] -> a1
  
  
  ### check tips
  
  checkPhylo(t1, a1) -> c1
  checkPhylo(t2, a1) -> c2
  c1$tree -> t1.t
  c2$tree -> t2.t
  
  phyDat(a1) -> a1
  
  ### Get delta - site-wise log-likeli (SH-test)
  
  fit1 <- pml(t1.t, a1)
  fit2 <- pml(t2.t, a1)
  fit1 <- optim.pml(fit1) # optimize edge weights
  fit2 <- optim.pml(fit2)
  
  SH.test(fit1, fit2, B=10000) -> sh.out
  sh.out[2,3] - sh.out[1,3] -> deltaln
  min(sh.out[,4]) -> pvalue
  dat$delta[i] <- deltaln
  dat$pvalue[i] <- pvalue
  cat("\r", i)
}

setwd(wd)


######################################
### Plot
######################################

cols <- c("gray70", "gray5")

dat[order(dat$delta, decreasing = T),] -> dat.o

hist(dat.o$delta)
hist(dat.o$pvalue)

dat.o$col <- NA
dat.o$col[dat.o$type == "raw"] <- cols[1]
dat.o$col[dat.o$type == "gblock"] <- cols[2]
unique(dat.o$col)

range(dat.o$delta)
ylim = c(-30, 30)

matrix(ncol=2, nrow=2) -> m

barplot(dat.o$delta, col = dat.o$col, ylim = ylim, xlab="All", border=dat.o$col)

######################################
### Export
######################################

### check p-values

dat.o[dat.o$pvalue < 0.05,] -> dat.s
dat.s[dat.s$delta > 0,] -> dat.s1
dat.s[dat.s$delta < 0,] -> dat.s2

nrow(dat.s)
nrow(dat.s1)
nrow(dat.s2)

hist(dat.s1$delta, breaks=100)
hist(dat.s2$delta, breaks=100)

### 

dat.o$support <- "Both"
dat.o$support[dat.o$delta > 1 ] <- "T1"
dat.o$support[dat.o$delta < -1] <- "T2"
dat.o$support[dat.o$pvalue > 0.05] <- "Both"

table(dat.o$support)

unique(dat.o$support)

write.csv(dat.o, "Gene_topology_results_delta_1.csv", row.names=F)

### include number of variables

dat.o$aligned_bp <- NA
dat.o$var_sites <- NA
dat.o$PIS <- NA

for (i in 1:nrow(dat.o)) {
  dat.o[i,] -> d0
  paste(d0$region, d0$type, sep="_") -> lab0
  match(lab0, names(aligns)) -> x
  aligns[[x]] -> a0
  alignStats(a0, include.amb = F) -> s0
  s0[2] -> d0$aligned_bp
  s0[3] -> d0$var_sites
  s0[4] -> d0$PIS
  d0 -> dat.o[i,]
  cat("\r", i)
}

### export

write.csv(dat.o, "Gene_topology_results_delta_1_with_alignstats.csv", row.names=F)


### table

as.factor(dat.o$support) -> dat.o$support

aggregate(dat.o$support, by=list(dat.o$type), FUN=table) -> summary
summary -> summary.all

cuts <- round(quantile(dat.o$PIS))
paste(cuts[-5], cuts[-1], sep=" < PIS =< ") -> labs

vector("list", length=length(labs)) -> sums
names(sums) <- labs

for (i in 2:length(cuts)) {
  cuts[i-1] -> start 
  cuts[i] -> end
  dat.o[which(dat.o$PIS > start ),] -> d0
  d0[which(d0$PIS <= end),] -> d0
  aggregate(d0$support, by=list(d0$type), FUN=table) -> sums[[i-1]]
}

do.call(rbind, sums) -> summary.comps

### igual 0

dat.o[which(dat.o$PIS < 1),] -> d0
aggregate(d0$support, by=list(d0$type), FUN=table) -> sum0

### merge

rbind(sum0, summary.comps, summary.all) -> sum.f
rownames(sum.f)[1:2] <- c("PIS = 0.1", "PIS = 0.2")
rownames(sum.f)[11:12] <- c("All.1", "All.2")

unlist(lapply(strsplit(rownames(sum.f), ".", fixed=T), "[", 1)) -> x
data.frame(PIS=x, type=sum.f$Group.1, sum.f$x) -> sum.f
sum.f[order(sum.f$type),] -> sum.f

write.csv(sum.f, "Gene_topology_summary_delta_1.csv", row.names = F)
sum.f


### plot

summary$x -> bpdat
rownames(bpdat) <- summary$Group.1

cols2 <- c("gray", "white", "black")

dat.o[dat.o$support != "Both",] -> dat.s
dat.s[dat.s$delta > 0,] -> dat.t1
dat.s[dat.s$delta < 0,] -> dat.t2
rbind(dat.t1,dat.t2) -> dat.t

### plastome order

which(is.na(match(dat.o$region, gene.order$region))) -> x
dat.o$region[x]

match(dat.o$region, gene.order$region) -> dat.o$order
dat.o[order(dat.o$order),] -> dat.ord


### Plot

pdf("Gene_topology_delta_1.pdf", width=11, height=8.5)
barplot(dat.o$delta, col = dat.o$col, ylim = ylim, xlab="Markers", border=dat.o$col)
legend("topright", legend=c("raw", "gblocked"), fill=cols, box.col = NA)
title("Delta 1")
barplot(dat.ord$delta, col = dat.ord$col, ylim = ylim, xlab="Markers - Plastome order", border=dat.ord$col)
legend("topright", legend=c("raw", "gblocked"), fill=cols, box.col = NA)
title("Delta 1")
barplot(dat.t$delta, col = dat.t$col, ylim = ylim, xlab="Markers (delta > 0.1)", border=dat.t$col)
legend("topright", legend=c("raw", "gblocked"), fill=cols, box.col = NA)
title("Delta 1")
barplot(bpdat, col = cols)
legend("topright", legend=c("raw", "gblocked"), fill=cols, box.col = NA)

barplot(t(bpdat), col=cols2)
plot(c(1,1), axes=F, col="white", xlab="", ylab="")
legend("center", legend=colnames(bpdat), fill=cols2)
dev.off()





