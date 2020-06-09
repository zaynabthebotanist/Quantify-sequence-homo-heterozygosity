
library(Biostrings)
library(ips)
library(stringr)
library(strataG)
library(adegenet)
library(dplyr)

list.files()
filenames <- as.list(list.files(pattern = "*.fasta"))
n <- c(1:length(filenames))

#the sequences as DNAbin
l <- list()
for (i in n) {
  l[[i]] <- as.DNAbin(as.matrix(readDNAStringSet(filenames[[i]])))
}

#number of pis
vars <- list()
for (i in n) {
  vars[[i]] <- pis(l[[i]], what = "abs", use.ambiguities = FALSE)
}

#positions of pis 
index <- list()
for (i in n) {
  index[[i]] <- pis(l[[i]], what = "ind", use.ambiguities = FALSE)
}

#the sequences as DNAStringSet
f <- list()
for (i in n) {
  f[[i]] <- readDNAStringSet(filenames[[i]], use.names = TRUE)
}

##as gtype
gt <- list()
for (i in n) {
  gt[[i]] <- genind2gtypes(DNAbin2genind(l[[i]]))
}

##look at the gtype variants
df <- list()
for (i in n) {
  df[[i]] <- data.frame(gt[[i]]@data)
}  

dfa <- list()
for (i in n) {
  dfa[[i]] <- data.frame(gsub( ".1", "", as.character(df[[i]][,3])))
}

for (i in n) {
  df[[i]][,3] <- dfa[[i]]
}

##Filter using the indexing list (index) above
n <- c(1:length(gt))
fin <- list()
for (i in n) {
  fin[[i]] <- subset(df[[i]], df[[i]][,3] %in% index[[i]])
}

##split by index
fina <- list()
for (i in n) {
  fina[[i]] <- split(fin[[i]], fin[[i]][,3])
}

#remove empty dfs in the list
final <- list()
for (i in n) {
  final[[i]] <- fina[[i]][sapply(fina[[i]], function(x) dim(x)[1]) > 0]
}

finale <- final
for (i in n) {
  finale[[i]] <- bind_rows(finale[[i]])
}

##number of rows in each
dims <- list()
for (i in n) {
  dims[[i]] <- dim(finale[[i]])[1]
}

names <- list()
for (i in n) {
  names[[i]] <- data.frame(rep(filenames[[i]], dims[[i]]))
}

##merge the scaffold names and the df of alleles
mg <- list() 
for (i in n) {
  mg[[i]] <- cbind(finale[[i]], names[[i]])
}

mgsplit <- bind_rows(mg)
mgsplit <- split(mgsplit, list(mgsplit[,3], mgsplit[,5]))

mgs <- mgsplit[sapply(mgsplit, function(x) dim(x)[1]) > 0]

n <- c(1:length(mgs))
tab <- list()
for (i in n) {
  tab[[i]] <- data.frame(table(mgs[[i]]$allele))
  tab[[i]] <- tab[[i]][order(-tab[[i]]$Freq),]
}

##calculate homo- and heterozygosity
##some loci have more than one alternate allele

ref.allele <- list()
for (i in n) {
  ref.allele[[i]] <- tab[[i]][1,1]
}

alt.allele1 <- list()
for (i in n) {
  alt.allele1[[i]] <- tab[[i]][2,1]
}

alt.allele2 <- list()
for (i in n) {
  if (nrow(tab[[i]]) > 2)
  alt.allele2[[i]] <- tab[[i]][[3,1]]
}


##now a series loops to figure out the states 0, 1, 2

# haplotype1
mgshap1 <- list()
for (i in n) {
  mgshap1[[i]] <- mgs[[i]][endsWith(mgs[[i]][,1], "_hap1") == TRUE,]
  mgshap1[[i]][,1] <- gsub("_hap1", "", mgshap1[[i]][,1])
}

#haplotype 2
mgshap2 <- list()
for (i in n) {
  mgshap2[[i]] <- mgs[[i]][endsWith(mgs[[i]][,1], "_hap1") == FALSE,]
  mgshap2[[i]][,1] <- gsub("_hap2", "", mgshap2[[i]][,1])
}

library(tidyverse)

merg <- list()
for (i in n) {
  merg[[i]] <- full_join(mgshap1[[i]], mgshap2[[i]], by="id")
}

##make allele states
for (i in n) {
  merg[[i]][,10] <- ref.allele[[i]] 
  merg[[i]][,11] <- alt.allele1[[i]]
}

for (i in n) {
  merg[[i]][,12] <- merg[[i]]$allele.x == merg[[i]]$V10
  merg[[i]][,13] <- merg[[i]]$allele.y == merg[[i]]$V10
}

for (i in n) {
  merg[[i]][,14] <- gsub("TRUE", "0", merg[[i]][,12])
  merg[[i]][,14] <- gsub("FALSE", "1", merg[[i]][,14])
  merg[[i]][,15] <- gsub("TRUE", "0", merg[[i]][,13])
  merg[[i]][,15] <- gsub("FALSE", "1", merg[[i]][,15])
}

only <- list()
for (i in n) {
  only[[i]] <- merg[[i]][,14:15]
  only[[i]][,1] <- as.numeric(only[[i]][,1])
  only[[i]][,2] <- as.numeric(only[[i]][,2])
  only[[i]][,3] <- rowSums(only[[i]][,1:2])
}

ave <- list()
for (i in n) {
  ave[[i]] <- table(only[[i]][,3])
}

for (i in n) {
  if (sum(tab[[i]][,2]) > 72)
    ave[[i]] = NA
}

##count 0, 1 and 2

zeros <- list()
for (i in n) {
  zeros[[i]] <- ave[[i]][1]/sum(ave[[i]])
}

ones <- list()
for (i in n) {
  ones[[i]] <- ave[[i]][2]/sum(ave[[i]])
}

twos <- list()
for (i in n) {
  twos[[i]] <- ave[[i]][3]/sum(ave[[i]])
}

zero <- (sum(na.omit(unlist(zeros))))/length(tab)
one <- (sum(na.omit(unlist(ones))))/length(tab)
two <- (sum(na.omit(unlist(twos))))/length(tab)


