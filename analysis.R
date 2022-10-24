# Program: SarsCov2Origin
# Version: 1.0 alpha
# Date: 05/06/2022
# Author: Emiliano Vivas Rodríguez
# Contact: a01424732@tec.mx

install.packages(seqinr);
install.packages(ggplot2);
install.packages(RecordLinkage);
library(seqinr);
library(ggplot2);
library(RecordLinkage);

AMINOACIDS = c(UUU="FPhe", UUC="FPhe", UUA="LLeu", UUG="LLeu", UCU="SSer", UCC="SSer", UCA="SSer", UCG="SSer", UAU="YTyr", UAC="YTyr", UAA='-', UAG='-', UGU="CCys", UGC="CCys", UGA='-', UGG="WTrp", CUU="LLeu", CUC="LLeu", CUA="LLeu", CUG="LLeu", CCU="PPro", CCC="PPro", CCA="PPro", CCG="PPro", CAU="HHis", CAC="HHis", CAA="QGln", CAG="QGln", CGU="RArg", CGC="RArg", CGA="RArg", CGG="RArg", AUU="IIle", AUC="IIle", AUA="IIle", AUG="MMet", ACU="TThr", ACC="TThr", ACA="TThr", ACG="TThr", AAU="NAsn", AAC="NAsn", AAA="KLys", AAG="KLys", AGU="SSer", AGC="SSer", AGA="RArg", AGG="RArg", GUU="VVal", GUC="VVal", GUA="VVal", GUG="VVal", GCU="AAla", GCC="AAla", GCA="AAla", GCG="AAla", GAU="DAsp", GAC="DAsp", GAA="EGlu", GAG="EGlu", GGU="GGly", GGC="GGly", GGA="GGly", GGG="GGly");
SARSCov2 = read.fasta("sarscov2.txt");
RaTG13 = read.fasta("batcovratg13.txt");
MERSCov = read.fasta("merscov.txt");
finalVirus = "SARS-Cov 2";
originVirus = c("RaTG13", "MERS-Cov");
data = data.frame(virus = character(), gene = character(), protein = character(), length = integer(), initialLocation = integer(), finalLocation = integer(), nucleotids = character(), aminoacids = character());
mutations = data.frame(mutation = character(), nucleotide = numeric(), codon = character(), protein = character(), index = numeric(), structure = character(), gene = character(), change = logical());
mutations2 = data.frame(mutation = character(), nucleotide = numeric(), codon = character(), protein = character(), index = numeric(), structure = character(), gene = character(), change = logical());
similarity = data.frame(type=character(), originVirus=character(), finalVirus=character(), gene=character(), nucleotidePercentage = numeric(), aminoacidPercentage=numeric());

needleman = function(seq1, seq2, gap, mismatch, match){
  stopifnot(gap <= 0);
  stopifnot(mismatch <= 0);
  stopifnot(match >= 0);
  len1 = nchar(seq1); len2 = nchar(seq2);
  seq1 = unlist(strsplit(seq1, split = ""));
  seq2 = unlist(strsplit(seq2, split = ""));
  M = matrix(0, nrow = len1 + 1, ncol = len2 + 1);
  rownames(M) = c("-", seq1);
  colnames(M) = c("-", seq2);
  M[1, ] = cumsum(c(0, rep(gap, len2)));
  M[, 1] = cumsum(c(0, rep(gap, len1)));
  D = matrix(0, nrow = len1 + 1, ncol = len2 + 1);
  rownames(D) = c("-", seq1);
  colnames(D) = c("-", seq2);
  D[1, ] = rep("\u2190");
  D[, 1] = rep("\u2191");
  type = c("\u2190", "\u2191", "\u2b09");
  for (i in 2:(len1 + 1)){
    for (j in 2:(len2 + 1)){
      hor = M[i, j - 1] + gap;
      ver = M[i - 1, j] + gap;
      dia = ifelse(rownames(M)[i] == colnames(M)[j], M[i - 1, j - 1] + match, M[i - 1, j - 1] + mismatch);
      M[i, j] = max(hor, ver, dia);
      D[i, j] = type[which.max(c(hor, ver, dia))];
    }
  } 
  align1 = c();align2 = c();
  while(i > 1 && j > 1){
    if(D[i, j] == "\u2b09") {
      align1 = c(rownames(M)[i], align1);
      align2 = c(colnames(M)[j], align2);
      j = j - 1; i = i - 1;
    } else if (D[i, j] == "\u2191") {
      align1 = c(rownames(M)[i], align1);
      align2 = c("-", align2);
      i = i - 1;
    } else if (D[i, j] == "\u2190") {
      align1 = c("-", align1);
      align2 = c(colnames(M)[j], align2);
      j = j - 1;
    } 
  }
  return(list(aligned_seqs = matrix(c(align1, align2), byrow = TRUE, nrow = 2), score = M[nrow(M), ncol(M)], score_matrix = M, movement_matrix = D));
}

transformingToAminoacids = function(sequence){
  newSequence = "";
  for(index in seq(1, length(sequence), 3))
    newSequence = paste(newSequence, substr(AMINOACIDS[paste(sequence[index], sequence[index+1], sequence[index+2], sep = '', collapse = '')], 1, 1), sep = '', collapse = '');
  return(newSequence);
}

gettingMetadata = function(sequence){
  sequence = unlist(strsplit(attr(sequence,"Annot"), "\\[|\\]|:|=|\\.|join|\\(|\\)|,"));
  return(sequence[sequence!="" & sequence!=" "]);
}

registeringData = function(genomeDenomination, genome){
  for (index in 1:length(genome)){
    sequence = toupper(genome[[index]]);
    sequence[sequence == 'T'] = 'U';
    metadata = gettingMetadata(sequence);
    data[nrow(data) + 1, ] <<- list(genomeDenomination, metadata[which(metadata=="gene") + 1], metadata[which(metadata=="protein") + 1], as.integer(metadata[which(metadata=="location") + 2]) - as.integer(metadata[which(metadata=="location") + 1]) + 1, as.integer(metadata[which(metadata=="location") + 1]), as.integer(metadata[which(metadata=="location") + 2]), paste(sequence, sep = '', collapse = ''), transformingToAminoacids(sequence));
  }
}

analysis = function(indexOriginVirus, index, index2, final){
  if(data[index, ]$gene==data[index2, ]$gene){
    origin = data[index2, ]$nucleotids;
    solution = needleman(data[index, ]$aminoacids, data[index2, ]$aminoacids, gap = -1, mismatch = -1, match = 0);
    similarity[nrow(similarity) + 1, ] <<- list(paste(data[index, ]$gene, originVirus[indexOriginVirus], finalVirus, sep="_", collapse = ""), originVirus[indexOriginVirus], finalVirus, data[index, ]$gene, round(levenshteinSim(data[index, ]$nucleotids, data[index2, ]$nucleotids)*100,3), round(levenshteinSim(data[index, ]$aminoacids, data[index2, ]$aminoacids)*100, 3));
    similarity[nrow(similarity) + 1, ] <<- list(paste(data[index, ]$gene, originVirus[indexOriginVirus], finalVirus, sep="_", collapse = ""), originVirus[indexOriginVirus], "Alineamiento A", data[index, ]$gene, "", round(levenshteinSim(data[index, ]$aminoacids, paste(solution$aligned_seqs[1,], sep = '', collapse = ''))*100,3));
    similarity[nrow(similarity) + 1, ] <<- list(paste(data[index, ]$gene, originVirus[indexOriginVirus], finalVirus, sep="_", collapse = ""), finalVirus, "Alineamiento A", data[index, ]$gene, "", round(levenshteinSim(data[index2, ]$aminoacids, paste(solution$aligned_seqs[1,], sep = '', collapse = ''))*100, 3));
    similarity[nrow(similarity) + 1, ] <<- list(paste(data[index, ]$gene, originVirus[indexOriginVirus], finalVirus, sep="_", collapse = ""), originVirus[indexOriginVirus], "Alineamiento B", data[index, ]$gene, "", round(levenshteinSim(data[index, ]$aminoacids, paste(solution$aligned_seqs[2,], sep = '', collapse = ''))*100, 3));
    similarity[nrow(similarity) + 1, ] <<- list(paste(data[index, ]$gene, originVirus[indexOriginVirus], finalVirus, sep="_", collapse = ""), finalVirus, "Alineamiento B", data[index, ]$gene, "", round(levenshteinSim(data[index2, ]$aminoacids, paste(solution$aligned_seqs[2,], sep = '', collapse = ''))*100, 3));
    for(index3 in 1:nchar(final)){
      if(substr(final, index3, index3)!=substr(origin, index3, index3)){
        codon = ceiling(index3/3);
        finalCodon = substr(final, 3*codon-2, 3*codon);
        finalCodonDenomination = AMINOACIDS[paste(finalCodon, sep = '', collapse = '')];
        originCodon = substr(origin, 3*codon-2, 3*codon);
        originCodonDenomination = AMINOACIDS[paste(originCodon, sep = '', collapse = '')];
        if(indexOriginVirus==1) {
          mutations[nrow(mutations)+1, ] <<- list(paste(substr(origin, index3, index3), "to", substr(final, index3, index3), sep = '', collapse = ''), data[index, ]$initialLocation + index3, paste(paste(originCodon, sep = '', collapse = ''), "to", paste(finalCodon, sep = '', collapse = ''), sep = '', collapse = ''), paste(substr(originCodonDenomination, 1, 1), "to", substr(finalCodonDenomination, 1, 1), sep = '', collapse = ''), codon, data[index, ]$protein, data[index, ]$gene, is.na(substr(originCodonDenomination, 1, 1)) | substr(finalCodonDenomination, 1, 1) != substr(originCodonDenomination, 1, 1));
        } else {
          mutations2[nrow(mutations2)+1, ] <<- list(paste(substr(origin, index3, index3), "to", substr(final, index3, index3), sep = '', collapse = ''), data[index, ]$initialLocation + index3, paste(paste(originCodon, sep = '', collapse = ''), "to", paste(finalCodon, sep = '', collapse = ''), sep = '', collapse = ''), paste(substr(originCodonDenomination, 1, 1), "to", substr(finalCodonDenomination, 1, 1), sep = '', collapse = ''), codon, data[index, ]$protein, data[index, ]$gene, is.na(substr(originCodonDenomination, 1, 1)) | substr(finalCodonDenomination, 1, 1) != substr(originCodonDenomination, 1, 1));
        }
      }
    }
  }
}

registeringData(finalVirus, SARSCov2);
registeringData(originVirus[1], RaTG13);
registeringData(originVirus[2], MERSCov);
for(index in 1:length(SARSCov2)){
  for(index2 in (length(SARSCov2)+1):(length(SARSCov2)+length(RaTG13)))
    analysis(1, index, index2, data[index, ]$nucleotids);
  for(index2 in (length(SARSCov2)+length(RaTG13)+1):(length(SARSCov2)+length(RaTG13)+length(MERSCov)))
    analysis(2, index, index2, data[index, ]$nucleotids);
}

graphic1a = ggplot(mutations) + 
  aes(x=mutation, fill=mutation, label=..count..) +
  ggtitle(paste("Diferencias totales entre ", originVirus[1], " y ", finalVirus, '.', sep = '', collapse = '')) +
  labs(x="Diferencia (pseudomutación).", y="Frecuencia.", fill="Pseudomutaciones.") +
  geom_bar(stat = "count") +
  geom_text(stat = "count", aes(vjust=1)) +
  theme(axis.title=element_text(size=10,face="bold"));

graphic2a = ggplot(mutations[mutations$change==TRUE, ]) + 
  aes(x=mutation, fill=mutation, label=..count..) + 
  ggtitle(paste("Diferencias significativas entre ", originVirus[1], " y ", finalVirus, '.', sep = '', collapse = '')) + 
  labs(x="Diferencia (pseudomutación).", y="Frecuencia.", fill="Pseudomutaciones.") + 
  geom_bar(stat = "count") + 
  geom_text(stat = "count", aes(vjust=1)) + 
  theme(axis.title=element_text(size=10,face="bold"));

graphic3a = ggplot(similarity[similarity$nucleotidePercentage>0 & similarity$originVirus==originVirus[1], ], aes(x = gene, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud de genes estructurales entre ", originVirus[1], " y ", finalVirus, '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic4a = ggplot(similarity[similarity$originVirus==originVirus[1] & similarity$gene=='S', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen S de ", originVirus[1], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic5a = ggplot(similarity[similarity$originVirus==originVirus[1] & similarity$gene=='N', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen N de ", originVirus[1], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic6a = ggplot(similarity[similarity$originVirus==originVirus[1] & similarity$gene=='M', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen M de ", originVirus[1], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic7a = ggplot(similarity[similarity$originVirus==originVirus[1] & similarity$gene=='E', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen E de ", originVirus[1], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic1b = ggplot(mutations2) + 
  aes(x=mutation, fill=mutation, label=..count..) +
  ggtitle(paste("Diferencias totales entre ", originVirus[2], " y ", finalVirus, '.', sep = '', collapse = '')) +
  labs(x="Diferencia (pseudomutación).", y="Frecuencia.", fill="Pseudomutaciones.") +
  geom_bar(stat = "count") +
  geom_text(stat = "count", aes(vjust=1)) +
  theme(axis.title=element_text(size=10,face="bold"));

graphic2b = ggplot(mutations2[mutations$change==TRUE, ]) + 
  aes(x=mutation, fill=mutation, label=..count..) + 
  ggtitle(paste("Diferencias significativas entre ", originVirus[2], " y ", finalVirus, '.', sep = '', collapse = '')) + 
  labs(x="Diferencia (pseudomutación).", y="Frecuencia.", fill="Pseudomutaciones.") + 
  geom_bar(stat = "count") + 
  geom_text(stat = "count", aes(vjust=1)) + 
  theme(axis.title=element_text(size=10,face="bold"));

graphic3b = ggplot(similarity[similarity$nucleotidePercentage>0 & similarity$originVirus==originVirus[2], ], aes(x = gene, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud de genes estructurales entre ", originVirus[2], " y ", finalVirus, '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic4b = ggplot(similarity[similarity$originVirus==originVirus[2] & similarity$gene=='S', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen S de ", originVirus[2], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic5b = ggplot(similarity[similarity$originVirus==originVirus[2] & similarity$gene=='N', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen N de ", originVirus[2], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic6b = ggplot(similarity[similarity$originVirus==originVirus[2] & similarity$gene=='M', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen M de ", originVirus[2], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

graphic7b = ggplot(similarity[similarity$originVirus==originVirus[2] & similarity$gene=='E', ], aes(x = finalVirus, y = aminoacidPercentage, fill = type)) +
  ggtitle(paste("Porcentaje de similitud del gen E de ", originVirus[2], '.', sep = '', collapse = '')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = aminoacidPercentage), vjust = 1, hjust=2, colour = "white") +
  theme(legend.position = "bottom") +
  theme(axis.title=element_text(size=10,face="bold")) +
  coord_flip();

save.image("analysis.RData");
