library(Biostrings)
library(limma)
library(seqinr)

seq_align <- function(seqs_df, path_to_ref, overlap = FALSE)
{
  require(Biostrings)
  require(limma)
  require(seqinr)
  align_df <- data.frame()
  #pair_aln <- list()
  ref <- read.fasta(path_to_ref)
  ref_str <- toupper(sapply(ref, c2s))
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  if (overlap == FALSE)
    {
      split <- paste0(rep("N",10), collapse = "")
      seq_split <- strsplit2(x = seqs_df[,1], split = split)
      seq_all <- data.frame(sequence=paste0(seq_split[,1],seq_split[,2]), hapid = seqs_df[,2])
    } else {
      seq_all <- data.frame(sequence=seqs_df[,1], hapid = seqs_df[,2])
      }
  for (seq_1 in 1:length(seq_all$sequence))
  {
    aln <- pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
    num <- which.max(score(aln))
    patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
    dist <- adist(as.character(patt)[1],as.character(patt)[2])
    df <- data.frame(hapid = seq_all$hapid[seq_1], 
                     hapseq = as.character(patt)[2], 
                     refseq = as.character(patt)[1], 
                     refid = names(patt)[1], 
                     aln_score = score(aln[num]), 
                     dist = dist)
    align_df <- rbind(align_df,df)
    #pair_aln[seq_all$hapid[seq_1]] <- aln[num]
  }
  return(align_df)
}