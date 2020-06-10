### Analyzing DADA2 output

library(seqinr)
library(dada2)
library(data.table)
source("seq_align.R")

#runid="Run1"
dada2Rdata="~/Desktop/ampseq_val17_all_cutp_trim14.Rdata"
path_to_ref="~/ref3d7.fasta"
amplicons=c("CSP", "SERA2")
truth_set=TRUE
path_truth_set=list(CSP = "~/Desktop/dada2_analysis/ampseq_val17/csp_mock_sample_hap_ids.txt",
                    SERA2 = "~/Desktop/dada2_analysis/ampseq_val17/sera2_mock_sample_hap_ids.txt")
true_sample_info_path="~/Downloads/PASEC_mock_sample_hap_composition.txt"

# Only if sample information about true haplotype composition is present
if (truth_set)
{
  exp_sample_comp <- list()
  df <- fread(true_sample_info_path, sep = "\t", header = TRUE)
  df2 <- df[!is.na(df$Parasite_density),]
  df2$Exp_haps[is.na(df2$Exp_haps)] = ""
  for (amplicon in amplicons)
  {
    exp_sample_comp[[amplicon]] <- df2[df2$Amp == tolower(amplicon),c(1,3,4)]
    exp_hapmat <- strsplit2(exp_sample_comp[[amplicon]]$Exp_haps,";")
    exp_sample_comp[[amplicon]]$exp_count <- apply(exp_hapmat,1,function(x) sum(x != ""))
  }
}

# Load saved DADA2 output
if (file.exists(dada2Rdata))
{
  load(dada2Rdata)
} else {
  stop(paste("file",dada2Rdata,"not found!"))
}

# Extract non-chimeric DADA2 Haplotypes
if (exists("seqtab.nochim"))
{
  seqs <- colnames(seqtab.nochim)
  hapid <- paste0("H",1:length(seqs))
  seqs_df <- data.frame(sequence = seqs, hapid = hapid)
  seqtab_haps <- seqtab.nochim
  colnames(seqtab_haps) <- hapid
  row.names(seqtab_haps) <- strsplit2(sample.names,"_trim")[,1]
} else {
  stop("cannot find DADA2 non-chimeric sequence table!")
}

# Map haplotypes to refernce pf3D7
align_df <- seq_align(seqs_df, path_to_ref, overlap = TRUE)

count_table <- list()
sample_info <- list()
sample_metrics <- list()
amplicon_metrics <- list()
amplicon_df <- list()
hap_info <- list()
for (amplicon in amplicons)
{
  amplicon_index <- (align_df$refid == amplicon)
  count_table[[amplicon]] <- seqtab_haps[,amplicon_index]
  amplicon_df[[amplicon]] <- align_df[amplicon_index,c(2,1,4,6)]
  # Truth set mapping
  if (truth_set)
  {
    # Haplotype info
    truehap_df <- seq_align(seqs_df[amplicon_index,], path_truth_set[[amplicon]], overlap = TRUE)
    read_count <- apply(count_table[[amplicon]],2,sum)
    amplicon_metrics[[amplicon]] <- data.frame(dada2_hap = truehap_df[,1],
                                             true_hap = truehap_df[,4],
                                             dist_from_true = truehap_df[,6],
                                             read_count = read_count,
                                             read_proportion = (read_count / sum(read_count)))
    # Sample info
    true_index <- (amplicon_metrics[[amplicon]]$dist_from_true == 0)
    true_count <- count_table[[amplicon]][,true_index]
    total_hap_count <- apply(count_table[[amplicon]],1,function(x) sum(x != 0))
    true_hap_count <- apply(true_count,1,function(x) sum(x != 0))
    sample_info[[amplicon]] <- data.frame(sample = row.names(seqtab_haps),
                                        total_hap_count = total_hap_count,
                                        true_hap_count = true_hap_count,
                                        total_read_count = apply(count_table[[amplicon]],1,sum),
                                        true_read_count = apply(true_count,1,sum))
    sample_metrics[[amplicon]] <- merge(sample_info[[amplicon]],
                                        exp_sample_comp[[amplicon]][,c(1,3,4)],
                                        by.x = "sample", by.y = "solexa_index",
                                        sort = FALSE)
    sample_metrics[[amplicon]]$FP <- (sample_metrics[[amplicon]]$total_hap_count - sample_metrics[[amplicon]]$true_hap_count)
    sample_metrics[[amplicon]]$FN <- (sample_metrics[[amplicon]]$exp_count - sample_metrics[[amplicon]]$true_hap_count)
    sample_metrics[[amplicon]]$sensitivity <- (sample_metrics[[amplicon]]$true_hap_count/sample_metrics[[amplicon]]$exp_count)
    sample_metrics[[amplicon]]$specificity <- (sample_metrics[[amplicon]]$true_hap_count/sample_metrics[[amplicon]]$total_hap_count)
  }
}

metric_summary <- function(sample_metrics, amplicons)
{
  summary_df <- data.frame()
  for (amplicon in amplicons)
  {
    sn_dist <- sample_metrics[[amplicon]]$sensitivity
    sp_dist <- sample_metrics[[amplicon]]$specificity
    errors <- sum(is.infinite(sn_dist)|(!is.nan(sn_dist))&sn_dist>1)
    mean_sn <- mean(sn_dist[(!is.infinite(sn_dist))&(!is.nan(sn_dist))&(sn_dist<=1)])
    mean_sp <- mean(sp_dist[(!is.infinite(sp_dist))&(!is.nan(sp_dist))&(sp_dist<=1)])
    summary_df <- rbind(summary_df,data.frame(amplicon = amplicon,
                                              total_hap = length(colnames(count_table[[amplicon]])),
                                              mean_sn = mean_sn,
                                              mean_sp = mean_sp,
                                              error_count = errors))
  }
  return(summary_df)
}

summary_df <- metric_summary(sample_metrics, amplicons)

# done!