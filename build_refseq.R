#!/bin/env R

require(seqinr)
require(vcfR)
require(data.table)


getvcfchr <- function(vcf, chr, start, end)
  {
    if (length(grep(chr,unique(vcf@fix[,"CHROM"]))) != 1) {
      stop(paste("chromosome",chr,"not found!"))
      }

    vcf_chr <- vcf[(vcf@fix[,"CHROM"] == chr)]
    gt <- extract.gt(vcf_chr)
    gt <- gt[!is.na(gt)]
    vcf_gt <- vcf_chr[gt != "0/0"]
    pos <- getPOS(vcf_gt)
    vcf_sub <- vcf_gt[(pos >= start & pos <= end)]

    if (dim(vcf_sub)[1] == 0) {
      warning("No Variant entries found for the specified positions and chromosome!")
      return(NA)
      }
    return(vcf_sub)
  }


annot <- function(vcf_sub, fasta, start, end)
  {
    location <- getPOS(vcf_sub)
    fasta_alt <- fasta
    for (loc in location) {
      if (loc >= start & loc <= end) {
        gt <- extract.gt(vcf_sub[(vcf_sub@fix[,"POS"] == loc)])
        if (gt == "0/0") {
          next
        }
        gt_num <- mean(as.numeric(strsplit(gt,"/")[[1]]))
        if ((gt_num %% 1) != 0) {
          warning(paste("Heterozygous site found at location", loc))
          next
        }
        ref_a <- vcf_sub@fix[(vcf_sub@fix[,"POS"] == loc),"REF"]
        alt_a <- vcf_sub@fix[(vcf_sub@fix[,"POS"] == loc),"ALT"]
        if (nchar(alt_a) > 1) {
          #Split alternate alleles by genotype
          alt_a <- strsplit(alt_a,",")[[1]][gt_num]
          if (alt_a == "*") {
            alt_a = ""
          }
        }
        if (nchar(ref_a) > 1) {
          refseq <- fasta[[1]][(loc):(loc+nchar(ref_a)-1)]
          if (toupper(c2s(refseq)) != ref_a) {
            stop(paste("Possible misalignment between location",loc, "and", loc+nchar(ref_a)-1))
          }
          altseq <- substring(alt_a,1:nchar(ref_a),c(1:(nchar(ref_a)-1),nchar(alt_a)))
          fasta_alt[[1]][(loc):(loc+nchar(ref_a)-1)] <- altseq
        } else {
          if (ref_a != toupper(fasta[[1]][loc])) {
            stop(paste("Possible misalignment at location",loc))
          }
          fasta_alt[[1]][loc] <- alt_a
        }
      } else {
        stop("variant location falling outside specified range!")
      }
    }
    fasta_sub <- toupper(fasta_alt[[1]][(start):(end)])
    return(fasta_sub)
  }

getseq <- function(chr, start, end, ref, id, fasta.out, vcf = NULL, vcf.out = NULL)
  {

    if (class(chr) != "character") {
      stop("argument 'chr' should be of type character!")
      }

    if (class(c(start,end)) != "integer") {
      stop("argument 'start' and 'end' must be of type integer")
      }

    ref.fasta <- read.fasta(ref)
    fasta <- ref.fasta[chr]
    fasta_sub <- toupper(fasta[[1]][(start):(end)])
    #seqinr::write.fasta(fasta_sub, names = id, file.out = fasta.out, open = "a", nbchar = 500)

    if (!is.null(vcf)) {
      vcf_sub <- getvcfchr(vcf, chr, start, end)
      if (is.na(vcf_sub)) {
        return(vcf_sub)
      }

      fasta_new <- annot(vcf_sub, fasta, start, end)
      if (is.null(vcf.out)) {
        stop("please specify secondary output filename if a vcf file is provided!")
      }
      seqinr::write.fasta(fasta_new, names = id, file.out = vcf.out, open = "a", nbchar = 500)
      vcf_tidy <- as.data.frame(vcfR2tidy(vcf_sub)$fix)
      return(vcf_tidy)
    } else {
      return(0)
    }
  }

# Main
mpath="all_amplicons_list.txt"
refpath="/seq/plasmodium/panchal/reference_genome/pf3d7_genome.fasta"
ampfile <- fread(mpath, sep = "\t", header = TRUE)
chr_list <- unique(ampfile$chr)
for (chrm in chr_list) {
  #vcf_path = paste0("/seq/plasmodium/data/pf3k/5.1/vcf/SNP_INDEL_",chr,".combined.filtered.vcf.gz")
  vcf_path = paste0(chrm,"_final.recode.vcf")
  vcf <- read.vcfR(vcf_path)
  ampfile_chr <- ampfile[ampfile$chr == chrm,]
  for (en in 1:dim(ampfile_chr)[1]) {
    vcf_tidy <- getseq(chr = chrm,
                    start = ampfile_chr$start[en],
                    end = ampfile_chr$end[en],
                    ref = refpath,
                    id = ampfile_chr$id[en],
                    fasta.out = NULL,
                    vcf = vcf,
                    vcf.out = "pfdd2_ref_snp_indel.fasta")
    if (is.na(vcf_tidy)) {
      paste("NA returned")
      next
    }
    if (vcf_tidy != 0) {
      write.table(vcf_tidy, paste0(ampfile_chr$id[en],"_output.txt"), sep = "\t", quote = F, row.names = F)
    }
  }
}
