#!/bin/env R

require(compiler)
enableJIT(3)
     

getvcfchr <- function(vcf, chr, start, end)
  {
    if (length(grep(chr,unique(vcf@fix[,"CHROM"]))) != 1) {
      stop(paste("chromosome",chr,"not found!"))
      }
      
    vcf_chr <- vcf[(vcf@fix[,"CHROM"] == chr)]
    pos <- getPOS(vcf_chr)
    vcf_sub <- vcf_chr[(pos >= (start-100) & pos <= (end+100))]
    
    if (dim(vcf_sub)[1] == 0) {
      warning("No Variant entries found for the specified positions and chromosome!")
      return(NA)
      }
    return(vcf_sub)
  }


annot <- function(vcf_sub, location, fasta, start, end, pos_indels = NULL)
  {
    for (loc in location) {
       if (loc %in% pos_indels) {
         fasta[[1]][loc] <- '[-/-]'
         next
         }
       ref_a <- vcf_sub@fix[(vcf_sub@fix[,"POS"] == loc),"REF"]
       alt_a <- vcf_sub@fix[(vcf_sub@fix[,"POS"] == loc),"ALT"]
       if (nchar(alt_a) > 1) {
         #Split alternate alleles into different nucleotides
         alt_a <- paste(strsplit(alt_a,",")[[1]], collapse = "/")
         }
       if (ref_a != toupper(fasta[[1]][loc])) {
         warning(paste("Possible misalignment at location",loc)) #Needs more refinement # Add option for accepting indels
         }
       if (loc >= start & loc <= end) {
         fasta[[1]][loc] <- paste0('[',fasta[[1]][loc],'/',alt_a,']') 
         } else {
             fasta[[1]][loc] <- 'N'
             }
       }
  return(fasta)
  }

getvcfchr <- cmpfun(getvcfchr)
annot <- cmpfun(annot)

#pf3d7 <- read.fasta("reference.fasta")
#load("vcf.Rdata")
#vcf <- read.vcfR("reference.vcf")

polyvar <- function(chr, start, end, ref, vcf, id, fasta.out, keep.indel = FALSE, fil)
  {
    require(seqinr)
    require(vcfR)
    
    if (class(chr) != "character") {
      stop("argument 'chr' should be of type character!")
      }
    
    if (class(c(start,end)) != "integer") {
      stop("argument 'start' and 'end' must be of type integer")
      }
    
    ref.fasta <- read.fasta(ref)
    fasta <- ref.fasta[chr]
    ann <- getAnnot(fasta)[[1]]    
    ann <- paste0(substr(ann,2,nchar(ann))," | start=",(start-100)," | end=",(end+100)," | ID=",id)
        
    vcf_sub <- getvcfchr(vcf, chr, start, end)
    if (is.na(vcf_sub)) {
      fasta_sub <- toupper(fasta[[1]][(start-100):(end+100)])
      seqinr::write.fasta(fasta_sub, names = ann, file.out = fasta.out, open = "a")
      return(vcf_sub)
      }

    gt <- extract.gt(vcf_sub)
    gt[is.na(gt)] = "0/0"
    hapc <- rowSums(gt != "0/0")
    vcf_sub <- vcf_sub[hapc >= fil]
    location <- getPOS(vcf_sub)
    
    if (length(location) < 1) {
      warning("No heterozygous locations found for the given sequence range!")
      vcf_sub <- NA
      fasta_sub <- toupper(fasta[[1]][(start-100):(end+100)])
      seqinr::write.fasta(fasta_sub, names = ann, file.out = fasta.out, open = "a")
      return(vcf_sub)
      } 
    
    if (!keep.indel) {
	  indels <- extract.indels(vcf_sub, return.indels = TRUE)
	  if (dim(indels@fix)[1] > 0) {
	    pos_indels <- getPOS(indels)
	    if (sum(!(pos_indels >= start & pos_indels <= end)) > 0) {
	      return(NULL)
	      }
	    fasta_new <- annot(vcf_sub, location, fasta, start, end, pos_indels)
	    }
	  }    
 
    #location <- as.numeric(vcf_sub@fix[(!is.na(vcf_sub@fix[,"ALT"])),"POS"])
        
    if (!exists("fasta_new")) {
      fasta_new <- annot(vcf_sub, location, fasta, start, end)
      }
      
    vcf_tidy <- vcfR2tidy(vcf_sub[(!is.na(vcf_sub@fix[,"ALT"]))])
    fasta_sub <- toupper(fasta_new[[1]][(start-100):(end+100)])
    
    seqinr::write.fasta(fasta_sub, names = ann, file.out = fasta.out, open = "a")
    
    return(vcf_tidy)
  }

polyvar <- cmpfun(polyvar)
