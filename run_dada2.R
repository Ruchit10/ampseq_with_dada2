#!/bin/R env

library(dada2)
library(limma)
library(seqinr)

# Universal parameters
path = "./fastq/filtered/cutPrimer"
filter = TRUE
run_name="ampseq_val17_all_cutp_trim14.Rdata"

# Parameters for filtering
maxEE = c(5,5)
truncLen = c(185,125)
truncQ = c(5,5)

# Parameters for Denoising
max_consist=10
randomize=TRUE
selfConsist=TRUE
omega_a=1e-120

parameter_df <- data.frame(maxEE=maxEE,
	truncLen=truncLen,
	truncQ=truncQ,
	max_consist=max_consist,
	randomize=randomize,
	selfConsist=selfConsist,
	OMEGA_A=omega_a)

# Begin

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
fnames <- sort(list.files(path, pattern="_R1.fastq.gz"))

if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
	stop("fastq files incomplete or not found")
}

sample.names <- strsplit2(fnames, "_R1")[,1]

filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt_R1.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_filt_R2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if (filter == TRUE) {
	print("filtering samples...")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
            maxN=0, maxEE=maxEE, truncLen=truncLen, truncQ=truncQ, rm.phix=TRUE,
            compress=TRUE, multithread=TRUE, verbose=TRUE)
	print("filtering done!")
} else {
	print("skipping filter except mandatory removal of N's... ")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=c(0,0), maxN=0, rm.phix=TRUE,
            compress=TRUE, multithread=TRUE, verbose=TRUE)
}


print("starting error model learning for forward reads...")

errF <- learnErrors(filtFs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)

print("starting error model learning for reverse reads...")

errR <- learnErrors(filtRs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)

print("starting dada2 for forward reads...")

dadaFs <- dada(filtFs, err=errF, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)

print("starting dada2 for reverse reads...")

dadaRs <- dada(filtRs, err=errR, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)

print("merging paird ends...")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

print("generating sequence table...")

seqtab <- makeSequenceTable(mergers)

print("removing chimeric sequences...")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

print("done!")

# Save Progress
save.image(run_name)