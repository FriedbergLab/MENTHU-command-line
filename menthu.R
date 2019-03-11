#!/usr/bin/Rscript

###############################

#### Usage ####
# Rscript menthu.R ["CRISPR Option"] ["PAM Sequence"] [Distance to DSB] [Overhang] ["TALEN Option"] ["TALEN scheme"] ["Gen Input Type"] ["Gen Input"] [Score Threshold] ["T7 opt"] [silent] [validate]
# Example SpCas9, Ensembl:
# Rscript menthu.R 'T' 'NGG' -3 0 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'
# Example Cas12a:
# Rscript menthu.R 'T' 'TTTN' 18 5 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'
# Example TALEN (left arm: 15 nts, spacer: 16 nts; right arm 18 nts):
# Rscript menthu.R 'F' NA 0 0 'T' '15/16/18' 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'
###
# Parameter Explanation
# Parameter:				Accepted Values:						Explanation:
# outFile						character string						The name of the file to output your results to. If using a fasta file with multiple sequences, multiple files will be created, using this as a prefix
# "CRISPR Option"		"T" or "F"									Flags the system to use CRISPR nuclease processing. If this option is "T" (true), "TALEN Option" must be "F" (false)  
# "PAM Sequence"		A PAM sequence							The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are allowed. Using "N" will scan every possible 
#																										cut site in the target sequence
# "Distance to DSB"	Integer											The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, use a negative value (e.g., "-3" for SpCa9); 
#																										for downstream, use a positive value (e.g., "18" for Cas12a.)
# "Overhang"				Integer >= 0								The length of 5' overhang produced by the nuclease (e.g., "5" for Cas12a). Use "0" for blunt-cutting nucleases. 
# "TALEN Option"		"T" or "F"									Flags the system to use TALEN processing. If this option is "T" (true), "CRISPR Option" must be "F" (false)
# "TALEN Scheme"		"15-18/14 or 16/15-18"			The left arm length, spacer length, and right arm length to use when searching for TALEN locations. E.g., for a TALEN
#																										with arms 15 nt long, with spacer 14 nt, use "15/14/15". TALEN arms can be 15-18 nt in length; the spacer should be 14 OR 16 nt
#																										in length (15 is not allowed for the spacer)
# "Gen Input Type"	"gb" "ens" "seq" "file"			Flags the system to get a GenBank/RefSeq ID ("gb"), Ensembl ID "ens", DNA sequence ("seq"), or to expect a FASTA file ("file")	
# "Gen Input"				See explanation							Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for "file" option, or DNA sequence for "seq"
# "Score Threshold"	Positive number							Only output results with MENTHU score above this threshold. Default is 1.0. Recommended only use sites with score >= 1.5
# "T7 opt"					"T" or "F"									If "T" (true), only displays results where the gRNA is compatible with T7-cloning. Default is "F" - displays all sites.
# "silent"					"T" or "F"									If "T" (true), does not output progress messages to the console. Default is "F" - output progress messages
# "validate"				"T" or "F"									If "T" (true), checks the command line arguments to make sure they are all valid (this may take some time); if "F", skip validation checks


args <- commandArgs(TRUE)

# Get command line arguments
outFile    <- args[1]
crisprFlag <- args[2]
pam        <- args[3]
distDSB    <- args[4]
oh         <- args[5]
talFlag    <- args[6]
tal        <- args[7]
glType     <- args[8]
gl         <- args[9]
threshold  <- args[10]
t7Flag     <- args[11]
silent     <- args[12]
validate   <- args[13]

# Source required supporting files
source("p_apeShiftFunctions.R")
source("p_genbankAccessoryFunctions.R")
source("p_menthu2.0AccessoryFunctions.R")
source("p_required2.0Functions_1.R")
source("p_targetAccessoryFunctions2.0.R")
source("p_ensemblAccessoryFunctions.R")
source("p_validationFunctions.R")

suppressMessages(library(stringr))
suppressMessages(library(stringi))
suppressMessages(library(Biostrings))
suppressMessages(library(rentrez))
suppressMessages(library(rlist))
suppressMessages(library(plyr))
suppressMessages(library(Rcpp))
suppressMessages(library(curl))
suppressMessages(library(httr))
suppressMessages(library(jsonlite))
suppressMessages(library(xml2))

#### Check all the inputs before proceeding ####
validate <- toupper(as.character(validate))
errMess  <- c()
warnMess <- c()

# If the "validate" parameter is true, validate the other inputs
if(as.logical(validate)) {
	allValid <- FALSE
	
	if(!silent) { print("Validating arguments...") } # Progress update
	
	# Check the CRISPR and TALEN information
	if(crisprFlag && talFlag) {
		warnMess <- c(warnMess, "Warning: Both CRISPR and TALEN nuclease flags were set to true; only the CRISPR option will be used.")
		talFlag <- FALSE
	}
	
	if(crisprFlag && !talFlag) {
		pam <- as.character(pam)
		
		# Check for non-IUPAC characters
		if(grepl("[^ACGTRYSWKMBDHVNUacgtryswkmbdhvnu]", pam, perl = TRUE)) {
			errMess <- c(errMess, "Error: Non-IUPAC nucleotide codes detected in PAM sequence.")
		}
		
		# Check that distDSB is an integer
		distDSB <- suppressWarnings(as.numeric(distDSB))
		if(is.na(distDSB)) {
			# Check that distDSB is a number
			errMess <- c(errMess, "Error: 'Distance to double-strand break' parameter must be a (whole) number.")
		} else if(distDSB %% 1 != 0) {
			# Check that it's a whole number
			errMess <- c(errMess, "Error: 'Distance to double-strand break' parameter must be a WHOLE number (integer).")
		}
		
		# Check that the overhang is a positive integer
		oh <- suppressWarnings(as.numeric(oh))
		if(is.na(oh)) {
			# Check that it is a number
			errMess <- c(errMess, "Error: 'Overhang' parameter must be a (whole) number greater than or equal to zero.")
		} else {
			# Check for whole number
			if(oh %% 1 != 0) {
				errMess <- c(errMess, "Error: 'Overhang' parameter must be a WHOLE number..")
			}
			
			# Check for >= 0
			if(oh < 0) {
				errMess <- c(errMess, "Error: 'Overhang' parameter must be greater than or equal to zero.")
			}
		}
		
	} else if(!crisprFlag && talFlag) {
		tal <- as.character(tal)
		
		# Check the formatting on the TALENs
		tals <- unlist(strsplit(tal, "/"))
		if(length(tals) != 3) {
			errMess <- c(errMess, "Error: 'TALEN' input parameter should be specified as 'left arm length/spacer length/right arm length', separated by slashes - e.g., '15/16/18'.")
		}
		
		suppressWarnings(tL  <- as.numeric(tals[1]))
		suppressWarnings(spa <- as.numeric(tals[2]))
		suppressWarnings(tR  <- as.numeric(tals[3]))
		
		tals <- c(tL, spa, tR)
		# Check valid inputs on the TALENs
		if(any(is.na(tL), is.na(spa), is.na(tR))) {
			errMess <- c(errMess, "Error: The TALEN arm lengths and spacer lengths should be whole numbers, e.g. '15/16/15'.")
			
		} else {
			if(any(tL %% 1 != 0, spa %% 1 != 0, tR %% 1 != 0)) {
				errMess <- c(errMess, "Error: The TALEN arm lengths and spacer lengths must be WHOLE numbers.")
				
			}
			
			if(any(tL > 18, tL < 15, tR > 18, tR < 15)) {
				errMess <- c(errMess, "Error: TALEN arm lengths should be between 15 and 18 nucleotides.")
			}
			
			if(any(spa != 14, spa != 16)) {
				errMess <- c(errMess, "Error: Spacer should be 14 or 16 nucleotides.")
			}
		}
		
		# Make sure at least one flag is specified
	} else if(!crisprFlag && !talFlag) {
		errMess <- c(errMess, paste0("Error: Nuclease input type is not recognized. Please indicate whether you are using a CRISPR ,",
																 "(argument 1 = 'T', argument 5 = 'F') or TALEN (argument 1 = 'F', argument 5 = 'T')"))
	}
	
	# Check the DNA sequence inputs
	if (glType != 'gb' || glType != 'ens' || glType != 'file' || glType != 'seq'){
		errMess <- c(errMess, paste0("Error: Unrecognized input type; use 'gb' for GenBank/RefSeq ID, 'ens' for Ensembl ID, ", 
																 "'file' for FASTA file, or 'seq' to specify the sequence in the command line."))
	}
	
	# Validate GenBank input
	if(glType == 'gb') {
		validGB <- validGenBankId(gl)
		
		if(!validGB[1]) {
			errMess <- c(errMess, validGB[2])
		}
		
	# Validate Ensembl input
	} else if(glType == 'ens') {
		validEns <- ensemblIdSpecies(gl, bool = TRUE)
		
		if(!validEns[1]) {
			errMess <- c(errMess, validEns[2])
		}
		
	# Validate file input	
	} else if(glType == 'file') {
		if(!file.exists(gl)) {
			errMess <- c(errMess, paste0("Error: File ", gl, " does not exist"))
		} else {
			fExt <- tools::file_ext(gl)
			if (toupper(fExt) != "FASTA" || toupper(fExt) != 'FA' || toupper(fExt) != 'FNA' || toupper(fExt) != "FFN" || toupper(fExt) != "FRN" || toupper(fExt) != "txt") {
				warnMess <- c(warnMess, paste0("Warning: File extension does not match common Fasta or text extensions."))
			}
		}
		
	# Validate DNA sequence
	} else if(glType == 'seq') {
		validDna <- validDnaSeq(gl)
		
		if(!validDna[1]) {
			errMess <- c(errMess, validDna[2])
		}
	}
	
	# Validate threshold
	suppressWarnings(threshold <- as.numeric(threshold))
	if(is.na(threshold)) {
		errMess <- c(errMess, "Error: Threshold must be a number")
	} else if(threshold < 0) {
		errMess <- c(errMess, "Error: Threshold must be >= 0")
	}
	
	# Validate T7 flag
	t7Flag <- toupper(t7Flag)
	if(!(as.logical(t7Flag) || suppressWarnings(as.logical(as.numeric(t7Flag))))) {
		errMess <- c(errMess, "Error: T7 filter must be 'TRUE' or 'FALSE'")
	} 
	
	# Validate silent flag
	silent <- toupper(silent)
	if(!(as.logical(silent) || suppressWarnings(as.logical(as.numeric(silent))))) {
		errMess <- c(errMess, "Error: 'Silent' parameter must be 'TRUE' or 'FALSE'")
	} 
	
	if (length(errMess) == 0) {
		allValid <- TRUE
	} else {
		allValid <- FALSE
	}
	
	
	# If the validate parameter if false, continue on with your life
} else if (!as.logical(validate)) {
	allValid <- TRUE
	crisprFlag <- as.logical(crisprFlag)
	pam        <- as.character(pam)
	distDSB    <- as.numeric(distDSB)
	oh         <- as.numeric(oh)
	talFlag    <- as.logical(talFlag)
	tal        <- as.character(tal)
	
	if(talFlag) {
		# Format tal
		tals <- unlist(strsplit(tal, "/"))
		
		suppressWarnings(tL  <- as.numeric(tals[1]))
		suppressWarnings(spa <- as.numeric(tals[2]))
		suppressWarnings(tR  <- as.numeric(tals[3]))
		
		tals <- c(tL, spa, tR)
	}

	glType     <- as.character(glType)
	gl         <- as.character(gl)
	threshold  <- as.numeric(threshold)
	t7Flag     <- as.logical(t7Flag)
	silent     <- as.logical(silent)
	
	# If the validate parameter is itself invalid, destroy the universe	
} else {
	allValid <- FALSE
	print("Error: Please use 'T' or 'F' to indicate the 'validate' parameter. (We couldn't validate whether you wanted the inputs validated...)")
}


if(!allValid){
	if(!silent) { 
		print(errMess) 
	} else {
		print(paste0("Error: ", paste(args, collapse = " ")))
	}
	
} else {
	

	
	if(glType == 'gb') {
		info <- suppressWarnings(getGenbankFile(gl))
		
		if(crisprFlag) {
			results <- calculateMENTHUGeneSeqGenBank(pam, distDSB, oh, wiggle = TRUE, wiggleRoom = 39, NULL, info, 1, silent)[[1]]
			
		} else if(talFlag) {
			results <- calculateMENTHUGeneSeqGenBank(NULL, NULL, NULL, wiggle = TRUE, wiggleRoom = 39, tals, info, 1, silent)[[1]]
			
			
		} else {
			
			
		}
		results <- results[which(results$MENTHU_Score >= threshold),]
		if(t7Flag && crisprFlag){
			results <- results[which(!grepl("/", results$Tool_Type, fixed = TRUE)), ]
			results <- results[which(unlist(lapply(1:nrow(results), 
																						 function(x) grepl("((^[G]{1,2})|(^[ACGT][G]))", results$Target_Sequence[x], perl=TRUE, ignore.case=TRUE)))), ]
		}
		
		results <- results[order(-results$MENTHU_Score), ]
		write.table(results, outFile, append = FALSE, quote = FALSE, sep = ",", na = "NOT APPLICABLE", row.names = FALSE, col.names = TRUE)
		
		
	} else if(glType == 'ens') { 
		ensemblInfo <- handleEnsemblInput(gl, wiggle = TRUE, wiggleRoom = 39)
		
		# If the entry is NOT an exon
		if(getEnsemblIdType(gl, check = TRUE) != "exon") {
			
			# Figure out how many exons there are in the transcript/protein
			ensemblInfo$rank <- as.numeric(ensemblInfo$rank)
			
			# Determine how many exons are present
			maxR <- max(ensemblInfo$rank)
			
			exonStuff <- seq(1, maxR)
			
		# If the entry is an exon	
		} else {
			exonStuff <- 1
			
		}

		if(crisprFlag) {
			results <- calculateMENTHUEnsembl(pam, distDSB, oh, wiggle = TRUE, wiggleRoom = 39, NULL, ensemblInfo, exonStuff, silent)[[1]]
			
		} else if(talFlag) {
			results <- calculateMENTHUEnsembl(NULL, NULL, NULL, wiggle = TRUE, wiggleRoom = 39, tals, ensemblInfo, exonStuff, silent)[[1]]
			
		} else {
			
		}
		
		
		results <- results[which(results$MENTHU_Score >= threshold),]
		if(t7Flag && crisprFlag){
			results <- results[which(!grepl("/", results$Tool_Type, fixed = TRUE)), ]
			results <- results[which(unlist(lapply(1:nrow(results), 
																							 function(x) grepl("((^[G]{1,2})|(^[ACGT][G]))", results$Target_Sequence[x], perl=TRUE, ignore.case=TRUE)))), ]
		}
		results <- results[order(-results$MENTHU_Score), ]
		write.table(results, outFile, append = FALSE, quote = FALSE, sep = ",", na = "NOT APPLICABLE", row.names = FALSE, col.names = TRUE)
		
	} else if(glType == 'file') {
		inputFasta <- Biostrings::readDNAStringSet(gl, format = "fasta")
		
		inputSeqs <- as.character(inputFasta)
		
		
		if(crisprFlag) {
			results <- lapply(inputSeqs, function(x) {calculateMENTHUGeneSeq(pam, distDSB, oh, wiggle = TRUE, wiggleRoom = 39, NULL, x, silent)})
			
		} else if(talFlag) {
			results <- lapply(inputSeqs, function(x) {calculateMENTHUGeneSeq(NULL, NULL, NULL, wiggle = TRUE, wiggleRoom = 39, tals, x, silent)})
			
		} else {
			
		}
		
		# Fix the list stuff that happens if the file has > 1 fasta sequences in it
		if(length(results) > 1) {
			fileNames <- unlist(lapply(colnames(results), function(x) {paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", outFile), "_", x, ".csv")}))
		} else {
			fileName <- outFile
		}
		
		for(i in 1:length(fileNames)) {
			results <- results[which(results$MENTHU_Score >= threshold),]
			if(t7Flag && crisprFlag){
				results <- results[which(!grepl("/", results$Tool_Type, fixed = TRUE)), ]
				results <- results[which(unlist(lapply(1:nrow(results), 
																							 function(x) grepl("((^[G]{1,2})|(^[ACGT][G]))", results$Target_Sequence[x], perl=TRUE, ignore.case=TRUE)))), ]
			}
			results <- results[order(-results$MENTHU_Score), ]
			write.table(results[1, i], fileNames[i], append = FALSE, quote = FALSE, sep = ",", na = "NOT APPLICABLE", row.names = FALSE, col.names = TRUE)
		}
		
		
	} else if(glType == 'seq') {
		
		
		if(crisprFlag) {
			results <- calculateMENTHUGeneSeq(pam, distDSB, oh, wiggle = TRUE, wiggleRoom = 39, NULL, gl, silent)[[1]]
			
		} else if(talFlag) {
			results <- calculateMENTHUGeneSeq(NULL, NULL, NULL, wiggle = TRUE, wiggleRoom = 39, tals, stripWhiteSpace(toupper(gl)), silent)[[1]]
			
		} else {
			
		}
		results <- results[which(results$MENTHU_Score >= threshold),]
		if(t7Flag && crisprFlag){
			results <- results[which(!grepl("/", results$Tool_Type, fixed = TRUE)), ]
			results <- results[which(unlist(lapply(1:nrow(results), 
																						 function(x) grepl("((^[G]{1,2})|(^[ACGT][G]))", results$Target_Sequence[x], perl=TRUE, ignore.case=TRUE)))), ]
		}
		results <- results[order(-results$MENTHU_Score), ]
		write.table(results, outFile, append = FALSE, quote = FALSE, sep = ",", na = "NOT APPLICABLE", row.names = FALSE, col.names = TRUE)
	} else {
		
	}
	
	
}


