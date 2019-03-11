####Required Functions#######

#' calculateMENTHUGeneSeq
#'
#' @param casList 
#' @param wiggle 
#' @param wiggleRoom 
#' @param geneSeq 
#' @param threshold 
#' @param exonDF 
#' @param progress 
#' @param armin 
#' @param armax 
#' @param spamin 
#' @param spamax 
#'
#' @return
#' @export
#'
#' @examples

calculateMENTHUGeneSeq <- function(pamList, cutDistList, ohList, wiggle = TRUE, wiggleRoom = 39, tal, geneSeq, silent){
	require(Biostrings)
	require(plyr)
	
	# Set variables
	talFlag   <- FALSE
	noPamFlag <- FALSE
	
	geneSeq <- toupper(geneSeq)
	
	# In the case where exons are not specified, treat the whole gene sequence as a single exon
	exonSeqs <- geneSeq
	
	# If the user is using Cas:
	if(length(pamList) > 0) {
		
		if(!silent) {
			print("Scanning for CRISPR target sites...")
		}
		
		pamSites <- pamScan(pamList, 
												cutDistList, 
												ohList,
												exonSeqs,
												exonList   = 1,
												exonStarts = 1, 
												findCut    = TRUE, 
												type       = "cas9", 
												wiggle     = wiggle, 
												wiggleRoom = wiggleRoom)
		
		
		# Count the number of Cas target sites
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		if(!silent) {
			print("Pre-processing CRISPR target sites...")
		}
		
		tempExonDF <- data.frame(Exon_Num  = 1,
														 exonStart = 1,
														 exonEnd   = nchar(exonSeqs),
														 stringsAsFactors = FALSE)
		
		pamSites <- unique(suppressMessages(plyr::join(pamSites, tempExonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Count the number of targets satisfying context condition
		siteCountC <- nrow(pamSites)
		
		if(nrow(pamSites) >= 1){
			## Get the sequence context surrounding the cut site ##
			# Get the context for making the gRNA sequences
			forguide <- unlist(lapply(1:nrow(pamSites), function (x) substr(geneSeq, 
																																			pamSites$CutIndex[x] - 39, 
																																			pamSites$CutIndex[x] + 40)))
			# Get the context for dealing with overhangs
			context <- unlist(lapply(1:nrow(pamSites), function (x) paste0(substr(geneSeq,
																																						pamSites$CutIndex[x] - 39,
																																						pamSites$CutIndex[x]),
																																		 substr(geneSeq,
																																		 			 pamSites$CutIndex[x] + pamSites$ohLen[x] + 1,
																																		 			 pamSites$CutIndex[x] + pamSites$ohLen[x] + 40))))
			
			# Set the sequence context in the frame
			pamSites$seq   <- context
			pamSites$guide <- forguide
			
		} else {
			noPamFlag  <- TRUE
			siteCount  <- 0
			siteCountC <- 0
		}
		
	} else {
		#If the user is NOT using CRISPR/Cas system, set pamFlag to FALSE
		pamSites <- 0
		pamFlag  <- FALSE
		siteCount  <- 0
		siteCountC <- 0
		
	}
	
	if(length(tal) > 0){
		
		
		# If there are no exon inputs, make exon start null
		# Submit talen info to talPal
		talSites <- talPal(exonSeqs,
											 findCut    = TRUE,
											 wiggle     = TRUE,
											 wiggleRoom = 39,
											 tL         = tal[1], 
											 tR         = tal[3], 
											 spa        = tal[2], 
											 exonStarts = NULL,
											 exonList   = 1)
		
		
		# Update progress bar
		if(!silent) {
			print("Pre-processing TALEN target sites...")
		}
		
		
		tempExonDF <- data.frame(Exon_Num  = 1,
														 exonStart = 1,
														 exonEnd   = nchar(exonSeqs),
														 stringsAsFactors = FALSE)
		
		talSites <- unique(suppressMessages(plyr::join(talSites, tempExonDF, by = 'Exon_Num')))
		
		
		
		# Drop target sites where the cut site is not within the exon boundaries
		keepT     <- sapply(1:nrow(talSites), function(x) talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1))
		talSites  <- talSites[keepT, ]
		
		# Identify sites with enough sequence context to do calculations
		talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		talSites      <- talSites[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Get the sequence context surrounding the cut site
		contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(geneSeq, talSites$CutIndex[x] - 39, talSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		talSites$seq = contextT
		
	} else {
		# If TALENs are not used, set talSites to 0
		talSites <- 0
		
	}
	
	exonDF <- data.frame(exonStart        = 1, 
											 exonEnd          = nchar(geneSeq), 
											 stringsAsFactors = FALSE)
	
	if(pamFlag){
		if(!silent) {
			print("Calculating MENTHU v2.0 scores for CRISPR sites...")
		}
		
		# Function to calculate MENTHU v2.0 score
		menthuFrameFunc <- function(x){
			
			return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
		}
		
		# Get the menthu scores
		menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		# Update progress bar
		if(!silent) {
			print("Formatting CRISPR site results...")
		}
		
		row.names(menthuFrame)  <- c()
		colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
		
		# Get the critOne success count
		critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
		
		# Get the critTwo success count
		critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
		
		# Get both success count
		critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
		
		pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
		
		baseCrispr <- sapply(1:nrow(pamSites), 
												 function(x) if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 })
		
		# Generate the PAM sequence
		pam        <- sapply(1:nrow(pamSites), 
												 function(x) (if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])      
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	}
												 }))
		
		# Get the CRISPR target sequence required to target this site, and bold the PAM
		crispr <- sapply(1:length(baseCrispr), function(x) {
			if(pamSites$CutDist[x] < 0) {
				paste0(baseCrispr[x], pam[x])
			} else{
				paste0(pam[x], baseCrispr[x])
			}})
		
		# Create data frame of the current results
		pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																Frame_Shift      = pamSites$frameShift,
																Tool_Type        = pamSites$Target, 
																Strand           = pamSites$Orientation, 
																Exon_ID          = pamSites$Exon_Num, 
																DSB_Location     = pamSites$CutIndex,
																Microhomology    = pamSites$topMH,
																PreMA_Sequence   = pamSites$topDel,
																Context          = pamSites$seq,
																stringsAsFactors = FALSE)
	}
	
	if(talFlag){
		# Update progress bar
		if(!silent) {
			print("Calculating MENTHU v2.0 scores for TALEN sites...")
		}
		
		menthuFrameTFunc <- function(x){
			
			# Return calculation
			return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
		}
		
		# Calculate slope competition on all the context
		menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		# Update progress bar
		if(!silent) {
			print("Formatting TALEN site results...")
		}
		
		# Clean up the resulting data frame
		colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		rownames(menthuFrameT) <- c()
		
		# Merge slope frame and pamSites
		talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
		
		talenGenFunc <- function(talRow){
			dim  <- unlist(strsplit(talRow$Target, "/"))
			arm1 <- as.numeric(dim[1])
			spa1 <- as.numeric(dim[2])
			arm2 <- as.numeric(dim[3])
			
			armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
			spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
			armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
			
			return(paste0(armL, spac, armR))
		}
		
		# Generate the crispr target sequences
		talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
		
		# Clean the data frame
		row.names(talSites) <- c()
		
		# Create a new data frame of the results
		talFormFrame <- data.frame(Target_Sequence = talSites$talen,
															 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
															 Frame_Shift     = talSites$frameShift,
															 Tool_Type       = talSites$Target,
															 Strand          = talSites$Orientation,
															 Exon_ID         = talSites$Exon_Num,
															 DSB_Location    = talSites$CutIndex,
															 Microhomology   = talSites$topMH,
															 PreMA_Sequence  = talSites$topDel,
															 Context         = talSites$seq,
															 stringsAsFactors = FALSE)
		
		# Get the critOne success count
		critOneT <- nrow(talFormFrame[which(talFormFrame$frameShift != "NA"), ])
		
		# Get the critTwo success count
		critTwoT <- nrow(talFormFrame[which(talFormFrame$menthuScore >= 1.5), ])
		
		# Get both success count
		critBothT <- nrow(talFormFrame[which(talFormFrame$menthuScore >= 1.5 & talFormFrame$frameShift != "NA"), ])
	}
	
	if(pamFlag && talFlag){
		critOne   <- critOne  + critOneT
		critTwo   <- critTwo  + critTwoT
		critBoth  <- critBoth + critBothT
		
	} else if(talFlag){
		critOne   <- critOneT
		critTwo   <- critTwoT
		critBoth  <- critBothT
	}
	
	
	#Return frame
	if(pamFlag && talFlag){
		return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(pamFlag && !talFlag) {
		return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(!pamFlag && talFlag) {
		return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} 
}


#' calculateMENTHUGeneSeqGenBank
#'
#' @param pamList 
#' @param talenList 
#' @param gbFlag 
#' @param gbhFlag 
#' @param genbankInfo 
#' @param threshold 
#' @param firstExon 
#' @param exonTargetType 
#' @param exonStuff 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples

calculateMENTHUGeneSeqGenBank <- function(pamList, cutDistList, ohList, wiggle = TRUE, wiggleRoom = 39, talenList, genbankInfo, exonStuff, silent){
	require(plyr)
	require(Biostrings)
	require(curl)
	require(stringr)
	
	if(!silent) {
		print("Processing GenBank accession...")
	}
	
	# Get exon sequences and information through getExon
	exon <- getExon(genbankInfo, wiggle = TRUE, wiggleRoom = 39, exonStuff)
	
	# Get exon indices
	exonInfo <- exon[[1]]
	# Get the exon sequences
	exonSeq  <- exon[[2]]
	# Get the gene sequence
	geneSeq  <- exon[[3]]
	
	exonDF <- data.frame(Exon_Num         = exonInfo$exonNum,
											 exonStart        = exonInfo$start, 
											 exonEnd          = exonInfo$end, 
											 stringsAsFactors = FALSE)
	
	# If the user is using Cas:
	if(length(pamList) > 0){
		# Update progress
		if(!silent) {
			print("Scanning for target sites...")
		}
		
		if(length(exonInfo) > 0){
			# If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
			pamSites <- pamScan(pamList, 
													cutDistList,
													ohList,
													exonSeq, 
													exonList   = exonInfo$exonNum, 
													exonStarts = exonInfo$start, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = TRUE, 
													wiggleRoom = 39)
		} else {
			pamSites <- pamScan(pamList, 
													cutDistList, 
													ohList,
													exonSeq,
													exonList   = "1",
													exonStarts = NULL, 
													findCut    = TRUE, 
													type       = "cas9", 
													wiggle     = wiggle, 
													wiggleRoom = wiggleRoom)
		}
		
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		# Update progress bar
		if(!silent) {
			print("Pre-processing CRISPR target sites...")
		}

		pamSites <- unique(suppressMessages(plyr::join(pamSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		siteCountC <- nrow(pamSites)
		
		# Get the sequence context surrounding the cut site for the gRNA
		forguide <- unlist(lapply(1:nrow(pamSites), function (x) substr(geneSeq, 
																																		pamSites$CutIndex[x] - 39, 
																																		pamSites$CutIndex[x] + 40)))
		# Get the sequence context surrounding the cut site for the overhang score calculation
		context  <- unlist(lapply(1:nrow(pamSites), function (x) paste0(substr(geneSeq, 
																																					 pamSites$CutIndex[x] - 39, 
																																					 pamSites$CutIndex[x]), 
																																		substr(geneSeq, 
																																					 pamSites$CutIndex[x] + pamSites$ohLen[x] + 1, 
																																					 pamSites$CutIndex[x] + pamSites$ohLen[x] + 40))))
		
		# Set the sequence context in the frame
		pamSites$seq   <- context
		pamSites$guide <- forguide
		
	} else {
		# If the user is NOT using Cas, set pamFlag to FALSE
		pamSites   <- 0
		pamFlag    <- FALSE
		
		siteCount  <- 0
		siteCountC <- 0
	}
	
	# If there are TALEN inputs
	if(length(talenList) > 0){
		
		if(!silent) {
			print("Idetifying TALEN targets...")
		}
		# Set all exon starts to the exon starts in the input frame
		# Submit talen info to talPal
		# If there are exon inputs
		if(length(exonInfo > 0)){
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom = 39,
												 tL         = talenList[1], 
												 tR         = talenList[3], 
												 spa        = talenList[2], 
												 exonList   = exonInfo$exonNum,
												 exonStarts = exonInfo$start)
			
		} else {
			talSites <- talPal(exonSeq,
												 findCut    = TRUE,
												 wiggle     = TRUE,
												 wiggleRoom = 39,
												 tL         = talenList[1], 
												 tR         = talenList[3], 
												 spa        = talenList[2],
												 exonList   = "1",
												 exonStarts = NULL)
		}
		
		# Update progress bar
		if(!silent) {
			print("Pre-processing TALEN target sites...")
		}
		progress$inc(0.01, detail = "Pre-processing TALEN target sites...")
		talSites <- unique(suppressMessages(plyr::join(talSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keepT     <- sapply(1:nrow(talSites), function(x) talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1))
		talSites  <- talSites[keepT, ]
		
		# Identify sites with enough sequence context to do calculations
		talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq)))  ] <- TRUE
		talSites      <- talSites[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= nchar(geneSeq))), ]
		
		# Get the sequence context surrounding the cut site
		contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(geneSeq, talSites$CutIndex[x] - 39, talSites$CutIndex[x] + 40)))
		
		# Set the sequence context in the frame
		talSites$seq = contextT
		
	} else {
		# If TALENs are not used, set talSites list to empty
		talSites <- 0
	}
	
	# Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence  = as.character(), 
														MENTHU_Score     = as.numeric(), 
														Frame_Shift      = as.character(), 
														Tool_Type        = as.character(), 
														Strand           = as.character(), 
														Exon_ID          = as.numeric(), 
														DSB_Location     = as.integer(),
														stringsAsFactors = FALSE)

	if(pamFlag){
		if(!silent) {
			print("Calculating MENTHU v2.0 scores for CRISPR sites...")
		}
		
		# Function to calculate MENTHU v2.0 score
		menthuFrameFunc <- function(x){
			
			return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
		}
		
		# Get the menthu scores
		menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		if(!silent) {
			print("Formatting CRISPR site reults...")
		}
		# Update progress bar
		progress$inc(0.01, detail = "Formatting CRISPR site results...")
		
		row.names(menthuFrame)  <- c()
		colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
		
		# Get the critOne success count
		critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
		
		# Get the critTwo success count
		critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
		
		# Get both success count
		critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
		
		pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
		
		# Drop 0s
		pamSites <- pamSites[which(pamSites$menthuScore > 0), ]
		

		baseCrispr <- sapply(1:nrow(pamSites), 
												 function(x) if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 })
		
		# Generate the PAM sequence
		pam        <- sapply(1:nrow(pamSites), 
												 function(x) (if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])      
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	}
												 }))
		
		# Get the CRISPR target sequence required to target this site, and bold the PAM
		crispr <- sapply(1:length(baseCrispr), function(x) {
			if(pamSites$CutDist[x] < 0) {
				paste0(baseCrispr[x], pam[x])
			} else{
				paste0(pam[x], baseCrispr[x])
			}})
		
		# Create data frame of the current results
		pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																Frame_Shift      = pamSites$frameShift,
																Tool_Type        = pamSites$Target, 
																Strand           = pamSites$Orientation, 
																Exon_ID          = pamSites$Exon_Num, 
																DSB_Location     = pamSites$CutIndex,
																Microhomology    = pamSites$topMH,
																PreMA_Sequence   = pamSites$topDel,
																Context          = pamSites$seq,
																stringsAsFactors = FALSE)
	}
	
	if(talFlag){
		if(!silent) {
			print("Calculating MENTHU v2.0 scores for TALEN sites...")
		}
		
		menthuFrameTFunc <- function(x){
			
			# Return calculation
			return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
		}
		
		# Calculate competition on all the context
		menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		if(!silent) {
			print("Formatting TALEN site results...")
		}
		
		# Clean up the resulting data frame
		colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		rownames(menthuFrameT) <- c()
		
		# Merge slope frame and pamSites
		talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
		
		talenGenFunc <- function(talRow){
			dim  <- unlist(strsplit(talRow$Target, "/"))
			arm1 <- as.numeric(dim[1])
			spa1 <- as.numeric(dim[2])
			arm2 <- as.numeric(dim[3])
			
			armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
			spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
			armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
			
			return(paste0(armL, spac, armR))
		}
		
		# Generate the crispr target sequences
		talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
		
		# Clean the data frame
		row.names(talSites) <- c()
		
		# Create a new data frame of the results
		talFormFrame <- data.frame(Target_Sequence = talSites$talen,
															 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
															 Frame_Shift     = talSites$frameShift,
															 Tool_Type       = talSites$Target,
															 Strand          = talSites$Orientation,
															 Exon_ID         = talSites$Exon_Num,
															 DSB_Location    = talSites$CutIndex,
															 Microhomology   = talSites$topMH,
															 PreMA_Sequence  = talSites$topDel,
															 Context         = talSites$seq,
															 stringsAsFactors = FALSE)
	}
	
	# Return frame
	if(pamFlag && talFlag){
		return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(pamFlag && !talFlag){
		return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(!pamFlag && talFlag){
		return(list(talFormFrame, 0, 0, 0, 0, 0))
	} else {
		return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	}
}


#' calculateMENTHUEnsembl
#'
#' @param pamList 
#' @param talenList 
#' @param gbFlag 
#' @param gbhFlag 
#' @param genbankInfo 
#' @param threshold 
#' @param firstExon 
#' @param exonTargetType 
#' @param exonStuff 
#' @param progress 
#'
#' @return
#' @export
#'
#' @examples

calculateMENTHUEnsembl <- function(pamList, cutDistList, ohList, wiggle = TRUE, wiggleRoom = 39, talenList, ensemblInfo, exonStuff, silent){
	require(plyr)

	# Update progress
	if(!silent) {
		print("Processing Ensembl sites...")
	}
	
	# Generate a subset of exons
	exonSubset <- ensemblInfo[which(as.numeric(ensemblInfo$rank) %in% as.numeric(exonStuff)), ]
	exonSubset <- exonSubset[order(as.numeric(exonSubset$rank)), ]
	
	# Get the exon sequences
	exonSeq  <- exonSubset$sequence
	
	# Create a data frame to hold information about the exon sequence
	exonDF <- data.frame(Exon_Num         = as.numeric(exonSubset$rank),
											 absStart         = as.numeric(exonSubset$contextStart),
											 absEnd           = as.numeric(exonSubset$contextEnd),
											 exonStart        = 0 + as.numeric(exonSubset$exp5),
											 exonEnd          = nchar(exonSubset$sequence) - as.numeric(exonSubset$exp3),
											 seq              = exonSeq,
											 stringsAsFactors = FALSE)
	
	# If the user is using Cas:
	if(length(pamList) > 0){
		
		if(!silent) {
			print("Scanning for CRISPR target sites...")
		}
		
		# If there is exon information, use it to correct indexing, otherwise, exonStarts is NULL
		pamSites <- pamScan(pamList, 
												cutDistList, 
												ohList,
												exonSeq, 
												exonList   = exonDF$Exon_Num, 
												exonStarts = NULL, 
												findCut    = TRUE, 
												type       = "cas9", 
												wiggle     = wiggle, 
												wiggleRoom = wiggleRoom)
		
		siteCount <- nrow(pamSites)
		
		# Set pamFlag TRUE - PAMs are used
		pamFlag <- TRUE
		
		# Update progress bar
		if(!silent) {
			print("Pre-processing CRISPR target sites...")
			
		}

		pamSites        <- unique(suppressMessages(plyr::join(pamSites, exonDF, by = 'Exon_Num')))
		
		# THIS NEEDS FIXING
		pamSites$absCut <- pamSites$CutIndex + pamSites$absStart
		
		# Drop target sites where the cut site is not within the exon boundaries
		keep     <- sapply(1:nrow(pamSites), function(x) pamSites$CutIndex[x] %in% seq(from = pamSites$exonStart[x], to = pamSites$exonEnd[x], by = 1))
		pamSites <- pamSites[keep, ]
		
		# Identify sites with enough sequence context to do calculations
		pamSites$contextCondition[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(pamSites$seq)))  ] <- TRUE
		pamSites      <- pamSites[intersect(which(pamSites$CutIndex >= 40), which((pamSites$CutIndex + 40) <= nchar(pamSites$seq))), ]
		
		siteCountC <- nrow(pamSites)
		
		# Get the sequence context surrounding the cut site
		forguide <- unlist(lapply(1:nrow(pamSites), function (x) substr(pamSites$seq[x], 
																																		pamSites$CutIndex[x] - 39, 
																																		pamSites$CutIndex[x] + 40)))
		context  <- unlist(lapply(1:nrow(pamSites), function (x){ paste0(substr(pamSites$seq[x], 
																																						pamSites$CutIndex[x] - 39, 
																																						pamSites$CutIndex[x]), 
																																		 substr(pamSites$seq[x], 
																																		 			 pamSites$CutIndex[x] + pamSites$ohLen[x] + 1, 
																																		 			 pamSites$CutIndex[x] + pamSites$ohLen[x] + 40))}
		))
		# Set the sequence context in the frame
		pamSites$seq   <- context
		pamSites$guide <- forguide
		
	} else {
		# If the user is NOT using Cas, set pamFlag to FALSE
		pamSites   <- 0
		pamFlag    <- FALSE
		
		siteCount  <- 0
		siteCountC <- 0
	}
	
	# If there are TALEN inputs
	if(length(talenList) > 0){
		# Set the range flag to true
		rFlag <- TRUE
		
		# Set all exon starts to the exon starts in the input frame
		# Submit talen info to talPal
		# If there are exon inputs
		talSites <- talPal(exonSeq,
											 findCut    = TRUE,
											 wiggle     = TRUE,
											 wiggleRoom = 39,
											 range      = rFlag, 
											 tL         = talenList[1], 
											 tR         = talenList[3], 
											 spa        = talenList[2],
											 exonList   = exonDF$Exon_Num,
											 exonStarts = NULL)

		
		# Update progress bar
		if(!silent) {
			print("Pre-processing TALEN target sites...")
		}
		talSites <- unique(suppressMessages(plyr::join(talSites, exonDF, by = 'Exon_Num')))
		
		# Drop target sites where the cut site is not within the exon boundaries
		keepT     <- sapply(1:nrow(talSites), function(x){
			talSites$CutIndex[x] %in% seq(from = talSites$exonStart[x], to = talSites$exonEnd[x], by = 1)
		})
		talSites  <- talSites[keepT, ]
		
		# Identify sites with enough sequence context to do calculations
		talSites$contextCondition[intersect(which(talSites$CutIndex >= 40), which((talSites$CutIndex + 40) <= (talSites$exonEnd + 39)))] <- TRUE
		talSites <- talSites[which(talSites$contextCondition), ]
		
		# Get the sequence context surrounding the cut site
		contextT <- unlist(lapply(1:nrow(talSites), function (x) substr(exonSeq[as.numeric(talSites$Exon_Num[x])], 
																																		talSites$CutIndex[x] - 39, 
																																		talSites$CutIndex[x] + 40)))
		
		
		# Set the sequence context in the frame
		talSites$seq = contextT
		
	} else {
		# If TALENs are not used, set talSites list to empty
		talSites <- 0
	}
	
	# Create data frame to hold results
	menthuFrame <- data.frame(Target_Sequence  = as.character(), 
														MENTHU_Score     = as.numeric(), 
														Frame_Shift      = as.character(), 
														Tool_Type        = as.character(), 
														Strand           = as.character(), 
														Exon_ID          = as.numeric(), 
														DSB_Location     = as.integer(),
														stringsAsFactors = FALSE)
	
	if(pamFlag){
		if(!silent) {
			print("Calculating MENTHUv2.0 scores for CRISPR sites...")
		}
		
		# Function to calculate MENTHU v2.0 score
		menthuFrameFunc <- function(x){
			return(calculateMenthu2(as.character(x), cutSite = 40, weight = 20, maxdbm = 5))
		}
		
		# Get the menthu scores
		menthuFrame <- as.data.frame(matrix(unlist(sapply(context, menthuFrameFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		# Update progress
		if(!silent) {
			print("Formatting CRISPR site results...")
		}
		
		row.names(menthuFrame)  <- c()
		colnames(menthuFrame)   <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		menthuFrame$menthuScore <- as.numeric(menthuFrame$menthuScore)
		
		# Get the critOne success count
		critOne <- nrow(menthuFrame[which(menthuFrame$frameShift != "NA"), ])
		
		# Get the critTwo success count
		critTwo <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5), ])
		
		# Get both success count
		critBoth <- nrow(menthuFrame[which(menthuFrame$menthuScore >= 1.5 & menthuFrame$frameShift != "NA"), ])
		
		pamSites <- unique(suppressMessages(plyr::join(pamSites, menthuFrame)))
		
		# Drop 0s
		pamSites <- pamSites[which(pamSites$menthuScore > 0), ]
		
		baseCrispr <- sapply(1:nrow(pamSites), 
												 function(x) if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] - 19, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 20)
												 	}
												 })
		
		# Generate the PAM sequence
		pam        <- sapply(1:nrow(pamSites), 
												 function(x) (if(pamSites$Orientation[x] == "forward"){
												 	if(pamSites$CutDist[x] < 0){
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))      
												 	} else {
												 		substr(pamSites$guide[x], 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])      
												 	}
												 } else {
												 	if(pamSites$CutDist[x] < 0){
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1, 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + nchar(pamSites$Target[x]))
												 	} else {
												 		substr(reverseComplement(pamSites$guide[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x] + 1 - nchar(pamSites$Target[x]), 
												 					 (nchar(pamSites$guide[x]) / 2) - pamSites$CutDist[x])
												 	}
												 }))
		
		# Get the CRISPR target sequence required to target this site, and bold the PAM
		crispr <- sapply(1:length(baseCrispr), function(x) {
			if(pamSites$CutDist[x] < 0) {
				paste0(baseCrispr[x], pam[x])
			} else{
				paste0(pam[x], baseCrispr[x])
			}})
		
		# Create data frame of the current results
		pamFormFrame  <- data.frame(Target_Sequence  = crispr, 
																MENTHU_Score     = round(as.numeric(pamSites$menthuScore), digits = 2), 
																Frame_Shift      = pamSites$frameShift,
																Tool_Type        = pamSites$Target, 
																Strand           = pamSites$Orientation, 
																Exon_ID          = pamSites$Exon_Num, 
																DSB_Location     = pamSites$CutIndex,
																Microhomology    = pamSites$topMH,
																PreMA_Sequence   = pamSites$topDel,
																Context          = pamSites$seq,
																stringsAsFactors = FALSE)
	}
	
	if(talFlag){
		# Update progress bar
		if(!silent) {
			print("Calculating MENTHUv2.0 scores for TALEN sites...")
		}
		
		menthuFrameTFunc <- function(x){

			# Return calculation
			return(calculateMenthu2(as.character(x), weight = 20, maxdbm = 5))
		}
		
		# Calculate slope competition on all the context
		menthuFrameT <- as.data.frame(matrix(unlist(sapply(contextT, menthuFrameTFunc)), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
		
		# Update progress bar
		if(!silent){
			print("Formatting TALEN site results...")
			
		}
		
		# Clean up the resulting data frame
		colnames(menthuFrameT) <- c("seq", "menthuScore", "frameShift", "topDel", "topMH")
		rownames(menthuFrameT) <- c()
		
		# Merge slope frame and pamSites
		talSites <- unique(suppressMessages(plyr::join(talSites, menthuFrameT)))
 		
		talenGenFunc <- function(talRow){
			dim  <- unlist(strsplit(talRow$Target, "/"))
			arm1 <- as.numeric(dim[1])
			spa1 <- as.numeric(dim[2])
			arm2 <- as.numeric(dim[3])
			
			armL <- substr(talRow$seq, start = 40 - (spa1 / 2) - arm1, stop = 40 - (spa1 / 2) - 1)
			spac <- substr(talRow$seq, start = 40 - (spa1 / 2),        stop = 40 + (spa1 / 2) - 1)
			armR <- substr(talRow$seq, start = 40 + (spa1 / 2),        stop = 40 + (spa1 / 2) - 1 + arm2)
			
			return(paste0(armL, spac, armR))
		}
		
		# Generate the crispr target sequences
		talSites$talen <- sapply(1:nrow(talSites), function(x) talenGenFunc(talSites[x, ]))
		
		# Clean the data frame
		row.names(talSites) <- c()
		
		# Create a new data frame of the results
		talFormFrame <- data.frame(Target_Sequence = talSites$talen,
															 MENTHU_Score    = round(abs(as.numeric(talSites$menthuScore)), digits = 2),
															 Frame_Shift     = talSites$frameShift,
															 Tool_Type       = talSites$Target,
															 Strand          = talSites$Orientation,
															 Exon_ID         = talSites$Exon_Num,
															 DSB_Location    = talSites$CutIndex,
															 Microhomology   = talSites$topMH,
															 PreMA_Sequence    = talSites$topDel,
															 Context         = talSites$seq,
															 stringsAsFactors = FALSE)
		
		# Get the critOne success count
		critOneT <- nrow(talFormFrame[which(talFormFrame$frameShift != "NA"), ])
		
		# Get the critTwo success count
		critTwoT <- nrow(talFormFrame[which(talFormFrame$menthuScore >= 1.5), ])
		
		# Get both success count
		critBothT <- nrow(talFormFrame[which(talFormFrame$menthuScore >= 1.5 & talFormFrame$frameShift != "NA"), ])
		
	}
	
	if(pamFlag && talFlag){
		critOne   <- critOne  + critOneT
		critTwo   <- critTwo  + critTwoT
		critBoth  <- critBoth + critBothT
		
	} else if(talFlag){
		critOne   <- critOneT
		critTwo   <- critTwoT
		critBoth  <- critBothT
	}

	if(pamFlag && talFlag){
		return(list(rbind(pamFormFrame, talFormFrame), siteCount, siteCountC, critOne, critTwo, critBoth))
	} else if(pamFlag && !talFlag){
		return(list(pamFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	} else {
		return(list(talFormFrame, siteCount, siteCountC, critOne, critTwo, critBoth))
	}
}


#' convertToNumeric
#'
#' This function takes an input string of comma- or whitespace-separated numbers and converts to a numeric vector
#' @param characterStringOfNumbers 
#'
#' @return
#' @export
#'
#' @examples

convertToNumeric <- function(characterStringOfNumbers){
	# Split the string by ','
	splittedString <- unlist(strsplit(characterStringOfNumbers, "[\\, |\\,| ]+"))
	
	# If there is more than one object in the list
	if(length(splittedString) > 1){
		# Find instances where a range was specified
		rangeLocs <- grep("-", splittedString, fixed = TRUE)
		
		# If more than one range
		if(length(rangeLocs) > 1){
			rangeEnds <- strsplit(splittedString[rangeLocs], "-")
			numberList <- splittedString[-rangeLocs]
			
			# Generate a sequence in the range
			for(i in 1:length(rangeEnds)){
				numberList <- c(numberList, seq(from = rangeEnds[[i]][1], to = rangeEnds[[i]][2]))
			}
		} else {
			numberList <- splittedString
		}
		
		numberList <- as.numeric(unlist(numberList))
		numberList <- sort(numberList)
		
	} else {
		numberList <- as.numeric(splittedString)
	}
	
	return(numberList)
}

####Exon Handler for Custom Exon Input####
exonHandler <- function(exonRHandsonTable){
	exonTable <- hot_to_r(exonRHandsonTable)
	exonDF <- exonTable[apply(exonTable, MARGIN = 1, function(x) any(x > 0)),]
	
	return(exonDF)
}


#' window
#'  
#' This function gets the necessary sequence context to do MENTHU calculations
#' 
#' @examples 
#' 
#' @export 
#' 

window <- function(sequence, position, winSize = 80) {
	
	if(position < (winSize / 2)){
		return("")
		
	} else {
		return(substring(sequence, ((position - ((winSize + (winSize %% 2)) / 2) - 1)), (position + ((winSize + (winSize %% 2)) / 2))))
	}
	
}

# Filter Results
filterResults <- function(results, opT7, opThresh){
	rResults <- results
	
	if(opT7){
		rResults <- rResults[which(!grepl("/", rResults$Tool_Type, fixed = TRUE)), ]
		rResults <- rResults[which(unlist(lapply(1:nrow(rResults), 
																						 function(x) grepl("((^[G]{1,2})|(^[ACGT][G]))", rResults$Target_Sequence[x], perl=TRUE, ignore.case=TRUE)))), ]
	}
	
	if(opThresh){
		rResults <- rResults[which(rResults$MENTHU_Score >= 1.5), ]
	}
	
	if(nrow(rResults) < 1){
		return(list(FALSE, ""))
		
	} else {
		return(list(TRUE, rResults))
		
	}
}

#' reverse
#'
#' This function takes a string as input and reverses the order of the string
#'
#' @param seq A string to reverse
#'
#' @return revSeq The seq string in reverse
#'
#' @examples
#' reverse("123456")
#'
#'
#' @export

reverse <- function(seq){
	UseMethod("reverse", seq)
}

reverse.default <- function(seq){
	stop("Error: Cannot reverse objects that are not character strings or integers. Please check input sequence.")
}


reverse.character <- function(seq){
	revSeq <- seq
	
	for(i in 1:nchar(seq)){
		curL <- substr(seq, i, i)
		substr(revSeq, (nchar(seq) + 1 - i), (nchar(seq) + 1 - i)) <- curL
	}
	
	return(revSeq)
}

reverse.numeric <- function(seq){
	charSeq <- as.character(seq)
	revSeq <- charSeq
	revSeq <- as.numeric(reverse.character(charSeq))
	return(revSeq)
}

#' complement
#'
#' This function takes a DNA or RNA sequence as input (along with a parameter specifying the type of sequence) 
#' and outputs the complement of the input sequence. E.g., "ATTG" will return "TAAC" if type = "DNA" and "UAAC" if type = "RNA"
#'
#' @param seq A DNA or RNA sequence from which to generate a complement string
#' @param type (Now defunct; currently, all IUPAC 1-letter nucleotide codes are supported.) Default is "DNA"; a DNA sequence can only contain 
#' "A", "C", "G", or "T" for the purposes of complement(). The other option is "RNA"; an RNA sequence can only contain 
#' "A", "C", "G", or "U" for the purposes of complement().
#'
#' @return compSeq The complement of the input sequence
#' @export
#'
#' @examples
#' seq <- "AAAATGGCGAAG"
#' type <- "DNA"
#' complement(seq, type)

complement <- function(seq, type){
	UseMethod("complement", seq)
}

complement.default <- function(seq, type){
	#Prevent attempts to complement objects that are not sequences
	stop("Error: Cannot complement a non-character vector object. Please check input sequence.")
}

complement.character <- function(seq, type){
	compSeq <- seq
	
	fromVal <- c("A", "C", "G", "T", "a", "c", "g", "t", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v", "n")
	toVal   <- c("T", "G", "C", "A", "t", "g", "c", "a", "Y", "R", "S", "W", "M", "K", "V", "H", "D", "B", "N", "y", "r", "s", "w", "m", "k", "v", "h", "d", "b", "n")
	
	
	compSeq <- plyr::mapvalues(unlist(strsplit(compSeq, split = "")),
														 from         = fromVal,
														 to           = toVal,
														 warn_missing = FALSE)
	compSeq <- paste(compSeq, collapse = "")
	return(compSeq)
}

complement.list <- function(seq, type){
	retList <- list()
	for(i in seq){
		retList <- c(retList, complement(i, type))
	}
	return(retList)
}

#' reverseComplement
#'
#' This function takes a DNA or RNA sequence as input and outputs the reverse complement of the sequence.
#'
#' @param seq A character vector from which to generate a reverse complement.
#' @param type This is deprecated, as 'complement' now deals with all IUPAC 1-letter codes. 
#' Default is "DNA"; allowed characters are "A", "C", "G", and "T" (case insensitive). 
#' Other option is "RNA"; allowed characters are "A", "C", "G", and "U" (case insensitive.)
#'
#' @return seqRevComp The reverse complement of the input sequence
#' @export
#'
#' @examples
#' dnaSeq <- "AATGCC"
#' reverseComplement(dnaSeq)
#' rnaSeq <- "UUAGCC"
#' reverseComplement(rnaSeq, type = "RNA")

reverseComplement <- function(seq, type = "DNA"){
	UseMethod("reverseComplement", seq)
}

reverseComplement.default <- function(seq, type = "DNA"){
	stop("Error: Input sequence is not a character vector. Please check input sequence.")
}

reverseComplement.character <- function(seq, type = "DNA"){
	#Reverse the sequence
	seqRev <- reverse(seq)
	#Get the complement of the reversed sequence
	seqRevComp <- complement(seqRev, type = "DNA")
	return(seqRevComp)
}

reverseComplement.list <- function(seq, type = "DNA"){
	retList <- list()
	for(i in seq){
		retList <- c(retList, reverseComplement.character(i))
	}
	return(unlist(retList))
}

#'
#' This function removes all white space from character vectors. If it is passed a data frame, it will remove all white space from all columns with character data types.
#'
#' @param wsco An object to be stripped of white space
#'
#' @return stripped The object, stripped of white space
#' @export
#'
#' @examples
#' stripWhiteSpace(c("red", " gre en ", "  ", " blue "))
#' V1 <- c("  1", " 2 ", "    3")
#' V2 <- c(1, 2, 3)
#' dummyDF <- data.frame(V1, V2, stringsAsFactors=FALSE)
#' dummyDF
#' stripWhiteSpace(dummyDF)

stripWhiteSpace <- function(wsco) UseMethod("stripWhiteSpace")

#Default method; will not run the method unless it is passed a data frame with at least one column of characters or a character vector
stripWhiteSpace.default <- function(wsco){
	stop("Object is not a character vector")
}

#For handling character vectors
stripWhiteSpace.character <- function(wsco){
	stripped <- gsub('\\s+', '', wsco)
	return(stripped)
}

#For handling data frames with at least one character column type
stripWhiteSpace.data.frame <- function(wsco){
	count <- 0
	dDF <- wsco
	for(i in 1:ncol(wsco)){
		if(class(wsco[,i])=="character"){
			dDF[,i] <- stripWhiteSpace(wsco[,i])
			count <- count + 1
		}
	}
	#Determine if any of the columns are character vectors; stop execution if they are not
	if(count == 0){
		stop("Data frame has no character columns")
	}
	return(dDF)
}

