##### Validation Functions ####
# Check for a valid GenBank function
validGenBankId <- function(id){
	errMess <- c()
	
	if(id != ""){
		#Let RefSeq RNA accessions through
		if(       stringr::str_detect(id, regex("^(NM|NR|XM|XR)_[0-9]{6}",       ignore_case = TRUE))){
			
			#Catch RefSeq protein accesssions
		} else if(stringr::str_detect(id, regex("^(AP|NP|YP|XP|WP)_[0-9]{6}",    ignore_case = TRUE))){
			errMess <- c(errMess, "Error: RefSeq protein accession detected; please use a NUCLEOTIDE ID.")
			
			# Catch ginormous genomic region files; disabled for running locally
		} else if(stringr::str_detect(id, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}", ignore_case = TRUE))){
			
			# Catch GenBank protein accessions
		} else if(stringr::str_detect(id, regex("^[A-Z]{3}[0-9]{5}",             ignore_case = TRUE))){
			errMess <- c(errMess, "Error: GenBank protein accession detected; please use a NUCLEOTIDE ID.")
			
			# Catch anything else not conforming to input types
		} else {
			if(any(stringr::str_detect(id, regex("^[a-zA-Z]{2}[0-9]{6}",                      ignore_case = TRUE)), 
						 stringr::str_detect(id, regex("^[a-zA-Z]{1}[0-9]{5}",                  ignore_case = TRUE)),
						 stringr::str_detect(id, regex("^(AC|NC|NG|NT|NW|NZ)_[0-9]{6}\\.[0-9]", ignore_case = TRUE)),
						 stringr::str_detect(id, regex("^[a-zA-Z]{4}[0-9]{8,10}",               ignore_case = TRUE)))) {
				
			} else {
				errMess <- c(errMess, "Error: ID doesn't match a known GenBank or RefSeq NUCLEOTIDE ID format.")
			}
		}
	} 
	
	if(length(errMess) > 0) {
		return(c(FALSE, errMess))
	} else {
		return(c(TRUE, ""))
	}
}


# Check for a valid Ensembl function


ensemblIdSpecies <- function(id, ensIdList = ensIds, bool = FALSE){
	# Pull out the mouse matches
	if(substr(id, 1, 3) == "MGP"){
		prefix <- unlist(stringr::str_extract_all(id, "[Mm][Gg][Pp][_][A-Za-z0-9]+[_]"))
	} else if(substr(id, 1, 2) == "FB"){
		# Pull out the drosophila matches, which for some inexplicable reason have their own type formatting
		prefix <- "FB"
		
	} else {
		# Get the ID number and the preceding two letters
		suffix <- unlist(stringr::str_extract_all(id, "[A-Za-z][A-Za-z][0-9.]+"))
		suffixTwoLetter <- substr(suffix, 1, 2)
		suffixOneLetter <- substr(suffix, 2, 2)
		
		if(grepl("FM|GT|fm|gt|Fm|Gt|fM|gT", suffixTwoLetter, perl = TRUE)){
			prefix <- gsub(suffix, "", id)
			
		} else if(grepl("E|e|G|g|P|p|R|r|T|t", suffixOneLetter, perl = TRUE)){
			prefix <- gsub(unlist(stringr::str_extract_all(id, "[A-Za-z][0-9.]+")), "", id)
			
		} else {
			
			if(bool){
				return(c(FALSE, "Error: Input ID does not match a recognized Ensembl ID format"))
				
			} else {
				return("")
				
			}
		}
	}
	
	# Match Ensembl id to list of Ensembl Ids
	match <- prefix %in% ensIdList$Id
	
	if(match){
		species <- ensIdList[which(ensIdList$Id == prefix), ]
		
		if(bool){
			return(c(TRUE, ""))
			
		} else {
			return(species)
			
		}
	} else {
		if(bool){
			return(c(FALSE, "Error: Input ID prefix does not match any known Ensembl prefixes."))
			
		} else {
			return("unknown")
		}
	}
}


###### Validate DNA sequence input #######
validDnaSeq <- function(seq){
	errMess <- c()
	if(grepl("[^ACGTacgt]", seq, perl = TRUE)) {
		errMess <- c(errMess, "Error: The DNA sequence contains non-standard nucleotides; allowed nucleotides are A, C, G, and T")
	}
	
	if(nchar(seq) < 80 ) {
		errMess <- c(errMess, "Error: MENTHU requires a DNA sequence of >= 80 nucleotides")
	}
	
	if (length(errMess) > 0 ){
		return(c(FALSE, errMess))
	} else {
		return(c(TRUE, ""))
	}
}

##### 

