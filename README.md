# MENTHU-command-line

This repository contains scripts for running the MENTHU analysis at large scale, offline.

Detailed instructions forthcoming.

#### Usage ####
Rscript menthu.R ["outFile"] ["CRISPR Option"] ["PAM Sequence"] [Distance to DSB] [Overhang] ["TALEN Option"] ["TALEN scheme"] ["Gen Input Type"] ["Gen Input"] [Score Threshold] ["T7 opt"] [silent] [validate]

Example SpCas9, Ensembl:
Rscript menthu.R 'EnsemblExample.csv' 'T' 'NGG' -3 0 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'

Example Cas12a:
Rscript menthu.R 'EnsemblCas12aExample.csv' 'T' 'TTTN' 18 5 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'

Example TALEN (left arm: 15 nts, spacer: 16 nts; right arm 18 nts):
Rscript menthu.R 'GenBankTalenExample.csv' 'F' NA 0 0 'T' '15/16/18' 'gb' 'AY214391.1' 1.5 'F' 'F' 'F'

Parameter Explanation
Parameter:				Accepted Values:						Explanation:
"outFile"					character string	file name	The name of the file to output your results to. If using a fasta file with multiple sequences, multiple files will be created, using this as a prefix
"CRISPR Option"		"T" or "F"									Flags the system to use CRISPR nuclease processing. If this option is "T" (true), "TALEN Option" must be "F" (false)  
"PAM Sequence"		A PAM sequence							The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are allowed. Using "N" will scan every possible cut site in the target sequence
"Distance to DSB"	Integer											The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, use a negative value (e.g., "-3" for SpCa9); for downstream, use a positive value (e.g., "18" for Cas12a.)
"Overhang"				Integer >= 0								The length of 5' overhang produced by the nuclease (e.g., "5" for Cas12a). Use "0" for blunt-cutting nucleases. 
"TALEN Option"		"T" or "F"									Flags the system to use TALEN processing. If this option is "T" (true), "CRISPR Option" must be "F" (false)
"TALEN Scheme"		"15-18/14 or 16/15-18"			The left arm length, spacer length, and right arm length to use when searching for TALEN locations. E.g., for a TALEN	with arms 15 nt long, with spacer 14 nt, use "15/14/15". TALEN arms can be 15-18 nt in length; the spacer should be 14 OR 16 nt	in length (15 is not allowed for the spacer)
"Gen Input Type"	"gb" "ens" "seq" "file"			Flags the system to get a GenBank/RefSeq ID ("gb"), Ensembl ID "ens", DNA sequence ("seq"), or to expect a FASTA file ("file")	
"Gen Input"				See explanation							Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for "file" option, or DNA sequence for "seq"
"Score Threshold"	Positive number							Only output results with MENTHU score above this threshold. Default is 1.0. Recommended only use sites with score >= 1.5
"T7 opt"					"T" or "F"									If "T" (true), only displays results where the gRNA is compatible with T7-cloning. 
"silent"					"T" or "F"									If "T" (true), does not output progress messages to the console. 
"validate"				"T" or "F"									If "T" (true), checks the command line arguments to make sure they are all valid (this may take some time); if "F", skip validation checks

