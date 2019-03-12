# MENTHU-command-line

This repository contains scripts for running the MENTHU analysis at large scale, offline.

You can run MENTHU online through a web interface here: http://genesculpt.org/menthu/

You can run MENTHU through a GUI locally using the instructions here: https://github.com/Dobbs-Lab/MENTHU#run-menthu-locally

## Installation
MENTHU-command-line requires an R installation in order to work. 

If you already have R installed, then you simply need to download and unzip the files in this repsository, change directory to the unzipped files, and

If you do not have R installed: 

Download R for your appropriate operating system:

Windows: https://mirror.las.iastate.edu/CRAN/bin/windows/

You should select the "base" option, or click "install R for the first time".
 

Mac OS: https://mirror.las.iastate.edu/CRAN/bin/macosx/

Scroll down to the "Files" section, and find the R pkg file that lists your operating system (El Capitan, Mavericks, Snow Leopard, etc). Select the R-3.x.x.pkg file corresponding to your system - pay special attention to the "Important" section under R-3.4.3.pkg if you have "El Capitan"; you may want to consider using R-3.3.3.pkg if you don't want to install additional tools to support running R 3.4.3 on "El Capitan".


Linux/Unix: https://mirror.las.iastate.edu/CRAN/bin/linux/

Find your Unix distro, and folow the instructions in the directory.
 

Once you have downloaded the R installer, run it to install R. You may be required to enter administrator credentials; if you do not have these credentials, talk to your institution's IT department to have them install the software.


If you need additional help installing R, please check the installation instructions for your operating system:

Windows:    https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Windows

Mac OS:     https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS

Linux/Unix: https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes


### Adding R to PATH
MENTHU-command-line requires that the path to the R executable/binary be present in the environment PATH variable.
Instructions for doing this are below.
You should replace [text in brackets] with the specified information for your system.

Linux:
```
export PATH="[/path/to/R/]:$PATH"
```
Mac:
```
export PATH="[/path/to/R/]:$PATH"
```
Windows, CMD:
```
setx path "%path%;[\path\to\R\]"
```
Windows, Cygwin:
```
export PATH="[/path/to/R/]:$PATH"
```

## Usage 
Rscript menthu.R ["outFile"] ["CRISPR Option"] ["PAM Sequence"] [Distance to DSB] [Overhang] ["TALEN Option"] ["TALEN scheme"] ["Gen Input Type"] ["Gen Input"] [Score Threshold] ["T7 opt"] [silent] [validate]

**Example SpCas9, Ensembl:**
```
Rscript menthu.R 'EnsemblExample.csv' 'T' 'NGG' -3 0 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'
```
**Example Cas12a:**
```
Rscript menthu.R 'EnsemblCas12aExample.csv' 'T' 'TTTN' 18 5 'F' NA 'ens' 'ENSDART00000011520.8' 1.5 'F' 'F' 'F'
```
**Example TALEN (left arm: 15 nts, spacer: 16 nts; right arm 18 nts):**
```
Rscript menthu.R 'GenBankTalenExample.csv' 'F' NA 0 0 'T' '15/16/18' 'gb' 'AY214391.1' 1.5 'F' 'F' 'F'
```
### Parameter Explanation
### Parameter:				Accepted Values:						Explanation:
**"outFile"**					character string	file name	The name of the file to output your results to. If using a fasta file with multiple sequences, multiple files will be created, using this as a prefix

**"CRISPR Option"**		"T" or "F"									Flags the system to use CRISPR nuclease processing. If this option is "T" (true), "TALEN Option" must be "F" (false)  

**"PAM Sequence"**		A PAM sequence							The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are allowed. Using "N" will scan every possible cut site in the target sequence

**"Distance to DSB"**	Integer											The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, use a negative value (e.g., "-3" for SpCa9); for downstream, use a positive value (e.g., "18" for Cas12a.)

**"Overhang"**				Integer >= 0								The length of 5' overhang produced by the nuclease (e.g., "5" for Cas12a). Use "0" for blunt-cutting nucleases. 

**"TALEN Option"**		"T" or "F"									Flags the system to use TALEN processing. If this option is "T" (true), "CRISPR Option" must be "F" (false)

**"TALEN Scheme"**		"15-18/14 or 16/15-18"			The left arm length, spacer length, and right arm length to use when searching for TALEN locations. E.g., for a TALEN	with arms 15 nt long, with spacer 14 nt, use "15/14/15". TALEN arms can be 15-18 nt in length; the spacer should be 14 OR 16 nt	in length (15 is not allowed for the spacer)

**"Gen Input Type"**	"gb" "ens" "seq" "file"			Flags the system to get a GenBank/RefSeq ID ("gb"), Ensembl ID "ens", DNA sequence ("seq"), or to expect a FASTA file ("file")	

**"Gen Input"**				See explanation							Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for "file" option, or DNA sequence for "seq"

**"Score Threshold"**	Positive number							Only output results with MENTHU score above this threshold. Default is 1.0. Recommended only use sites with score >= 1.5

**"T7 opt"**					"T" or "F"									If "T" (true), only displays results where the gRNA is compatible with T7-cloning. 

**"silent"**					"T" or "F"									If "T" (true), does not output progress messages to the console. 

**"validate"**				"T" or "F"									If "T" (true), checks the command line arguments to make sure they are all valid (this may take some time); if "F", skip validation checks


