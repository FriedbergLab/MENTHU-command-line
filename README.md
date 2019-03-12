# MENTHU-command-line v1.0.0-alpha

This repository contains scripts for running the MENTHU analysis at large scale, offline. It is compatible with Linux, Mac, Windows (using Cygwin), and Windows CMD.

This tool is intended for users who want to identify MENTHU scores for large sections of DNA (e.g., chromosomes, rather than genes).

If you only intend to determine MENTHU scores on a small scale, we highly recommend that you:

Run MENTHU online through a web interface here: http://genesculpt.org/menthu/

Run MENTHU through a GUI locally using the instructions here: https://github.com/Dobbs-Lab/MENTHU#run-menthu-locally

## MENTHU-command-line Installation
MENTHU-command-line requires an R installation in order to work. 

If you already have R installed, then you can skip the R installation instructions and go straight to [Package Installation](https://github.com/Dobbs-Lab/MENTHU-command-line#package-installation); see Step 2 ([Add R to Path](https://github.com/Dobbs-Lab/MENTHU-command-line#.


### [1. Download and Install R](#1-download-and-install-r)

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


### [2. Add R to PATH](#add-r-to-path)
MENTHU-command-line requires that the path to the R executable/binary be present in the environment PATH variable.
You can check if R is installed in your PATH with the following command (on Linux, Mac, Windows CMD, and Windows Cygwin):

```
R --version
```
If R is in your path, you will receive a message that looks a lot like below:

```
R version 3.4.4 <2018-03-15> -- "Someone to Lean On"
Copyright <C> 2018 THe R Foundation for Statistical Computing
```

You can then go directly to [Package Installation](https://github.com/Dobbs-Lab/MENTHU-command-line#package-installation).

If you receive ```-bash: R: command not found"``` on Unix/Cygwin or ```"'R' is not recognized as an internal or external command, operable program or batch file```, then R is not in your PATH, and you should follow the instructions below.

Replace [text in brackets] with the specified information for your system.

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

### [3. Download MENTHU-command-line](#download-mcl)

To download MENTHU-command-line, you can either:

(1) Click the green "Clone or download" button in the upper right corner and choose "Download ZIP"

OR

(2) From a Linux, Mac, or Windows Cygwin command line, navigate to the location you want to keep the MENTHU-command-line script, and enter the following command:

```
git clone https://github.com/Dobbs-Lab/MENTHU-command-line.git
```

### [4. Package Installation](#package-installation)
MENTHU-command-line requires a few additional R packages in order to function. You may need administrator or root privileges to install these packages.

MENTHU-command-line includes a script, packageInstaller.R, to perform easy 1-step installation of the needed R packages.

Use the following command the FIRST TIME you run MENTHU-command-line:

```
Rscript packageInstaller.R
```

This script will check if you have the appropriate packages installed, and then install them. This step may take several minutes, but does not need to be run subsequent times you run MENTHU-command-line (unless MENTHU-command-line is updated).

### [5. Updating MENTHU-command-line](#updating)
We will update the MENTHU-command-line code as needed (and there are a few additional convenience features under development.) Check the changelog for update information.

To update from the command line, navigate into the directory containing your MENTHU-command-line scripts, and use the following command:

```
git pull origin master
```

Otherwise you will need to re-download the zip file as in [Download MENTU-command-line](#download-mcl).

You should also re-run the packageInstaller.R script, in case new packages are now required.


## Running MENTHU-command-line

MENTHU-command-line can be run from Unix-like command lines (Linux, Mac, Cygwin on Windows) or Windows CMD. However, due to end-of-line conversion differences between Unix and Dos, **you must use menthu.R on Linux, Mac, and Cygwin, and menthu-cmd.R in Windows CMD**. 

Check to make sure you're using the correct command for your system!

MENTHU-command-line can be run using the following syntax:

```
Rscript menthu.R [outFile] [CRISPR Option] [PAM Sequence] [Distance to DSB] [Overhang] [TALEN Option] [TALEN scheme] [Gen Input Type] [Gen Input] [Score Threshold] [T7 opt] [verbose] [validate]
```

"Rscript" tells the system to use ```R``` to execute ```menthu.R``` (or ```menthu-cmd.R``` for Windows CMD).

### Parameter Explanation

The parameters are explained below. Each parameter is delimited by a space. Parameter values should not have spaces; if you want to put spaces in the output file name, the name should be in quotes, e.g. "outputFile.csv". Parameter values (including strings) do not have to be in quotes, except for the output file name exception.

### Parameter:				Accepted Values:						Explanation:
**"outFile"**					character string	file name	The name of the file to output your results to. If using a fasta file with multiple sequences, multiple files will be created, using this as a prefix

**"CRISPR Option"**		T or F									Flags the system to use CRISPR nuclease processing. If this option is T (true), "TALEN Option" must be F (false)  

**"PAM Sequence"**		A PAM sequence							The PAM sequence for the CRISPR sequence. Ambiguous nucleotides are allowed. Using N will scan every possible cut site in the target sequence. This parameter must be present, but is not used, if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)

**"Distance to DSB"**	Integer											The distance from the PAM sequence to the DSB site. For DSBs upstream of a PAM, use a negative value (e.g., -3 for SpCa9); for downstream, use a positive value (e.g., 18 for Cas12a.) This parameter must be present, but is not used, if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)

**"Overhang"**				Integer >= 0								The length of 5' overhang produced by the nuclease (e.g., 5 for Cas12a). Use 0 for blunt-cutting nucleases. This parameter must be present, but is not used, if "CRISPR Option" is false (i.e., you can put a 0 or NA in this spot.)

**"TALEN Option"**		T or F									Flags the system to use TALEN processing. If this option is T (true), "CRISPR Option" must be F (false)

**"TALEN Scheme"**		15-18/14 or 16/15-18			The left arm length, spacer length, and right arm length to use when searching for TALEN locations. E.g., for a TALEN	with arms 15 nt long, with spacer 14 nt, use 15/14/15. TALEN arms can be 15-18 nt in length; the spacer should be 14 OR 16 nt	in length (15 is not allowed for the spacer) This parameter must be present, but is not used, if "TALEN Option" is false (i.e., you can put a 0 or NA in this spot.)

**"Gen Input Type"**	gb ens seq file			Flags the system to get a GenBank/RefSeq ID (gb), Ensembl ID (ens), DNA sequence (seq), or to expect a FASTA file (file)	

**"Gen Input"**				See explanation							Provide the accession for GenBank/RefSeq/Ensembl inputs, file name for "file" option, or DNA sequence for "seq". If the file name has spaces in it, put this parameter in quotes.

**"Score Threshold"**	Positive number							Only output results with MENTHU score above this threshold. Default is 1.0. We recommend to only use sites with score >= 1.5

**"T7 opt"**					T or F									If T (true), only displays results where the gRNA is compatible with T7-cloning. 

**"Verbose"**					T or F								If T (true), outputs progress messags to the console.

**"validate"**				T or F								If T (true), checks the command line arguments to make sure they are all valid (this may take some time); if F, skip validation checks




**Example SpCas9, Ensembl:**
```
Rscript menthu.R EnsemblExample.csv T NGG -3 0 F NA ens ENSDART00000011520.8 1.5 F F F
```
**Example Cas12a:**
```
Rscript menthu.R EnsemblCas12aExample.csv T TTTN 18 5 F NA ens ENSDART00000011520.8 1.5 F F F
```
**Example TALEN (left arm: 15 nts, spacer: 16 nts; right arm 18 nts):**
```
Rscript menthu.R GenBankTalenExample.csv F NA 0 0 T 15/16/18 gb AY214391.1 1.5 F F F
```
