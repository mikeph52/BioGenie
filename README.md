# BioGenie
BioGenie is a complete bioinformatics command line tool for macOS, GNU Linux and MS Windows, written in C++.

> [!IMPORTANT]
> _BioGenie is now available on Conda. [Visit the conda page](https://anaconda.org/channels/mikeph52/packages/biogenie/overview)._

<img width="912" height="740" alt="Image" src="https://github.com/user-attachments/assets/3a5cc84b-1197-47ec-9638-068abe844c39" />

[![Anaconda-Server Badge](https://anaconda.org/mikeph52/biogenie/badges/version.svg)](https://anaconda.org/mikeph52/biogenie)
[![Anaconda-Server Badge](https://anaconda.org/mikeph52/biogenie/badges/platforms.svg)](https://anaconda.org/mikeph52/biogenie)
[![Anaconda-Server Badge](https://anaconda.org/mikeph52/biogenie/badges/downloads.svg)](https://anaconda.org/mikeph52/biogenie)
[![Anaconda-Server Badge](https://anaconda.org/mikeph52/biogenie/badges/latest_release_date.svg)](https://anaconda.org/mikeph52/biogenie)
[![Anaconda-Server Badge](https://anaconda.org/mikeph52/biogenie/badges/license.svg)](https://anaconda.org/mikeph52/biogenie)

It currently supports FASTA formats(.fasta, .faa, .fa) and FASTQ.
- To run the app, simply type:
```
biogenie <function> <filename>
```
- For example, to calculate GC percentage:
```
biogenie -gc example.fasta
```

## Documentation
BioGenie uses functions to execute different tools for different applications.
Read [Documentation](Documentation/documentation.md) for more information (*References included*).

### Pipelines:
- Nucleostats ---> "-nucleo". Statistics and metrics for DNA. Ideal for Primer design.
- Proteostats ---> "-prot". Statistics and protein structural properties.
- Assemblystats ---> "-asmbl". Statistics and metrics for genome assembly.
- Pairwise Global Sequence Alignment with Needleman-Wunsch ---> "-nw".
- Pairwise Local Sequence Alignment with Smith-Waterman ---> "-sw".

> [!TIP]
> _Assemblystats is in experimental phase. Always check your data first before evaluating your results._


### Single commands:
- Get the complement DNA sequence ---> "-c".
- Get the reverse complement DNA sequence ---> "-rc".
- Get the codon number ---> "-nc".
- Get the mRNA ---> "-t".
- Motif search ---> "-mf".
- Ambiguous bases statistics ---> "-amb".
- GC percentage calculation ---> "-gc".
- Generate the aminoacids(Protein chain) ---> "-p".
- Generate the Protein chain with color ---> "-pc".
- Color the protein chain from a FASTA file ---> "-pca".
- Separate different sequencies in a FASTA file ---> "-ss".
- Print the different sequence headers from a FASTA file ---> "-sh".
- Trim DNA sequence ---> "-tr".
- Calculate the Number of Base Pairs(bp) ---> "-bp".
- Get the purine/pyrimidine ratio --> "-pp".
- Calculate melting temperature (Tm) of DNA sequences using the Wallace Rule(only valid for oligos <20bp) ---> "-mt1".
- Calculate melting temperature (Tm) of DNA sequences using the SantaLucia 1998 nearest-neighbor method ---> "-mt2".
- Calculate the Isoelectric Point of a protein ---> "-pi".
- Calculate the molecular weight of a protein(kDa) ---> "-mw".
- Calculate the Extinction Coefficient of a protein ---> "-ec".
- Coloured cDNA sequence ---> "-cc".
- Coloured DNA sequence ---> "-sc".
- Get the Open Reading Frame(ORF) ---> "-orf".
- Generate cDNA sequence FASTA ---> "-cw".
- Generate Reverse cDNA sequence FASTA ---> "-rcw".
- Generate mRNA sequence FASTA ---> "-tw".
- Calculate Codon Usage Bias(CUB) ---> "-cub".
- Export Codon Usage Bias(CUB) to CSV file ---> "-wcub".
- Calculate Hydrogen Bonds of dsDNA ---> "-hb".
- Caclulate the Genome coverage(x) ---> "-cx".
- Calculate FASTQ quality control ---> "-fqc".
- K-mer Frequency Analysis ---> "-kmer".

More functions will be added in the future.

> [!NOTE]
> _If you have any suggestions for new features or a bug encountered, create an Issue or send me a message at: mikeph526@outlook.com. I'm happy to help._

## Testing
Various tests have been performed to ensure scientific accuracy and the optimal computing performance. Everything involving unit testing can be found in this [repository](https://github.com/mikeph52/biogenie-test).

### Examples
In order to test BioGenie functions, example files are provided in the [sample fasta](sample_fasta) folder. In the table bellow, are all the necessary information about the said files.

|Name|Description|NCBI Acc. Number|
|----|-----------|----------------|
|pst_dc3k.fa|Genome of _P.syringae DC3000_|ASM780v1|
|hopB1.fasta|Sequence of hopB1 gene|NC_004578.1|
|pgp_prot.faa|Protein sequence of Human PGP|NP_001035830.1|
|primer.fasta|A selection of well-known primers| _Various_|
|hemoglobin_pdec.fasta|_P.decipiens_ gene for hemoglobin|Z11681.1 |
|test.fastq| A fictional sequence to test FASTQ support| _Fictional_|
|pwa.fasta| A test file for PGSA testing|_Fictional_|
|pwa2.fasta| A test file for PLSA testing|_Fictional_|

## Installation
### GNU Linux
- Download BioGenie from Releases, or with wget:
```
wget https://github.com/mikeph52/BioGenie/releases/download/v.1.3.0/biogenie_linux_1.3.0
```
- Run "chmod +x" first.
```
sudo chmod +x biogenie_linux_1.3.0
```
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_linux_1.3.0 /usr/local/bin/biogenie
```
- If you need to build from source:
```
git clone https://github.com/mikeph52/BioGenie.git
clang++ main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### macOS
It supports Intel and Apple silicon processors.
- Download BioGenie from Releases, or with curl:
```
curl -l https://github.com/mikeph52/BioGenie/releases/download/v.1.3.0/biogenie_macos_1.3.0
```
- Run "chmod +x" first.
```
sudo chmod +x biogenie_macos_1.3.0
```
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_macos_1.3.0 /usr/local/bin/biogenie
```
- If you need to build from source(probably not):
```
git clone https://github.com/mikeph52/BioGenie.git
clang++ -std=c++17 main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### MS Windows
> [!IMPORTANT]
> _Windows is back(sort of). From version v.0.28.0, windows is supported. Still, this version is not going to be prioritized._

It requires Windows 10 or later.
- Download BioGenie from Releases or run:
```batch
wget.exe https://github.com/mikeph52/BioGenie/releases/download/v.1.3.0/biogenie_windows_1.3.0.exe
```
- Open **Start** and search for `Environment Variables`. 
- Click `Edit the system environment variables`. 
- In the System Properties window, click `Environment Variables`. 
- Under `User variables` (for current user) or `System variables` (for all users), select `Path` and click `Edit`. 
- Click `New` and paste the folder path where your _biogenie_ `.exe` lives (for example `C:\Bioinfo\CLI`). 
- Click `OK` on all open dialogs to save and close them. 
- Close and reopen any Command Prompt/PowerShell windows so they pick up the updated PATH. 

~There's also a GUI Version Available on Alpha testing~.
> [!CAUTION]
> This is an Alpha testing version. It is not functional. Not for scientific use. I don't think it is possible to make and maintain a gui version for windows. The code is complicated already. It's time to move on.

### Instal with conda
From version 1.1.0 BioGenie is available on conda, [visit the conda page](https://anaconda.org/channels/mikeph52/packages/biogenie/overview).

**Supported architectures:**
- MacOS: Intel x86-64, Apple Silicon
- Linux: x86-64, aarch64
- Windows: x64

- Run this command:
```bash
conda install mikeph52::biogenie
```
- It is highly advisable to create a conda environment first, even though BioGenie doesn't have any depedencies:
```bash
conda create -n biogenie
```
_Soon it will uploaded on bioconda._

## Changelog:
### 1.3.0 (April 20th, 2026)
- Conda support for multiple platforms.
- K-mer Frequency Analysis function added.
- Small formating bugs fixed.

### 1.2.2 (March 24th, 2026)
- Issue (https://github.com/mikeph52/BioGenie/issues/50) fixed.

### 1.2.0 (March 24th, 2026)
- Quality control for FASTQ files added.
- Conda support added.
- Minor bugs fixed.

### 1.1.0 (March 6th, 2026)
- Genome Coverage(x) added.
- Genome Coverage(x) added in "-asmbl" pipeline.
- FASTQ support added.
- FASTQ to FASTA function added.
- Major file-handling improvements made.
- Documentation enhanced.

### 1.0.0 (March 3rd, 2026)
- Issue (https://github.com/mikeph52/BioGenie/issues/42) fixed.
- Documentation enhanced. 
- Minor bugs and formatting issues fixed.

_For more info check out the [changelog](Documentation/changelog.md)._

