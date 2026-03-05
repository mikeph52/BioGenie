# BioGenie
BioGenie is a complete bioinformatics command line tool for macOS, GNU Linux and MS Windows, written in C++.

> [!IMPORTANT]
> _BioGenie 1.0.0 is finally out!!!_

<img width="1002" height="770" alt="Image" src="https://github.com/user-attachments/assets/4dfb77db-721d-42e0-ade9-73115f12057e" />


It currently supports fasta formats(.fasta, .faa, .fa).
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

More functions will be added in the future.

> [!NOTE]
> _If you have any suggestions for new features or a bug encountered, create an Issue or send me a message at: mikeph526@outlook.com. I'm happy to help._

## Testing
Various tests have been performed to ensure scientific accuracy and the optimal computing performance. Everything involving unit testing can be found in this [repository](https://github.com/mikeph52/biogenie-test).

## Installation
### GNU Linux
- Download BioGenie from Releases, or with wget:
```
wget https://github.com/mikeph52/BioGenie/releases/download/v.1.0.0/biogenie_linux_1.0.0
```
- Run "chmod +x" first.
```
sudo chmod +x biogenie_linux_1.0.0
```
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_linux_1.0.0 /usr/local/bin/biogenie
```
- If you need to build from source:
```
git clone https://github.com/mikeph52/BioGenie.git
g++ main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### macOS
It supports Intel and Apple silicon processors.
- Download BioGenie from Releases, or with curl:
```
curl -l https://github.com/mikeph52/BioGenie/releases/download/v.1.0.0/biogenie_macos_1.0.0
```
- Run "chmod +x" first.
```
sudo chmod +x biogenie_macos_1.0.0
```
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_macos_1.0.0 /usr/local/bin/biogenie
```
- If you need to build from source(probably not):
```
git clone https://github.com/mikeph52/BioGenie.git
g++ -std=c++17 main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### MS Windows
> [!IMPORTANT]
> _Windows is back(sort of). From version v.0.28.0, windows is supported. Still, this version is not going to be prioritized._

It requires Windows 10 or later.
- Download BioGenie from Releases.
- Add the executable to PATH(https://stackoverflow.com/questions/44272416/add-a-folder-to-the-path-environment-variable-in-windows-10-with-screenshots)
- Run from powershell.

~There's also a GUI Version Available on Alpha testing~.
> [!CAUTION]
> This is an Alpha testing version. It is not functional. Not for scientific use. I don't think it is possible to make and maintain a gui version for windows. The code is complicated already. It's time to move on.

## Changelog:
### 1.1.0 (March Xth, 2026)
- Genome Coverage(x) added.
- Genome Coverage(x) added in "-asmbl" pipeline.
- FASTQ support added.
- Major file-handling improvements made.

### 1.0.0 (March 3rd, 2026)
- Issue (https://github.com/mikeph52/BioGenie/issues/42) fixed.
- Documentation enhanced. 
- Minor bugs and formatting issues fixed.

_For more info checkout the [changelog](Documentation/changelog.md)._

