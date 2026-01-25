# BioGenie
BioGenie is a complete bioinformatics command line tool for macOS, GNU Linux and MS Windows, written in C++.

<img width="1002" height="770" alt="Image" src="https://github.com/user-attachments/assets/4dfb77db-721d-42e0-ade9-73115f12057e" />


It currently supports fasta formats(.fasta, .fa).
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
Read [Documentation](documentation.md) for more information (*References included*).
- Get the complement DNA sequence ---> "-c".
- Get the reverse complement DNA sequence ---> "-rc".
- Get the codon number ---> "-nc".
- Get the mRNA ---> "-t".
- Motif search ---> "-mf".
- Ambiguous bases statistics ---> "-amb".
- GC percentage calculation ---> "-gc".
- Generate the aminoacids(Protein chain) ---> "-p".
- Generate the Protein chain with color ---> "-pc".
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
- Custom preset pipeline 1 ---> "-pip1".
- Custom preset pipeline 2 ---> "-pip2".
- Custom preset pipeline 3 ---> "-pip3".


More functions will be added in the future.

> [!NOTE]
> _If you have any suggestions for new features or a bug encountered, create an Issue or send me a message at: mikeph526@outlook.com. I'm happy to help._

## Installation
### GNU Linux
- Download BioGenie from Releases, or with wget:
```
wget https://github.com/mikeph52/BioGenie/releases/download/v.0.27.0/biogenie_linux_0.27.0
``` 
- Run "chmod +x" first(Replace 0.x.x with the correct version).
```
sudo chmod +x biogenie_linux_0.x.x
``` 
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_linux_0.x.x /usr/local/bin/biogenie
```
- If you need to build from source:
```
git clone https://github.com/mikeph52/BioGenie.git
g++ main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### macOS
- Download BioGenie from Releases, or with curl:
```
curl -l https://github.com/mikeph52/BioGenie/releases/download/v.0.27.0/biogenie_macos_0.27.0
``` 
- Run "chmod +x" first(Replace 0.x.x with the correct version).
```
sudo chmod +x biogenie_macos_0.x.x
``` 
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_macos_0.x.x /usr/local/bin/biogenie
```
- If you need to build from source(probably not):
```
git clone https://github.com/mikeph52/BioGenie.git
g++ -std=c++17 main.cpp -o biogenie
sudo mv biogenie /usr/local/bin/
```
### MS Windows (Deprecated)
> [!IMPORTANT]
> _The Windows version is not being maintained and probably never will. Last available version: 0.14.0(https://github.com/mikeph52/BioGenie/releases/tag/v.0.14.0). Use WSL(Windows Subsystem for Linux) instead._
It requires Windows 10 or later.
- Download BioGenie from Releases.
- Add the executable to PATH(https://stackoverflow.com/questions/44272416/add-a-folder-to-the-path-environment-variable-in-windows-10-with-screenshots)
- Run from powershell.
There's also a GUI Version Available on Alpha testing.
> [!CAUTION]
> This is an Alpha testing version. It is not functional. Not for scientific use. I don't think it is possible to make and maintain a gui version for windows. The code is complicated already. It's time to move on.

## Changelog:
- 0.27.0:
Multi core support added for pip1, pip2 and pip3(https://github.com/mikeph52/BioGenie/issues/30), (https://github.com/mikeph52/BioGenie/issues/31). Ambiguous bases statistics added. Execution timer added. Minor bugs fixed.

- 0.26.0:
Hydrogen bonds calculator and protein chain with color added.
- 0.25.0:
Protein mollecular weight, isoelectric point, extinction coeficient and pipeline 3 added. Issues https://github.com/mikeph52/BioGenie/issues/25, https://github.com/mikeph52/BioGenie/issues/27 fixed. Major performance improvements made. Minor bugs fixed. 
- 0.24.0:
Motif search added. Issues https://github.com/mikeph52/BioGenie/issues/19, https://github.com/mikeph52/BioGenie/issues/21, https://github.com/mikeph52/BioGenie/issues/23, https://github.com/mikeph52/BioGenie/issues/24 fixed. Minor formatting issues fixed.
- 0.23.0:
DNA sequence with color added.
- 0.22.0:
Base Pair calculation added to Pipeline 2. Minor bugs fixed.
- 0.21.0:
Base Pair calculation added. Minor bugs fixed.
- 0.20.0:
Codon Usage Bias(CUB) calculation and export to csv function added.
- 0.19.0:
(https://github.com/mikeph52/BioGenie/issues/15)
Minor bugs fixed.
- 0.18.0:
Generate mRNA FASTA added. Minor format bugs fixed.
- 0.17.3:
(https://github.com/mikeph52/BioGenie/issues/13)
Issue Fixed.
- 0.17.2:
(https://github.com/mikeph52/BioGenie/issues/12)
Issue Fixed.
- 0.17.0:
Generate cDNA FASTA and reverse compliment FASTA function added.
- 0.16.0:
Open Reading Frames finder function added. Minor bugs fixed.
- 0.15.0:
(https://github.com/mikeph52/BioGenie/issues/9)
Fixed Melting Temperature Calculator functiion. SantaLucia 1998 nearest-neighbor method added as "-mt2". Pipeline 2 function also fixed. Now uses SantaLucia 1998 nearest-neighbor method for more accurate calculations. Minor bugs fixed.
- 0.14.0:
Coloured cDNA sequence. Minor adjustments made and bugs fixed.
- 0.13.0:
Purine/pyrimidine ratio and Melting temperature(Tm) calculator functions added. Preset Pipeline 2 added. Minor format fixes.
- 0.12.2:
Minor format fixes.
- 0.12.0:
(https://github.com/mikeph52/BioGenie/issues/7)
Protein function output fixed, "Pipeline 1" added and minor bugs fixed.
- 0.11.0:
FASTA file verification and integrity checker added.
- 0.10.0:
(https://github.com/mikeph52/BioGenie/issues/6)
FASTA sequence header print function added. DNA trimmer function added. Windows support added. Minor bugs fixed and quality of life improvements.
- 0.9.0:
Linux support added. FASTA sequencies separator added and minor bugs fixed.
- 0.8.0:
Protein chain option added and minor bugs fixed.
- 0.7.0:
(https://github.com/mikeph52/BioGenie/issues/4)
Documentation fixed.
- 0.6.0:
(https://github.com/mikeph52/BioGenie/issues/3)
Documentation added
- 0.5.0:
(https://github.com/mikeph52/BioGenie/issues/2)
Function operator added.
- 0.4.0:
(https://github.com/mikeph52/BioGenie/issues/1)   Reverse complement DNA function added.
- 0.3.0:
Number of codone calculator added.
- 0.2.3:
GC calculator, complimentary DNA and transcripted RNA functions added.
- 0.1.0:
First Version.
