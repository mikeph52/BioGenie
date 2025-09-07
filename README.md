# BioGenie
BioGenie is a complete bioinformatics command line tool for macOS, GNU Linux and MS Windows, written in C++.

<img width="874" height="724" alt="Image" src="https://github.com/user-attachments/assets/7f65593b-35c1-4092-bd2c-bf581854e6fa" />

> [!IMPORTANT]
> _Windows CLI version in Beta. It's not as quick as the other versions. GUI version is available in alpha, but it's still not usable. Available in Releases(https://github.com/mikeph52/BioGenie/releases/tag/testing)_

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
- Get the complement DNA sequence ---> "-c".
- Get the reverse complement DNA sequence ---> "-rc".
- Get the codon number ---> "-nc".
- Get the mRNA ---> "-t".
- GC percentage calculation ---> "-gc".
- Generate the aminoacids(Protein chain) ---> "-p".
- Separate different sequencies in a FASTA file ---> "-ss".
- Print the different sequence headers from a FASTA file ---> "-sh".
- Trim DNA sequence ---> "-tr".
- Get the purine/pyrimidine ratio --> "-pp".
- Calculate melting temperature (Tm) of DNA sequences using the Wallace Rule(only valid for oligos <20bp) ---> "-mt1".
- Calculate melting temperature (Tm) of DNA sequences using the SantaLucia 1998 nearest-neighbor method ---> "-mt2".
- Custom preset pipeline 1 ---> "-pip1".
- Custom preset pipeline 2 ---> "-pip2".
- Coloured cDNA sequence ---> "-cc".
- Get the Open Reading Frame(ORF) ---> "-orf".

More functions will be added in the future.

## Installation
### GNU Linux
- Download BioGenie from Releases, or with wget:
```
wget https://github.com/mikeph52/BioGenie/releases/download/v.0.16.0/biogenie_linux_0.16.0
``` 
- Run "chmod +x" first(Replace 0.x.x with the correct version).
```
sudo chmod +x biogenie_linux_0.x.x
``` 
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_linux_0.x.x /usr/local/bin/biogenie
```
### macOS
- Download BioGenie from Releases, or with curl:
```
curl -l https://github.com/mikeph52/BioGenie/releases/download/v.0.16.0/biogenie_macos_0.16.0
``` 
- Run "chmod +x" first(Replace 0.x.x with the correct version).
```
sudo chmod +x biogenie_macos_0.x.x
``` 
- Move it to bin folder by executing the following command:
```
sudo mv biogenie_macos_0.x.x /usr/local/bin/biogenie
```
### MS Windows
It requires Windows 10 or later.
- Download BioGenie from Releases.
- Add the executable to PATH(https://stackoverflow.com/questions/44272416/add-a-folder-to-the-path-environment-variable-in-windows-10-with-screenshots)
- Run from powershell.
There's also a GUI Version Available on Alpha testing.
> [!CAUTION]
> This is an Alpha testing version. It is not functional. Not for scientific use.

> [!IMPORTANT]
> _Windows version in Beta. Has a lot of room for improvements. In MS Windows Home Edition, there's a reported bug, causing an exception handle(https://github.com/mikeph52/BioGenie/issues/8)._
## Changelog:
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
