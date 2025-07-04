# BioGenie
BioGenie is a complete bioinformatics command line tool for macOS, GNU Linux and MS Windows (Windows 10 or later), written in C++.
> [!IMPORTANT]
> Windows version in Beta. It's not as quick as the other versions.
<img width="762" alt="Image" src="https://github.com/user-attachments/assets/fb9d961e-e7e4-4a99-a401-95281a9f96e0" />

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
- Generate the aminoacids(Protein chain) ---> '-p'.
- Separate different sequencies in a FASTA file ---> '-ss'.
- Print the different sequence headers from a FASTA file ---> '-sh'.
- Trim DNA sequence ---> 'tr'.

More functions will be added in the future.

## Installation
### GNU Linux
- Download BioGenie from Releases, or with wget:
```
wget https://github.com/mikeph52/BioGenie/releases/download/v.0.10.0/biogenie_linux_0.10.0
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
curl -l https://github.com/mikeph52/BioGenie/releases/download/v.0.10.0/biogenie_macos_0.10.0
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
- Download BioGenie from Releases.
- Add the executable to PATH(https://stackoverflow.com/questions/44272416/add-a-folder-to-the-path-environment-variable-in-windows-10-with-screenshots)
- Run from powershell.
> [!IMPORTANT]
> Windows version in Beta. Has a lot of room for improvements.
## Changelog:
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
