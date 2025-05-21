# BioGenie
BioGenie is a complete bioinformatics command line tool, written in C++.

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

### Documentation
BioGenie uses functions to execute different tools for different applications.
- Get the complement DNA sequence --> "-c".
- Get the reverse complement DNA sequence --> "-rc".
- Get the codon number --> "-nc".
- Get the mRNA --> "-t".
- GC percentage calculation --> "-gc".

More functions will be added in the future.

### Installation
- Download BioGenie from Releases.
- Run "chmod +x" first.
```
sudo chmod +x biogenie
``` 
- Move it to bin folder by executing the following command:
```
sudo mv biogenie /usr/local/bin/biogenie
```

### Changelog:
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
