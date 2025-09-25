# Documentation
BioGenie uses functions to execute different tools for different applications. References are provided.
- ###  Get the complement DNA sequence.
Biogenie returns the complement DNA sequence (cDNA) from the whole FASTA sequence according to the rules of Complementarity (Watson et al. 2014) , by typing the "-c" flag. N nucleotides are trimmed in order to save space. Complementary DNA (cDNA) refers to a laboratory-synthesized DNA molecule that is a complementary copy of the messenger RNA (mRNA). It is distinct from genomic DNA as it lacks promoters and introns, representing the actively expressed genes at the time of harvesting (Pelley 2007). Also, with "-cw" flag the cDNA sequence FASTA can be generated and with the "-cc" flag, a coloured cDNA sequence can be printed using the common ab1 file colours for each nucleotide (NWABR 2012).

- ### Get the reverse complement DNA sequence.
Biogenie returns the reverse-complement DNA sequence from the whole FASTA sequence according to the rules of Complementarity (Watson et al. 2014), by typing the "-rc" flag.The reverse-complement of a sequence is usefull if the reverse strand contains an Open Reading Frame (ORF)(Jannik N et. al 2005).

- ### Get the mRNA sequence.
By typing the "-t" flag, the program returns the mRNA sequence from a FASTA File according to the rules of Complementarity (Watson et al. 2014). 

- ### GC Content calculation.
Genomic GC (Guanine-Cytosine) content is a fundamental molecular trait linked with many key genomic features such as codon and amino acid use (Mahajan 2022). Biogenie calculates the GC Content percentage by typing the "-gc" flag, according to the formula bellow (Madigan 2003):

![GC% formula](https://latex.codecogs.com/svg.latex?\bg_black\color{White}GC\%=\frac{Number\;of\;G+Number\;of\;C}{Total\;number\;of\;bases}\times100)

- ### Get the codon number.
Number of codons are calculated using the "-nc" flag (Yanofsky 2007).

- ### Generate the aminoacids(Protein chain).
Biogenie can generate a protein chain from a FASTA file, using the standard RNA codon table(Muto et al. 1987).

- ### Purine/Pyrimidine ratio.
Get the Purine/Pyrimidine ratio using the "-pp" flag according to the Chargaff's rules (Chargaff et al. 1952).

- ### Melting Temperature (Tm) of DNA sequences (Wallace Rule).
Calculate melting temperature (Tm) of DNA sequences using the Wallace Rule for short DNA oligonucleotides (approximately 14-20bp) with the "-mt1" flag. The formula is shown bellow(Thein 1986): 

![Tm formula](https://latex.codecogs.com/svg.latex?\bg_black\color{White}Tm\=\{2(A+T)}{4(G+C)})

*Note!: While it provides a quick estimate and is useful for primers 14–20 nucleotides in length, it is less accurate than more complex thermodynamic or nearest-neighbor methods and doesn't account for factors like salt concentration or Base stacking interactions.*

- ### Melting Temperature (Tm) of DNA sequences (Nearest-Neighbor).
Calculate melting temperature (Tm) of DNA sequences using the SantaLucia Nearest-Neighbor Thermodynamic model for for DNA oligonucleotides > 20bp (SantaLucia 1998):

![Tm formula](https://latex.codecogs.com/svg.latex?\bg_black\color{White}Tm=\frac{\Delta&space;H^\circ}{\Delta&space;S^\circ&plus;R\cdot\ln(C)}-273.15&plus;16.6\cdot\log_{10}[Na^&plus;])


## References
- Watson, James, Cold Spring Harbor Laboratory, Tania A. Baker, Massachusetts Institute of Technology, Stephen P. Bell, Massachusetts Institute of Technology, Alexander Gann, Cold Spring Harbor Laboratory, Michael Levine, University of California, Berkeley, Richard Losik, Harvard University ; with Stephen C. Harrison, Harvard Medical (2014). *Molecular biology of the gene (Seventh ed.)* . Boston: Benjamin-Cummings Publishing Company. ISBN 978-0-32176243-6.

- John W. Pelley,
*18 - Recombinant DNA and Biotechnology*,
Editor(s): John W. Pelley,
Elsevier's Integrated Biochemistry,
Mosby,
2007,
Pages 159-167,
ISBN 9780323034104,
https://doi.org/10.1016/B978-0-323-03410-4.50024-9.
(https://www.sciencedirect.com/science/article/pii/B9780323034104500249)

- Jannik N. Andersen, Robert L. Del Vecchio, Natarajan Kannan, James Gergel, Andrew F. Neuwald, Nicholas K. Tonks,
*Computational analysis of protein tyrosine phosphatases: practical guide to bioinformatics and data resources,
Methods,
Volume 35, Issue 1*,
2005,
Pages 90-114,
ISSN 1046-2023,
https://doi.org/10.1016/j.ymeth.2004.07.012.
(https://www.sciencedirect.com/science/article/pii/S1046202304001756)

- Mahajan, S., & Agashe, D. (2022). *Evolutionary jumps in bacterial GC content.* G3 (Bethesda, Md.), 12(8), jkac108. https://doi.org/10.1093/g3journal/jkac108

- Madigan, MT. and Martinko JM. (2003). *Brock biology of microorganisms (10th ed.)*. Pearson-Prentice Hall. ISBN 978-84-205-3679-8.

- Yanofsky C (9 March 2007). *Establishing the Triplet Nature of the Genetic Code*. Cell. 128 (5): 815–818. doi:10.1016/j.cell.2007.02.029. PMID 17350564. S2CID 14249277.

- Muto, A.; Osawa, S. (January 1987). *The guanine and cytosine content of genomic DNA and bacterial evolution*. Proc Natl Acad Sci USA. 84 (1): 166–9. Bibcode:1987PNAS...84..166M. doi:10.1073/pnas.84.1.166. PMC 304163. PMID 3467347.

- Northwest Association for Biomedical Research. (2012, August 14). ©Northwest Association for Biomedical Research—updated August 14, 2012 1. *Analyzing A DNA Sequence Chromatogram*. https://www.nwabr.org/sites/default/files/Analyzing-A-DNA-Sequence-Chromatogram.pdf 

- Elson, D., Chargaff, E. *On the desoxyribonucleic acid content of sea urchin gametes*. Experientia 8, 143–145 (1952). https://doi.org/10.1007/BF02170221

- Thein S. L., Lynch J. R., Weatherall D. J., et al. *Direct detection of haemoglobin E with synthetic oligonucleotides.* The Lancet, 1986, 327(8472): 93.

- SantaLucia J., Jr (1998). *A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.* Proceedings of the National Academy of Sciences of the United States of America, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460
