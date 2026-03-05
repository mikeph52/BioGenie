// BioGenie by mikeph_ 2025-2026
// Current version 1.0.0 3/3/2026
#include <iostream>
#include <string>
#include <cctype>
#include <iomanip> 
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <functional> 
#include <thread> //new stuff-multi thread
#include <mutex> //for mutual execution
#include <queue> //for thread pool
#include <condition_variable> //synchronization primitive
#include <map>
// Public Functions 
void title(){
    std::cout << "-----------------------\n";
    std::cout << "BioGenie 1.0.0 \nby mikeph_ 2025-2026\n\n";    
}
void helpme() {
    std::cout << "\n\nNAME\n";
    std::cout << "       biogenie - A complete bioinformatics command line tool, written in C++. \n\n";
    std::cout << "USAGE\n";
    std::cout << "       biogenie [function] [FASTA_file_path]\n\n";
    std::cout << "DESCRIPTION\n";
    std::cout << "       BioGenie uses functions to execute different tools for different applications.\n";
    std::cout << "       Read Documentation for more information (References included).\n\n";
    std::cout << "PIPELINES\n";
    std::cout << "       -nucleo  Nucleostats: DNA statistics (primer design)\n";
    std::cout << "       -prot    Proteostats: Protein statistics and structural properties\n";
    std::cout << "       -asmbl   Assemblystats: Genome assembly statistics\n\n";
    std::cout << "DNA SEQUENCE UTILITIES\n";
    std::cout << "       -c       Get the complement DNA sequence\n";
    std::cout << "       -rc      Get the reverse complement DNA sequence\n";
    std::cout << "       -nc      Get the codon number\n";
    std::cout << "       -t       Get the RNA sequence (transcription)\n";
    std::cout << "       -gc      Calculate GC percentage (and AT percentage)\n";
    std::cout << "       -pp      Get the purine/pyrimidine ratio\n";
    std::cout << "       -amb     Ambiguous base statistics\n";
    std::cout << "       -bp      Calculate the number of base pairs (bp)\n";
    std::cout << "       tr       Trim DNA (0-based indexing)\n";
    std::cout << "       -sc      Get DNA sequence with colour\n";
    std::cout << "       -cc      Get complement DNA sequence with colour\n\n";
    std::cout << "PROTEIN UTILITIES\n";
    std::cout << "       -p       Generate amino acids (protein chain)\n";
    std::cout << "       -pc      Generate amino acids with colour\n";
    std::cout << "       -pca     Color protein chain from FASTA file\n";
    std::cout << "       -pi      Calculate isoelectric point (pI)\n";
    std::cout << "       -mw      Calculate molecular weight (kDa)\n";
    std::cout << "       -ec      Calculate extinction coefficient\n";
    std::cout << "       -prot    Proteostats: protein statistics and structural properties\n\n";
    std::cout << "FASTA UTILITIES\n";
    std::cout << "       -ss      Separate sequences in a FASTA file\n";
    std::cout << "       -sh      Print FASTA sequence headers\n";
    std::cout << "       -cw      Generate cDNA FASTA\n";
    std::cout << "       -rcw     Generate reverse cDNA FASTA\n";
    std::cout << "       -tw      Generate mRNA FASTA\n\n";
    std::cout << "CODON & MOTIF ANALYSIS\n";
    std::cout << "       -cub     Calculate Codon Usage Bias (CUB)\n";
    std::cout << "       -wcub    Export Codon Usage Bias (CUB) to CSV\n";
    std::cout << "       -mf      Search motifs\n\n";
    std::cout << "THERMODYNAMICS\n";
    std::cout << "       -mt1     Melting temperature (Wallace Rule, <20bp)\n";
    std::cout << "       -mt2     Melting temperature (SantaLucia 1998 method)\n";
    std::cout << "       -hb      Calculate hydrogen bonds of dsDNA\n\n";
    std::cout << "SEQUENCE ALIGNMENT\n";
    std::cout << "       -nw      Global alignment (Needleman-Wunsch)\n";
    std::cout << "       -sw      Local alignment (Smith-Waterman)\n\n";
    std::cout << "GENOME ANALYSIS\n";
    std::cout << "       -nucleo  Nucleostats: DNA statistics (primer design)\n";
    std::cout << "       -asmbl   Assemblystats: genome assembly statistics\n";
    std::cout << "       -orf     Identify Open Reading Frames (ORFs)\n";
    std::cout << "       -cx      Genome Coverage(x)(Depth)\n\n";
    std::cout << "AUTHOR\n";
    std::cout << "       BioGenie, developed and maintained by Mike Philippakis, Github:mikeph52,\n";
    std::cout << "       under GNU GENERAL PUBLIC LICENSE Version 3, 2025-2026.\n\n";
    std::cout << "CONTACT\n";
    std::cout << "       If you have any suggestions for new features or a bug encountered,\n";
    std::cout << "       create an Issue or send me a message at:mikeph526@outlook.com.\n";
    std::cout << "                             I'm happy to help.\n";
    std::cout << "       Github page:https://github.com/mikeph52/BioGenie\n\n";
    std::cout << "---------------------------------------------------------------------\n";
}
void message(){
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[Pipelines]:\n";
        std::cerr << "[-nucleo Nucleostats][-prot Proteostats][-asmbl Assemblystats][-nw Needleman-Wunsch][-sw Smith-Waterman]\n\n";
        std::cerr << "[Single commands]:\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence][-t RNA][-bp Base Pairs]\n";
        std::cerr << "[-gc GC percentage calculator][-p protein chain][-pp purine/pyrimidine ratio][-nc codon number]\n";
        std::cerr << "[-mt1 melting temp.(Wallace rule)][-mt2 melting temp.(Nearest-neighbour)][-orf ORF Finder]\n";
        std::cerr << "[-ec Extinction Coefficient][-hb Hydrogen Bonds][-amb Ambiguous stats][-mw prot kDa]\n";
        std::cerr << "[-cub Codon Usage Bias][-mf Find MOTIFs][-pi Isoelectric Point][-cx Coverage(x)]\n\n";
        std::cerr << "[Utilities]:\n";
        std::cerr << "[-ss FASTA sequencies separator][-sh FASTA sequencies headers][-tr DNA Trimmer][-sc colour sequence]\n";
        std::cerr << "[-wcub Codon Usage Bias to CSV][-cw generate cDNA fasta][-rcw Reverse cDNA fasta][-tw mRNA fasta]\n";
        std::cerr << "[-pc protein chain w color][-pca protein fasta w color][-cc cDNA coloured]\n";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
}
//Genetic code
const std::unordered_map<std::string, char>codonTable = {
    {"ATA",'I'}, {"ATC",'I'}, {"ATT",'I'}, {"ATG",'M'},
    {"ACA",'T'}, {"ACC",'T'}, {"ACG",'T'}, {"ACT",'T'},
    {"AAC",'N'}, {"AAT",'N'}, {"AAA",'K'}, {"AAG",'K'},
    {"AGC",'S'}, {"AGT",'S'}, {"AGA",'R'}, {"AGG",'R'},
    {"CTA",'L'}, {"CTC",'L'}, {"CTG",'L'}, {"CTT",'L'},
    {"CCA",'P'}, {"CCC",'P'}, {"CCG",'P'}, {"CCT",'P'},
    {"CAC",'H'}, {"CAT",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
    {"CGA",'R'}, {"CGC",'R'}, {"CGG",'R'}, {"CGT",'R'},
    {"GTA",'V'}, {"GTC",'V'}, {"GTG",'V'}, {"GTT",'V'},
    {"GCA",'A'}, {"GCC",'A'}, {"GCG",'A'}, {"GCT",'A'},
    {"GAC",'D'}, {"GAT",'D'}, {"GAA",'E'}, {"GAG",'E'},
    {"GGA",'G'}, {"GGC",'G'}, {"GGG",'G'}, {"GGT",'G'},
    {"TCA",'S'}, {"TCC",'S'}, {"TCG",'S'}, {"TCT",'S'},
    {"TTC",'F'}, {"TTT",'F'}, {"TTA",'L'}, {"TTG",'L'},
    {"TAC",'Y'}, {"TAT",'Y'}, {"TAA",'*'}, {"TAG",'*'},
    {"TGC",'C'}, {"TGT",'C'}, {"TGA",'*'}, {"TGG",'W'}
};
// ANSI ESCAPE CODES FOR COLOR
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define BLUE    "\033[34m"
#define YELLOW  "\033[33m"
#define CYAN    "\033[36m"
#define GREEN   "\033[32m"
#define WHITE   "\033[37m"
#define BRED    "\033[41m"
//File Verifiers
class FastaVerifier {
    public:
        explicit FastaVerifier(const std::string& filePath) : filename(filePath) {}
        bool verify() {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Cannot open file: " << filename << std::endl;
                return false;
            }
            std::string line;
            bool inSequence = false;
            int lineNumber = 0;

            while (std::getline(fastaFile, line)) {
                ++lineNumber;
                if (line.empty()) continue;
                if (line[0] == '>') {
                    inSequence = true;
                    if (line.length() == 1) {
                        reportError(lineNumber, "Empty header (no sequence identifier).");
                        return false;
                    }
                } else {
                    if (!inSequence) {
                        reportError(lineNumber, "Sequence data appears before a header.");
                        return false;
                    }
                    if (!isValidSequenceLine(line)) {
                        reportError(lineNumber, "Invalid characters found in sequence line.");
                        return false;
                    }
                }
            }
            return true;
        }
    private:
        std::string filename;
        bool isValidSequenceLine(const std::string& line) {
            for (char ch : line) {
                if (std::isspace(ch)) continue; // fix for validating whitespaces
                char upper = std::toupper(ch);
                if (!(upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' || upper == 'N')) {
                    return false;
                }
            }
            return true;
        }
        void reportError(int lineNumber, const std::string& message) {
            std::cerr << "Error at line " << lineNumber << ": " << message << std::endl;
        }
};
class FastqVerifier {
    public:
        struct Result {
            bool valid = true;
            size_t records = 0;
            size_t line_number = 0;
            std::string error_message;
        };
        Result verify(const std::string& filename) {
            std::ifstream file(filename);
            Result result;
            if (!file) {
                result.valid = false;
                result.error_message = "Cannot open file.";
                return result;
            }
            std::string line1, line2, line3, line4;
            while (true) {
                if (!std::getline(file, line1)) break;
                if (!std::getline(file, line2) ||
                    !std::getline(file, line3) ||
                    !std::getline(file, line4)) {
                    result.valid = false;
                    result.error_message = "Incomplete FASTQ record.";
                    return result;
                }
                result.line_number += 4;
                if (line1.empty() || line1[0] != '@') {
                    result.valid = false;
                    result.error_message = "Header line does not start with '@'";
                    return result;
                }
                if (line3.empty() || line3[0] != '+') {
                    result.valid = false;
                    result.error_message = "Separator line does not start with '+'";
                    return result;
                }
                if (line2.length() != line4.length()) {
                    result.valid = false;
                    result.error_message = "Sequence and quality lengths differ.";
                    return result;
                }
                if (!valid_sequence(line2)) {
                    result.valid = false;
                    result.error_message = "Invalid character in sequence.";
                    return result;
                }
                result.records++;
            }
            return result;
        }
    private:
        bool valid_sequence(const std::string& seq) {
            for (char c : seq) {
                switch (c) {
                    case 'A': case 'C': case 'G': case 'T': case 'N':
                    case 'a': case 'c': case 'g': case 't': case 'n':
                        break;
                    default:
                        return false;
                }
            }
            return true;
        }
};
void FileVerifier(const std::string& filename) {
    std::string line;
    int lineNumber = 0;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file\n";
        return;
    }
    while (std::getline(file, line)) {
        ++lineNumber;
        if (line.empty()) continue;
        if (line[0] == '>') {
            // FASTA verifier
            FastaVerifier fastaverifier(filename);
            if (fastaverifier.verify()) {
                std::cout << "FASTA file status  [OK]\n";
            } else {
                std::cerr << "FASTA file status  [FAULT]\n";
            }
            return;
        } 
        else if (line[0] == '@') {
            // FASTQ verifier
            FastqVerifier verifier;
            auto result = verifier.verify(filename);
            if (result.valid) {
                std::cout << "FASTQ file status  [OK]\n";
                std::cout << "Records: " << result.records << "\n";
            } else {
                std::cout << "FASTQ file status  [FAULT]\n"
                          << result.error_message
                          << " at line " << result.line_number << "\n";
            }
            return;
        }
    }
    std::cerr << "Unknown file format\n";
}
// Functions
class GCCalc {
  private:
        double GCContent(const std::string& sequence) {
            int gcCount = 0;
            int validBases = 0;

            for (char base : sequence) {
                char upperBase = std::toupper(base);
                if (upperBase == 'G' || upperBase == 'C') {
                    gcCount++;
                    validBases++;
                } else if (upperBase == 'A' || upperBase == 'T') {
                    validBases++;
                }
            }

            if (validBases == 0) return 0.0;

            return (static_cast<double>(gcCount) / validBases) * 100.0;
        }
  public:
      void FASTA_loader(const std::string& filename) {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file: " << filename << "\n";
                exit(1);
            }
            std::string line, Header, sequence;
            std::cout << "\n-----------------------------------\n";
            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        double gcContent = GCContent(sequence);
                        std::cout << Header << ":\nGC Content = " << std::fixed << std::setprecision(2) << gcContent << "%\n";
                        std::cout << "AT Content = " << std::fixed << std::setprecision(2) << 100 - gcContent << "%\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    Header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
                double gcContent = GCContent(sequence);
                std::cout << Header << ":\nGC Content = " << std::fixed << std::setprecision(2) << gcContent << "%\n";
                std::cout << "AT Content = " << std::fixed << std::setprecision(2) << 100 - gcContent<< "%\n";
            }
            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class DNAcomplementary{
    private:
        char Complement(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return 'N'; // Unknown base
            }
        }
        //DNA complement strand init
        std::string ComplementStrand(const std::string& sequence) const {
            std::string complement;
            complement.reserve(sequence.size());
            for (char base : sequence) {
                complement += Complement(base);
            }
            return complement;
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string complement = ComplementStrand(sequence);
                        std::cout << ">" << header << " (complement)\n" << complement << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
    
            if (!sequence.empty()) {
                std::string complement = ComplementStrand(sequence);
                std::cout << ">" << header << " (complement)\n" << complement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class ReverseComplementDNA{
    private:
        char Complement(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return 'N'; // Unknown base
            }
        }
       
        std::string ReverseComplementStrand(const std::string& sequence) const {
            std::string revComplement;
            revComplement.reserve(sequence.size());
            for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
                revComplement += Complement(*it);
            }
            return revComplement;
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string revComplement = ReverseComplementStrand(sequence);
                        std::cout << ">" << header << " (reverse complement)\n" << revComplement << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
    
            if (!sequence.empty()) {
                std::string revComplement = ReverseComplementStrand(sequence);
                std::cout << ">" << header << " (reverse complement)\n" << revComplement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }

};
class Transcription{
    private:
        char transRNA(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'A';
                case 'T': return 'U';
                case 'G': return 'G';
                case 'C': return 'C';
                default:  return 'N'; // Unknown base
            }
        }
    
        std::string RNAStrand(const std::string& sequence) const {
            std::string RNAseq;
            RNAseq.reserve(sequence.size());
            for (char base : sequence) {
                RNAseq += transRNA(base);
            }
            return RNAseq;
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string complement = RNAStrand(sequence);
                        std::cout << ">" << header << "RNA:\n" << complement << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
    
            if (!sequence.empty()) {
                std::string complement = RNAStrand(sequence);
                std::cout << ">" << header << "RNA:\n" << complement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class CodonNumber{
    private:
    int CodonCount(const std::string& sequence) const {
        int validBases = 0;

        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            if (upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G') {
                validBases++;
            }
        }
        return validBases / 3;
    }

    public:
    void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        int codons = CodonCount(sequence);
                        std::cout << ">" << header << "Codon count:\n" << codons << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
    
            if (!sequence.empty()) {
                int codons = CodonCount(sequence);
                std::cout << ">" << header << "Codon count:\n" << codons << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";

            fastaFile.close();
        }

};
class ProteinChain{
    private:
        std::string translateToAminoAcids(const std::string& sequence) {
            std::string protein;

            for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
                std::string codon = sequence.substr(i, 3);
                for (char& c : codon) c = std::toupper(c);

                if (codonTable.count(codon)) {
                    protein += codonTable.at(codon);
                } else {
                    protein += 'X';  // Unknown codon
                }
            }

        return protein;
    }

    public:
        void FASTA_loader(const std::string& filename)  {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }
            std::string line;
            std::string header;
            std::string sequence;
            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string protein = translateToAminoAcids(sequence);
                        std::cout << ">" << header << "Protein:\n" << protein<< "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
    
            if (!sequence.empty()) {
                std::string complement = translateToAminoAcids(sequence);
                std::cout << ">" << header << "Protein:\n" << complement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class ProteinColor {
    private:
        std::string translateToAminoAcids(const std::string& sequence) {
            std::string protein;

            for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
                std::string codon = sequence.substr(i, 3);
                for (char& c : codon) c = std::toupper(c);

                if (codonTable.count(codon))
                    protein += codonTable.at(codon);
                else
                    protein += 'X';   // Unknown codon
            }
            return protein;
        }
        std::string colorAA(char aa) {
        switch (aa) {
            case 'A': case 'V': case 'L': case 'I': case 'M':
            case 'F': case 'W': case 'Y':
                return std::string(YELLOW) + aa + RESET;
            case 'S': case 'T': case 'N': case 'Q':
                return std::string(CYAN) + aa + RESET;
            case 'K': case 'R': case 'H':
                return std::string(BLUE) + aa + RESET;
            case 'D': case 'E':
                return std::string(RED) + aa + RESET;
            case 'G': case 'P': case 'C':
                return std::string(GREEN) + aa + RESET;
            case 'X':
                return std::string(WHITE BRED) + aa + RESET;
            default:
                return std::string(WHITE) + aa + RESET;
            }
        }
    public:
        void FASTA_loader(const std::string& filename) {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line, header, sequence;
            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    // Print previous sequence if exists
                    if (!sequence.empty()) {
                        std::string protein = translateToAminoAcids(sequence);

                        std::string coloredProtein;
                        for (char aa : protein)
                            coloredProtein += colorAA(aa);

                        std::cout << ">" << header << " Protein:\n"
                                << coloredProtein << "\n\n";
                        std::cout << "-----------------------------------\n";

                        sequence.clear();
                    }
                    header = line.substr(1);
                }
                else {
                    sequence += line;
                }
            }

            // Print last sequence
            if (!sequence.empty()) {
                std::string protein = translateToAminoAcids(sequence);

                std::string coloredProtein;
                for (char aa : protein)
                    coloredProtein += colorAA(aa);

                std::cout << ">" << header << " Protein:\n"
                        << coloredProtein << "\n";
            }

            std::cout << "-----------------------------------\n";
            std::cout << "Process completed.\n";

            fastaFile.close();
        }
};
class ColorForAminoAcids {
    private:    
        std::string colorAA(char aa) {
            aa = std::toupper(aa);
            switch (aa) {
                case 'A': case 'V': case 'L': case 'I': case 'M':
                case 'F': case 'W': case 'Y':
                    return std::string(YELLOW) + aa + RESET;
                case 'S': case 'T': case 'N': case 'Q':
                    return std::string(CYAN) + aa + RESET;
                case 'K': case 'R': case 'H':
                    return std::string(BLUE) + aa + RESET;
                case 'D': case 'E':
                    return std::string(RED) + aa + RESET;
                case 'G': case 'P': case 'C':
                    return std::string(GREEN) + aa + RESET;
                case 'X':
                    return std::string(WHITE BRED) + aa + RESET;
                default:
                    return std::string(WHITE) + aa + RESET;
                }
        }
    public:
        void FASTA_loader(const std::string& filename) {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }
            std::string line, header, sequence;
            
            std::cout << "\n-----------------------------------\n";
            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;
                if (line[0] == '>') {
                    // Print previous sequence if exists
                    if (!sequence.empty()) {
                        std::cout << ">" << header << "\n";
                        int count = 0;
                        for (char aa : sequence) {
                            std::cout << colorAA(aa);
                            if (++count % 60 == 0)
                                std::cout << "\n";
                        }
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                }
                else {
                    sequence += line;
                }
            }
            // Print last sequence
            if (!sequence.empty()) {
                std::cout << ">" << header << "\n";
                int count = 0;
                for (char aa : sequence) {
                    std::cout << colorAA(aa);
                    if (++count % 60 == 0)
                        std::cout << "\n";
                }
            }
            std::cout << "\n-----------------------------------\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class FASTAChromosomeSeparator {
  public:
    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line;
        std::string header;
        std::string sequence;
        int count = 1;

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
               
                if (!header.empty()) {
                    WriteSequenceToFile(header, sequence, count++);
                    sequence.clear();
                }
                header = line.substr(1);  
            } else {
                sequence += line;
            }
        }

        if (!header.empty() && !sequence.empty()) {
            WriteSequenceToFile(header, sequence, count++);
        }

        std::cout << "\n-----------------------------------\n";
        std::cout << "Sequences from FASTA file successfully separated into " << count - 1 << " files.\n";
        std::cout << "-----------------------------------\n\n\n";
        std::cout << "Process completed.\n";

        fastaFile.close();
    }

  private:
  // generate title
    std::string SanitizeHeader(const std::string& header) const {
        std::string safeHeader = header;
        std::replace_if(safeHeader.begin(), safeHeader.end(),
                        [](char c) { return !std::isalnum(c) && c != '_'; }, '_');
        return safeHeader;
    }

    void WriteSequenceToFile(const std::string& header, const std::string& sequence, int index) const {
        std::string safeName = SanitizeHeader(header);
        std::string outputFilename = "chromosome_" + std::to_string(index) + "_" + safeName + ".fasta";

        std::ofstream outFile(outputFilename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to create output file: " << outputFilename << "\n";
            return;
        }

        outFile << ">" << header << "\n";

        const int lineLength = 80;
        for (size_t i = 0; i < sequence.size(); i += lineLength) {
            outFile << sequence.substr(i, lineLength) << "\n";
        }

        outFile.close();
        std::cout << "Saved: " << outputFilename << "\n";
    }
};
class FASTASequenceHeader {
public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }
        std::string line;
        int total_heads = 0;
        std::cout << "\n--------- Sequence Headers ---------\n";
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                total_heads++;
                std::cout << line << "\n";
            }
        }
        std::cout << "-----------------------------------\n";
        std::cout << "Total FASTA Headers: " << total_heads << "\n\n\n";
        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
class DNATrimmer {
    private:
    void printTrimmed(const std::string& header, const std::string& sequence, int start, int end) const {
        int seqLen = sequence.size();
        int trimmedStart = std::max(0, start);
        int trimmedEnd = std::min(seqLen, end);

        if (trimmedStart >= trimmedEnd) {
            std::cerr << header << "\nInvalid trim range.\n\n";
            return;
        }

        std::string trimmed = sequence.substr(trimmedStart, trimmedEnd - trimmedStart);
        std::cout << header << " (trimmed: " << trimmedStart << "-" << trimmedEnd << ")\n";
        std::cout << trimmed << "\n\n";
    }

    public:
    void FASTA_loader(const std::string& filename, int start, int end) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line;
        std::string header;
        std::string sequence;

        std::cout << "\n-------- Trimmed DNA Sequences --------\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    printTrimmed(header, sequence, start, end);
                    sequence.clear();
                }
                header = line;
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            printTrimmed(header, sequence, start, end);
        }

        std::cout << "-----------------------------------\n\n\n";
        std::cout << "Process completed.\n";
        fastaFile.close();
    }
    
};
class PurinePyrimidineRatioAnalyzer {
    private:
    void calculateRatio(const std::string& header, const std::string& sequence) const{
        int purines = 0, pyrimidines = 0;
       
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            switch (upper) {
                case 'A':
                case 'G':
                    purines++;
                    break;
                case 'C':
                case 'T':
                    pyrimidines++;
                    break;
                default:
                    break; // Skip 'N' or unknown bases
            }
        }

        std::cout << ">" << header << "\n";
        std::cout << "Purines: " << purines << "\n";
        std::cout << "Pyrimidines: " << pyrimidines << "\n";

        if (pyrimidines == 0) {
            std::cout << "Purine/Pyrimidine Ratio: Undefined (pyrimidines = 0)\n";
        } else {
            double ratio = static_cast<double>(purines) / pyrimidines;
            std::cout << "Purine/Pyrimidine Ratio: " << std::fixed << std::setprecision(3) << ratio << "\n";
        }
        std::cout << "\n-----------------------------------\n";
    }
    public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n----- Purine/Pyrimidine Analysis -----\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    calculateRatio(header, sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            calculateRatio(header, sequence);
        }

        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
class MeltingTempCalculator1 {
private:
    void calculateTm(const std::string& header, const std::string& sequence) const {
        int a = 0, t = 0, g = 0, c = 0;

        for (char base : sequence) {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': ++a; break;
                case 'T': ++t; break;
                case 'G': ++g; break;
                case 'C': ++c; break;
                default: break; // Skip N
            }
        }

        int total = a + t + g + c;
        if (total == 0) {
            std::cout << ">" << header << "\nNo valid bases found. Skipping...\n";
            return;
        }

        int tm = 2 * (a + t) + 4 * (g + c);

        std::cout << ">" << header << "\n";
        std::cout << "A: " << a << ", T: " << t << ", G: " << g << ", C: " << c << "\n";
        std::cout << "Melting Temperature (Tm): " << tm << "°C\n";
        std::cout << "-----------------------------------\n";
    }
    public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n----- Melting Temperature Analysis -----\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    calculateTm(header, sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            calculateTm(header, sequence);
        }

        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
class MeltingTempCalculator2 {
    private:
        struct ThermoParams {
            double dH; // kcal/mol
            double dS; // cal/(mol*K)
        };

        // SantaLucia 1998 parameters
        const std::unordered_map<std::string, ThermoParams> nnParams = {
            {"AA", {-7.9, -22.2}}, {"TT", {-7.9, -22.2}},
            {"AT", {-7.2, -20.4}}, {"TA", {-7.2, -21.3}},
            {"CA", {-8.5, -22.7}}, {"TG", {-8.5, -22.7}},
            {"GT", {-8.4, -22.4}}, {"AC", {-8.4, -22.4}},
            {"CT", {-7.8, -21.0}}, {"AG", {-7.8, -21.0}},
            {"GA", {-8.2, -22.2}}, {"TC", {-8.2, -22.2}},
            {"CG", {-10.6, -27.2}},{"GC", {-9.8, -24.4}},
            {"GG", {-8.0, -19.9}}, {"CC", {-8.0, -19.9}}
        };

        double R = 1.987; // cal/(K*mol)

        double calculateTmNN(const std::string& sequence, double strandConc = 5e-7, double NaConc = 0.05) const {
            if (sequence.size() < 2) return 0.0;

            double dH = 0.0; // kcal/mol
            double dS = 0.0; // cal/(mol*K)

            // initiation correction
            dH += 0.2;
            dS += -5.7;

            // nearest-neighbor summation
            for (size_t i = 0; i < sequence.size() - 1; i++) {
                std::string pair;
                pair += std::toupper(sequence[i]);
                pair += std::toupper(sequence[i+1]);

                auto it = nnParams.find(pair);
                if (it != nnParams.end()) {
                    dH += it->second.dH;
                    dS += it->second.dS;
                }
            }

            // symmetry correction if self-complementary
            bool symmetric = true;
            for (size_t i = 0; i < sequence.size() / 2; i++) {
                if (std::toupper(sequence[i]) != complement(sequence[sequence.size() - 1 - i])) {
                    symmetric = false;
                    break;
                }
            }
            if (symmetric) dS -= 1.4;

            // salt correction
            dS += 0.368 * (sequence.size() - 1) * std::log(NaConc);

            // calculate Tm
            double tm = (1000 * dH) / (dS + R * std::log(strandConc/4.0)) - 273.15;
            return tm;
        }

        static char complement(char base) {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'G': return 'C';
                case 'C': return 'G';
                default: return 'N';
            }
        }

        void calculateTm(const std::string& header, const std::string& sequence) const {
            double tm = calculateTmNN(sequence);

            int a=0,t=0,g=0,c=0;
            for (char base : sequence) {
                switch (std::toupper(static_cast<unsigned char>(base))) {
                    case 'A': ++a; break;
                    case 'T': ++t; break;
                    case 'G': ++g; break;
                    case 'C': ++c; break;
                }
            }

            std::cout << ">" << header << "\n";
            std::cout << "A: " << a << ", T: " << t << ", G: " << g << ", C: " << c << "\n";
            std::cout << "Melting Temperature (NN model): " << tm << " °C\n";
            std::cout << "-----------------------------------\n";
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file: " << filename << "\n";
                return;
            }

            std::string line, header, sequence;
            std::cout << "\n----- Melting Temperature Analysis (Nearest Neighbor Model) -----\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        calculateTm(header, sequence);
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }

            if (!sequence.empty()) {
                calculateTm(header, sequence);
            }

            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class ORFFinder {
    private:
        const std::vector<std::string> stopCodons = {"TAA", "TAG", "TGA"};

        bool isStopCodon(const std::string& codon) const {
            return std::find(stopCodons.begin(), stopCodons.end(), codon) != stopCodons.end();
        }

        std::vector<std::tuple<int,int,std::string>> findORFsInFrame(const std::string& sequence, int frame) const {
            std::vector<std::tuple<int,int,std::string>> orfs;
            for (size_t i = frame; i + 2 < sequence.size(); i += 3) {
                std::string codon = sequence.substr(i, 3);
                for (char &c : codon) c = std::toupper(c);

                if (codon == "ATG") { // find start codon
                    size_t start = i;
                    size_t j = i + 3;
                    for (; j + 2 < sequence.size(); j += 3) {
                        std::string nextCodon = sequence.substr(j, 3);
                        for (char &c : nextCodon) c = std::toupper(c);
                        if (isStopCodon(nextCodon)) {
                            orfs.emplace_back(start, j+2, sequence.substr(start, j+3-start));
                            break;
                        }
                    }
                    // go to next ORF
                    i = j;
                }
            }
            return orfs;
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }

            std::string line, header, sequence;
            std::cout << "\n----- ORF Finder -----\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        for (int frame = 0; frame < 3; ++frame) {
                            auto orfs = findORFsInFrame(sequence, frame);
                            std::cout << ">" << header << " - Frame " << frame+1 << "\n";
                            if (orfs.empty()) {
                                std::cout << "No ORFs found.\n";
                            } else {
                                for (auto &[start, end, seq] : orfs) {
                                    std::cout << "Start: " << start << ", End: " << end
                                            << ", Length (codons): " << (end-start+1)/3 << "\n";
                                }
                            }
                            std::cout << "-----------------------------------\n";
                        }
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }

            if (!sequence.empty()) {
                for (int frame = 0; frame < 3; ++frame) {
                    auto orfs = findORFsInFrame(sequence, frame);
                    std::cout << ">" << header << " - Frame " << frame+1 << "\n";
                    if (orfs.empty()) {
                        std::cout << "No ORFs found.\n";
                    } else {
                        for (auto &[start, end, seq] : orfs) {
                            std::cout << "Start: " << start << ", End: " << end
                                    << ", Length (codons): " << (end-start+1)/3 << "\n";
                        }
                    }
                    std::cout << "-----------------------------------\n";
                }
            }

            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class DNAcompToFile {
private:
    char Complement(char base) const {
        switch (std::toupper(static_cast<unsigned char>(base))) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default:  return 'N'; 
        }
    }

    std::string ComplementStrand(const std::string& sequence) const {
        std::string complement;
        complement.reserve(sequence.size());
        for (char base : sequence) {
            complement += Complement(base);
        }
        return complement;
    }

    // Wrap sequence into lines of 'width' characters
    void WriteWrappedSequence(std::ofstream& outFile, const std::string& sequence, std::size_t width = 70) const {
        for (std::size_t i = 0; i < sequence.size(); i += width) {
            outFile << sequence.substr(i, width) << "\n";
        }
    }
public:
    void FASTA_writer(const std::string& inputFilename, const std::string& outputFilename) const {
        std::ifstream fastaFile(inputFilename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open input file " << inputFilename << "\n";
            return;
        }

        std::ofstream outFile(outputFilename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open output file " << outputFilename << "\n";
            return;
        }

        std::string line;
        std::string header;
        std::string sequence;

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    std::string complement = ComplementStrand(sequence);
                    outFile << ">" << header << " (complement)\n";
                    WriteWrappedSequence(outFile, complement);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            std::string complement = ComplementStrand(sequence);
            outFile << ">" << header << " (complement)\n";
            WriteWrappedSequence(outFile, complement);
        }

        fastaFile.close();
        outFile.close();
        std::cout << "Complement sequences written to " << outputFilename << "\n";
    }
};
class ReverseComplementDNAToFile{
    private:
        char Complement(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return 'N'; 
            }
        }
       
        std::string ReverseComplementStrand(const std::string& sequence) const {
            std::string revComplement;
            revComplement.reserve(sequence.size());
            for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
                revComplement += Complement(*it);
            }
            return revComplement;
        }

        // Wrap sequence into lines of 'width' characters
        void WriteWrappedSequence(std::ofstream& outFile, const std::string& sequence, std::size_t width = 70) const {
            for (std::size_t i = 0; i < sequence.size(); i += width) {
                outFile << sequence.substr(i, width) << "\n";
            }
        }

    public:
        void FASTA_writer(const std::string& inputFilename, const std::string& outputFilename) const {
            std::ifstream fastaFile(inputFilename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open input file " << inputFilename << "\n";
                return;
            }

            std::ofstream outFile(outputFilename);
            if (!outFile.is_open()) {
                std::cerr << "Error: Unable to open output file " << outputFilename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string revComplement = ReverseComplementStrand(sequence);
                        outFile << ">" << header << " (complement)\n";
                        WriteWrappedSequence(outFile, revComplement);
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
            std::string complement = ReverseComplementStrand(sequence);
            outFile << ">" << header << " (complement)\n";
            WriteWrappedSequence(outFile, complement);
             }

            fastaFile.close();
            outFile.close();
            std::cout << "Complement sequences written to " << outputFilename << "\n";
    }
};
class TranscriptionToFile{
    private:
        char transRNA(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'A';
                case 'T': return 'U';
                case 'G': return 'G';
                case 'C': return 'C';
                default:  return 'N'; 
            }
        }
    
        std::string RNAStrand(const std::string& sequence) const {
            std::string RNAseq;
            RNAseq.reserve(sequence.size());
            for (char base : sequence) {
                RNAseq += transRNA(base);
            }
            return RNAseq;
        }
        
        void WriteWrappedSequence(std::ofstream& outFile, const std::string& sequence, std::size_t width = 70) const {
            for (std::size_t i = 0; i < sequence.size(); i += width) {
                outFile << sequence.substr(i, width) << "\n";
            }
        }

    public:
        void FASTA_writer(const std::string& inputFilename, const std::string& outputFilename) const {
            std::ifstream fastaFile(inputFilename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open input file " << inputFilename << "\n";
                return;
            }

            std::ofstream outFile(outputFilename);
            if (!outFile.is_open()) {
                std::cerr << "Error: Unable to open output file " << outputFilename << "\n";
                return;
            }

            std::string line;
            std::string header;
            std::string sequence;

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string RNAseq = RNAStrand(sequence);
                        outFile << ">" << header << " (RNA)\n";
                        WriteWrappedSequence(outFile, RNAseq);
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
            std::string RNAseq = RNAStrand(sequence);
            outFile << ">" << header << " (RNA)\n";
            WriteWrappedSequence(outFile, RNAseq);
             }

            fastaFile.close();
            outFile.close();
            std::cout << "RNA sequence written to " << outputFilename << "\n";
    }
};
class CodonUsageBias {
private:
    std::unordered_map<std::string, std::string> codonTable;
    std::unordered_map<std::string, int> codonCounts;
    std::unordered_map<std::string, std::vector<std::string>> aaToCodons;
    //static const std::unordered_map<std::string, char> codonTable; <-- remove it later

    void initCodonTable() {
        codonTable = {
            {"TTT","F"},{"TTC","F"},{"TTA","L"},{"TTG","L"},{"CTT","L"},{"CTC","L"},{"CTA","L"},{"CTG","L"},
            {"ATT","I"},{"ATC","I"},{"ATA","I"},{"ATG","M"},{"GTT","V"},{"GTC","V"},{"GTA","V"},{"GTG","V"},
            {"TCT","S"},{"TCC","S"},{"TCA","S"},{"TCG","S"},{"CCT","P"},{"CCC","P"},{"CCA","P"},{"CCG","P"},
            {"ACT","T"},{"ACC","T"},{"ACA","T"},{"ACG","T"},{"GCT","A"},{"GCC","A"},{"GCA","A"},{"GCG","A"},
            {"TAT","Y"},{"TAC","Y"},{"TAA","*"},{"TAG","*"},{"CAT","H"},{"CAC","H"},{"CAA","Q"},{"CAG","Q"},
            {"AAT","N"},{"AAC","N"},{"AAA","K"},{"AAG","K"},{"GAT","D"},{"GAC","D"},{"GAA","E"},{"GAG","E"},
            {"TGT","C"},{"TGC","C"},{"TGA","*"},{"TGG","W"},{"CGT","R"},{"CGC","R"},{"CGA","R"},{"CGG","R"},
            {"AGT","S"},{"AGC","S"},{"AGA","R"},{"AGG","R"},{"GGT","G"},{"GGC","G"},{"GGA","G"},{"GGG","G"}
        };

        for (auto &kv : codonTable) {
            if (kv.second != "*")
                aaToCodons[kv.second].push_back(kv.first);
        }
    }
    void countCodons(const std::string& sequence) {
        for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
            std::string codon = sequence.substr(i, 3);
            for (auto &c : codon) c = std::toupper(c);
            if (codonTable.count(codon)) {
                codonCounts[codon]++;
            }
        }
    }
    void printResults() const {
        std::cout << "\n----- Codon Usage Bias -----\n";
        std::cout << "Codon\tAA\tCount\tRSCU\n";

        // Compute RSCU for each amino acid
        for (auto &aaEntry : aaToCodons) {
            const std::string& aa = aaEntry.first;
            const auto& codons = aaEntry.second;

            double total = 0.0;
            for (const auto& codon : codons) {
                auto it = codonCounts.find(codon);
                if (it != codonCounts.end())
                    total += it->second;
            }

            for (const auto& codon : codons) {
                int count = codonCounts.count(codon) ? codonCounts.at(codon) : 0;
                double rscu = (total > 0) ? (count * codons.size() / total) : 0.0;
                std::cout << codon << "\t" << aa << "\t" << count << "\t" << rscu << "\n";
            }
        }
    }
public:
    CodonUsageBias() {
        initCodonTable();
    }

    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    countCodons(sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            countCodons(sequence);
        }

        printResults();
    }
};
class CodonUsageBiasCSV {
private:
    std::unordered_map<std::string, std::string> codonTable;
    std::unordered_map<std::string, int> codonCounts;
    std::unordered_map<std::string, std::vector<std::string>> aaToCodons;

    void initCodonTable() {
        codonTable = {
            {"TTT","F"},{"TTC","F"},{"TTA","L"},{"TTG","L"},{"CTT","L"},{"CTC","L"},{"CTA","L"},{"CTG","L"},
            {"ATT","I"},{"ATC","I"},{"ATA","I"},{"ATG","M"},{"GTT","V"},{"GTC","V"},{"GTA","V"},{"GTG","V"},
            {"TCT","S"},{"TCC","S"},{"TCA","S"},{"TCG","S"},{"CCT","P"},{"CCC","P"},{"CCA","P"},{"CCG","P"},
            {"ACT","T"},{"ACC","T"},{"ACA","T"},{"ACG","T"},{"GCT","A"},{"GCC","A"},{"GCA","A"},{"GCG","A"},
            {"TAT","Y"},{"TAC","Y"},{"TAA","*"},{"TAG","*"},{"CAT","H"},{"CAC","H"},{"CAA","Q"},{"CAG","Q"},
            {"AAT","N"},{"AAC","N"},{"AAA","K"},{"AAG","K"},{"GAT","D"},{"GAC","D"},{"GAA","E"},{"GAG","E"},
            {"TGT","C"},{"TGC","C"},{"TGA","*"},{"TGG","W"},{"CGT","R"},{"CGC","R"},{"CGA","R"},{"CGG","R"},
            {"AGT","S"},{"AGC","S"},{"AGA","R"},{"AGG","R"},{"GGT","G"},{"GGC","G"},{"GGA","G"},{"GGG","G"}
        };

        for (auto &kv : codonTable) {
            if (kv.second != "*")
                aaToCodons[kv.second].push_back(kv.first);
        }
    }

    void countCodons(const std::string& sequence) {
        for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
            std::string codon = sequence.substr(i, 3);
            for (auto &c : codon) c = std::toupper(c);
            if (codonTable.count(codon)) {
                codonCounts[codon]++;
            }
        }
    }

    void printResults(const std::string& outFile) const {
        std::ofstream csv(outFile);
        if (!csv.is_open()) {
            std::cerr << "Error: Could not open " << outFile << " for writing.\n";
            return;
        }

        std::cout << "\n----- Codon Usage Bias -----\n";
        std::cout << "Codon\tAA\tCount\tRSCU\n";

        csv << "Codon,AA,Count,RSCU\n"; // CSV header

        // Compute RSCU for each amino acid
        for (auto &aaEntry : aaToCodons) {
            const std::string& aa = aaEntry.first;
            const auto& codons = aaEntry.second;

            double total = 0.0;
            for (const auto& codon : codons) {
                auto it = codonCounts.find(codon);
                if (it != codonCounts.end())
                    total += it->second;
            }

            for (const auto& codon : codons) {
                int count = codonCounts.count(codon) ? codonCounts.at(codon) : 0;
                double rscu = (total > 0) ? (count * codons.size() / total) : 0.0;

                std::cout << codon << "\t" << aa << "\t" << count << "\t" 
                          << std::fixed << std::setprecision(2) << rscu << "\n";

                csv << codon << "," << aa << "," << count << "," 
                    << std::fixed << std::setprecision(2) << rscu << "\n";
            }
        }

        csv.close();
        std::cout << "\nResults saved to " << outFile << "\n";
    }

public:
    CodonUsageBiasCSV() {
        initCodonTable();
    }

    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    countCodons(sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            countCodons(sequence);
        }

        fastaFile.close();

        // export cub as csv
        std::string outFile = filename + "_cub.csv";
        printResults(outFile);
    }
};
class BasePairCounter {
private:
    int CountBases(const std::string& sequence) const {
        int count = 0;
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            if (upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G' || upper == 'N') {
                count++;
            }
        }
        return count;
    }

public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n----- Base Pair Counter -----\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    int bpCount = CountBases(sequence);
                    std::cout << ">" << header << "\nBase pairs: " << bpCount << "\n";
                    std::cout << "-----------------------------------\n";
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            int bpCount = CountBases(sequence);
            std::cout << ">" << header << "\nBase pairs: " << bpCount << "\n";
            std::cout << "-----------------------------------\n";
        }

        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
class MotifFinder {
private:
    std::string toUpper(const std::string& str) {
        std::string result = str;
        std::transform(result.begin(), result.end(), result.begin(), ::toupper);
        return result;
    }
    // reverse complement function
    std::string reverseComplement(const std::string& seq) {
        std::string rc;
        rc.reserve(seq.size());

        for (int i = seq.size() - 1; i >= 0; --i) {
            switch (toupper(seq[i])) {
                case 'A': rc.push_back('T'); break;
                case 'T': rc.push_back('A'); break;
                case 'C': rc.push_back('G'); break;
                case 'G': rc.push_back('C'); break;
                default:  rc.push_back('N'); break;
            }
        }
        return rc;
    }
    // Find all motifs
    std::vector<size_t> findPositions(const std::string& seq, const std::string& motif) {
        std::vector<size_t> positions;
        std::string seqUpper = toUpper(seq);
        std::string motifUpper = toUpper(motif);

        size_t pos = seqUpper.find(motifUpper, 0);
        while (pos != std::string::npos) {
            positions.push_back(pos + 1);  // 0-based indexing, +1 for biological indexing
            pos = seqUpper.find(motifUpper, pos + 1);
        }
        return positions;
    }
    void searchMotif(const std::string& header, const std::string& sequence, const std::string& motif) {
        std::cout << header << std::endl;
        // Forward motif
        auto forwardHits = findPositions(sequence, motif);
        // Reverse-complement motif
        std::string rc = reverseComplement(motif);
        auto rcHits = findPositions(sequence, rc);

        if (forwardHits.empty() && rcHits.empty()) {
            std::cout << "Motif not found (forward or reverse complement)." << std::endl;
        } else {
            if (!forwardHits.empty()) {
                std::cout << "Forward motif found at positions: ";
                for (auto p : forwardHits) std::cout << p << " ";
                std::cout << std::endl;
            }
            if (!rcHits.empty()) {
                std::cout << "Reverse-complement motif (" << rc << ") found at positions: ";
                for (auto p : rcHits) std::cout << p << " ";
                std::cout << std::endl;
            }
        }
        std::cout << "-----------------------------------" << std::endl;
    }
public:
    void FASTAloader(const std::string& filename, const std::string& motif) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        if (motif.empty()) {
            std::cerr << "Error: Motif string is empty." << std::endl;
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n----- Motif Search -----" << std::endl;

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    searchMotif(header, sequence, motif);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            searchMotif(header, sequence, motif);
        }

        std::cout << "----- Process completed. -----" << std::endl;
        fastaFile.close();
    }
};
class MolecularWeightCalculator {
private:
    const std::unordered_map<char, double> aaMasses = {
        {'A', 71.0788}, {'R', 156.1875}, {'N', 114.1039}, {'D', 115.0886},
        {'C', 103.1388}, {'E', 129.1155}, {'Q', 128.1307}, {'G', 57.0519},
        {'H', 137.1411}, {'I', 113.1594}, {'L', 113.1594}, {'K', 128.1741},
        {'M', 131.1986}, {'F', 147.1766}, {'P', 97.1167}, {'S', 87.0782},
        {'T', 101.1051}, {'W', 186.2132}, {'Y', 163.1760}, {'V', 99.1326},
        {'X', 110.0} 
    };
    std::string translateToAminoAcids(const std::string& sequence) {
        std::string protein;
        for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
            std::string codon = sequence.substr(i, 3);
            for (char& c : codon) c = std::toupper(c);
            if (codonTable.count(codon)) {
                protein += codonTable.at(codon);
            } else {
                protein += 'X';
            }
        }
        return protein;
    }
    double calculateMolecularWeight(const std::string& protein) {
        double totalMass = 0.0;
        for (char aa : protein) {
            char upperAA = std::toupper(aa);
            auto it = aaMasses.find(upperAA);
            totalMass += (it != aaMasses.end()) ? it->second : 110.0;
        }
        return totalMass;
    }
public:
    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        std::string line, header, sequence;
        std::cout << "----- Protein Molecular Weight (kDa) -----" << std::endl;
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    std::string protein = translateToAminoAcids(sequence);
                    double mw = calculateMolecularWeight(protein);
                    std::cout << header << "\n\nProtein: " << protein.size() << " AA" << std::endl;
                    std::cout << "Molecular Weight: " << std::fixed << std::setprecision(3) << (mw - 100)/1000 << " kDa" << std::endl;
                    std::cout << "-----------------------------------" << std::endl;
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            std::string protein = translateToAminoAcids(sequence);
            double mw = calculateMolecularWeight(protein);
            std::cout << header << "\n\nProtein: " << protein.size() << " AA" << std::endl;
            std::cout << "Molecular Weight: " << std::fixed << std::setprecision(3) << (mw - 100)/1000 << " kDa" << std::endl;
            std::cout << "-----------------------------------" << std::endl;
        }
        std::cout << "Process completed." << std::endl;
        fastaFile.close();
    }
};
class ProteinIsoelectricPoint {
private:
    // aprox pKa values
    const double pKa_N_term = 9.6;
    const double pKa_C_term = 2.4;
    const double pKa_K = 10.5;
    const double pKa_R = 12.5;
    const double pKa_H = 6.0;
    const double pKa_D = 3.9;
    const double pKa_E = 4.1;
    const double pKa_C = 8.3;
    const double pKa_Y = 10.1;

    std::string translateToAminoAcids(const std::string& sequence) {
        std::string protein;
        for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
            std::string codon = sequence.substr(i, 3);
            for (char& c : codon) c = std::toupper(static_cast<unsigned char>(c));
            if (codonTable.count(codon)) {
                protein += codonTable.at(codon);
            } else {
                protein += 'X'; 
            }
        }
        return protein;
    }

    struct Counts {
        int nterm_len = 0; 
        int cterm_len = 0;
        int D = 0, E = 0, C = 0, Y = 0, H = 0, K = 0, R = 0;
    };

    Counts countIonizable(const std::string& protein) const {
        Counts c;
        c.nterm_len = c.cterm_len = protein.empty() ? 0 : 1;
        for (char aa : protein) {
            switch (std::toupper(static_cast<unsigned char>(aa))) {
                case 'D': c.D++; break;
                case 'E': c.E++; break;
                case 'C': c.C++; break;
                case 'Y': c.Y++; break;
                case 'H': c.H++; break;
                case 'K': c.K++; break;
                case 'R': c.R++; break;
                default: break;
            }
        }
        return c;
    }

    double netChargeAtPH(const Counts& c, double pH) const {
        const double ten = 10.0;

        auto pos = [&](double pKa, int n) {
            if (n == 0) return 0.0;
            double term = 1.0 / (1.0 + std::pow(ten, pH - pKa));
            return n * term;  // +1 
        };

        auto neg = [&](double pKa, int n) {
            if (n == 0) return 0.0;
            double term = 1.0 / (1.0 + std::pow(ten, pKa - pH));
            return -n * term; // -1 
        };

        double charge = 0.0;
        // Termini
        if (c.nterm_len > 0)
            charge += pos(pKa_N_term, 1);
        if (c.cterm_len > 0)
            charge += neg(pKa_C_term, 1);

        // Side chains
        charge += pos(pKa_K, c.K);
        charge += pos(pKa_R, c.R);
        charge += pos(pKa_H, c.H);
        charge += neg(pKa_D, c.D);
        charge += neg(pKa_E, c.E);
        charge += neg(pKa_C, c.C);
        charge += neg(pKa_Y, c.Y);

        return charge;
    }

    double computePI(const std::string& protein) const {
        if (protein.empty()) return 0.0;
        Counts c = countIonizable(protein);

        double low = 0.0, high = 14.0;
        double mid = 7.0;
        // binary search
        for (int iter = 0; iter < 50; ++iter) { 
            mid = (low + high) / 2.0;
            double q = netChargeAtPH(c, mid);
            if (q > 0.0) {
                low = mid;
            } else {
                high = mid;
            }
        }
        return (low + high) / 2.0;
    }

public:
    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        std::string line, header, sequence;
        std::cout << "----- Protein Isoelectric Point (pI) -----" << std::endl;
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    std::string protein = translateToAminoAcids(sequence);
                    double pI = computePI(protein);
                    std::cout << header << "\nProtein: " << protein.size() << " AA" << std::endl;
                    std::cout << "Isoelectric point (pI): "
                              << std::fixed << std::setprecision(2) << pI << std::endl;
                    std::cout << "-----------------------------------" << std::endl;
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            std::string protein = translateToAminoAcids(sequence);
            double pI = computePI(protein);
            std::cout << header << "\nProtein: " << protein.size() << " AA" << std::endl;
            std::cout << "Isoelectric point (pI): "
                      << std::fixed << std::setprecision(2) << pI << std::endl;
            std::cout << "-----------------------------------" << std::endl;
        }
        std::cout << "Process completed." << std::endl;
        fastaFile.close();
    }
};
class ProteinExtinctionCoefficient {
private:
    std::string translateToAminoAcids(const std::string& sequence) {
        std::string protein;
        for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
            std::string codon = sequence.substr(i, 3);
            for (char& c : codon) c = std::toupper(static_cast<unsigned char>(c));
            if (codonTable.count(codon)) {
                protein += codonTable.at(codon);
            } else {
                protein += 'X';
            }
        }
        return protein;
    }

    double calculateExtinction(const std::string& protein) const {
        int C = 0, W = 0, Y = 0;
        for (char aa : protein) {
            char upper = std::toupper(static_cast<unsigned char>(aa));
            if (upper == 'C') C++;
            else if (upper == 'W') W++;
            else if (upper == 'Y') Y++;
        }

        // Gill & von Hippel method: singles + all pairwise interactions
        double epsilon = 0.0;
        
        // Single residues
        epsilon += C * 120.0;
        epsilon += W * 5500.0;
        epsilon += Y * 1490.0;
        
        // WW pairs
        epsilon += W * (W - 1) * 11000.0 / 2.0;
        // WY + YW pairs  
        epsilon += W * Y * 6990.0;
        // WC + CW pairs
        epsilon += W * C * 5620.0;
        // YY pairs
        epsilon += Y * (Y - 1) * 2980.0 / 2.0;
        // YC + CY pairs
        epsilon += Y * C * 2410.0;
        // CC pairs
        epsilon += C * (C - 1) * 120.0 / 2.0;

        return epsilon;
    }

public:
    void FASTA_loader(const std::string& filename) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }

        std::string line, header, sequence;
        std::cout << "----- Extinction Coefficient (ε280) -----" << std::endl;
        std::cout << "Units: M^-1 cm^-1 (Gill & von Hippel)" << std::endl;

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    std::string protein = translateToAminoAcids(sequence);
                    double epsilon = calculateExtinction(protein);
                    std::cout << header << "\nProtein: " << protein.size() << " AA" << std::endl;
                    std::cout << "ε280: " << std::fixed << std::setprecision(0) << epsilon 
                              << " M^-1 cm^-1" << std::endl;
                    std::cout << "-----------------------------------" << std::endl;
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            std::string protein = translateToAminoAcids(sequence);
            double epsilon = calculateExtinction(protein);
            std::cout << header << "\nProtein: " << protein.size() << " AA" << std::endl;
            std::cout << "ε280: " << std::fixed << std::setprecision(0) << epsilon 
                      << " M^-1 cm^-1" << std::endl;
            std::cout << "-----------------------------------" << std::endl;
        }
        std::cout << "Process completed." << std::endl;
        fastaFile.close();
    }
};
class cDNA_colour{
    private:
        char Complement(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default:  return 'N'; // Unknown base
            }
        }
        std::string ColorBase(char base) const {
        switch (base) {
            case 'A': return "\033[42mA\033[0m"; // Green background
            case 'T': return "\033[41mT\033[0m"; // Red background
            case 'G': return "\033[44mG\033[0m"; // Blue background
            case 'C': return "\033[40mC\033[0m"; // Black background
            default:  return "\033[47mN\033[0m"; // Gray background
            }
        }
        //DNA complement strand init
        std::string ComplementStrandColored(const std::string& sequence) const {
        std::string result;
        for (char base : sequence) {
            char comp = Complement(base);
            result += ColorBase(comp);
        }
        return result;
        }

    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }
            std::string line;
            std::string header;
            std::string sequence;
            std::cout << "\n-----------------------------------\n";
            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string complement = ComplementStrandColored(sequence);
                        std::cout << ">" << header << " (complement)\n" << complement << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
                std::string complement = ComplementStrandColored(sequence);
                std::cout << ">" << header << " (complement)\n" << complement << "\n";
            }
            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class Seq_colour{
    private:
        char Complement(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'A';
                case 'T': return 'T';
                case 'C': return 'C';
                case 'G': return 'G';
                default:  return 'N'; 
            }
        }
        std::string ColorBase(char base) const {
        switch (base) {
            case 'A': return "\033[42mA\033[0m"; // Green 
            case 'T': return "\033[41mT\033[0m"; // Red 
            case 'G': return "\033[44mG\033[0m"; // Blue 
            case 'C': return "\033[40mC\033[0m"; // Black 
            default:  return "\033[47mN\033[0m"; // Gray 
            }
        }
        //DNA complement strand init
        std::string ComplementStrandColored(const std::string& sequence) const {
        std::string result;
        for (char base : sequence) {
            char comp = Complement(base);
            result += ColorBase(comp);
        }
        return result;
        }
    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }
            std::string line;
            std::string header;
            std::string sequence;
            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::string complement = ComplementStrandColored(sequence);
                        std::cout << ">" << header << "\n" << complement << "\n\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
                std::string complement = ComplementStrandColored(sequence);
                std::cout << ">" << header << "\n" << complement << "\n";
            }
            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};
class HydrogenBondsCalc {
    private:
    void calculateHydrogenBonds(const std::string& header, const std::string& sequence) const{
        int sum_AT = 0, sum_GC = 0;
       
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            switch (upper) {
                case 'A':
                case 'T':
                    sum_AT++;
                    break;
                case 'C':
                case 'G':
                    sum_GC++;
                    break;
                default:
                    break; // Skip 'N' or unknown bases
            }
        }

        std::cout << ">" << header << "\n";
        std::cout << "A+T: " << sum_AT << "\n";
        std::cout << "G+C: " << sum_GC << "\n";

        if (sum_AT == 0) {
            std::cout << "Hydrogen bonds: Undefined (Hydrogen Bonds = 0)\n";
        } else {
            double hbonds = (sum_AT)*2 + sum_GC*3;
            std::cout << "Hydrogen Bonds: " << std::fixed << std::setprecision(0) << hbonds << "\n";
        }
        std::cout << "\n-----------------------------------\n";
    }
    public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n----- Hydrogen Bonds -----\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    calculateHydrogenBonds(header, sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            calculateHydrogenBonds(header, sequence);
        }

        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
class DNAambiguousStats {
    private:
        bool IsStandardBase(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A':
                case 'T':
                case 'G':
                case 'C':
                    return true;
                default:
                    return false;
            }
        }
        bool IsAmbiguousBase(char base) const {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'N': case 'R': case 'Y': case 'S': case 'W':
                case 'K': case 'M': case 'B': case 'D':
                case 'H': case 'V':
                    return true;
                default:
                    return false;
            }
        }
        void AnalyzeSequence(const std::string& sequence) const {
            std::map<char, size_t> counts;
            size_t total = 0;
            size_t ambiguousTotal = 0;

            for (char base : sequence) {
                char b = std::toupper(static_cast<unsigned char>(base));
                if (IsStandardBase(b) || IsAmbiguousBase(b)) {
                    counts[b]++;
                    total++;
                    if (IsAmbiguousBase(b)) ambiguousTotal++;
                }
            }
            std::cout << "Sequence length: " << total << " bp\n\n";
            std::cout << "Standard bases:\n";
            for (char b : {'A','T','G','C'}) {
                std::cout << b << ": " << counts[b] << "\n";
            }
            std::cout << "\nAmbiguous bases:\n";
            for (char b : {'N','R','Y','S','W','K','M','B','D','H','V'}) {
                std::cout << b << ": " << counts[b] << "\n";
            }
            double percent = total > 0 ? (ambiguousTotal * 100.0) / total : 0.0;
            std::cout << "\nTotal ambiguous bases: " << ambiguousTotal << "\n";
            std::cout << "Ambiguous percentage: "
                    << std::fixed << std::setprecision(2)
                    << percent << "%\n";
        }
    public:
        void FASTA_loader(const std::string& filename) const {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file " << filename << "\n";
                return;
            }
            std::string line;
            std::string header;
            std::string sequence;
            std::cout << "\n-----------------------------------\n";
            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        std::cout << ">" << header << " (ambiguous base statistics)\n";
                        AnalyzeSequence(sequence);
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    header = line.substr(1);
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty()) {
                std::cout << ">" << header << " (ambiguous base statistics)\n";
                AnalyzeSequence(sequence);
                std::cout << "\n-----------------------------------\n";
            }
            fastaFile.close();
            std::cout << "\nProcess completed.\n";
        }
};
class NeedlemanWunsch{
    private:
        struct AlignmentResult {
            std::string a1;
            std::string a2;
            int score;
        };

        bool isValidBase(char c) {
            c = std::toupper(c);
            return (c == 'A' || c == 'C' || c == 'G' ||
                    c == 'T' || c == 'N');
        }
        std::string sanitizeSequence(const std::string& line) {
            std::string clean;
            for (char c : line) {
                if (std::isspace(c)) continue;
                if (!isValidBase(c)) {
                    std::cerr << "Error: Invalid character '" << c << "' in FASTA sequence\n";
                    exit(1);
                }
                clean += std::toupper(c);
            }
            return clean;
        }   
        int match, mismatch, gap;
        int score(char a, char b){return(a == b) ? match : mismatch;}      
        AlignmentResult  align(const std::string& s1, const std::string& s2){
            int n = s1.size();
            int m = s2.size();
            std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));
            //algorithm init
            for (int i = 0; i <= n; i++)
                dp[i][0] = i * gap;
            for (int j = 0; j <= m; j++)
                dp[0][j] = j * gap;
            //fill matrix
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    int diag = dp[i - 1][j - 1] + score(s1[i - 1], s2[j - 1]);
                    int up   = dp[i - 1][j] + gap;
                    int left = dp[i][j - 1] + gap;
                    dp[i][j] = std::max({diag, up, left});
                }
            }
            std::string a1, a2;
            int i = n, j = m;
            while (i > 0 || j > 0){
                if(i > 0 && j > 0 && dp[i][j] == dp[i -1][j - 1] + score(s1[i - 1], s2[j - 1])){
                    a1 += s1[i - 1];
                    a2 += s2[j - 1];
                    i--; j--;
                }
                else if (i > 0 && dp[i][j] == dp[i - 1][j] + gap) {
                    a1 += s1[i - 1];
                    a2 += '-';
                    i--;
                }
                else {
                    a1 += '-';
                    a2 += s2[j - 1];
                    j--;
                }
            }
            std::reverse(a1.begin(), a1.end());
            std::reverse(a2.begin(), a2.end());
            AlignmentResult result;
            result.a1 = a1;
            result.a2 = a2;
            result.score = dp[n][m];
            
            std::cout << a1 << "\n";
            std::cout << a2 << "\n";
            std::cout << "\n\nAlignment score: " << result.score;
            return result;
        }   
        void calculateIdentity(const std::string& a1, const std::string& a2) {
            int identical = 0;
            int aligned = 0;

            for (size_t i = 0; i < a1.size(); i++) {
                if (a1[i] == '-' && a2[i] == '-') continue;
                aligned++;
                if (a1[i] == a2[i])
                    identical++;
            }
            double percentIdentity = (aligned > 0)
                ? (static_cast<double>(identical) / aligned) * 100.0
                : 0.0;
            std::cout << "\nPercentage identity: "
                    << std::fixed << std::setprecision(2)
                    << percentIdentity << "%\n";
        }
        int countGaps(const std::string& alignedSeq) {
            int gaps = 0;
            for (char c : alignedSeq) {
                if (c == '-') gaps++;
            }
            return gaps;
        }
    public:
        NeedlemanWunsch(int m = 1, int mm = -1, int g = -2)
            : match(m), mismatch(mm), gap(g) {}
        void FASTA_loader(const std::string& filename) {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file: " << filename << "\n";
                exit(1);
            }
            std::vector<std::string> headers;
            std::vector<std::string> sequences ;
            std::string line, seq;
            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;
                if (line[0] == '>') {
                    if (!seq.empty()) {
                        sequences.push_back(seq);
                        seq.clear();
                    }
                    headers.push_back(line.substr(1));
                } else {
                    seq += sanitizeSequence(line);
            }
        }
        if (!seq.empty())
            sequences.push_back(seq);
        fastaFile.close();
        if (sequences.size() != 2) {
            std::cerr << "Error: FASTA file must contain exactly two sequences\n";
            exit(1);
        }
        std::cout << "\nNeedleman-Wunsch Pairwise Sequence Alignment\n";
        std::cout << "-----------------------------------\n";
        std::cout << "Sequence 1: " << headers[0] << "\n";
        std::cout << "Sequence 2: " << headers[1] << "\n";
        std::cout << "-----------------------------------\n\n";
        AlignmentResult res = align(sequences[0], sequences[1]);
        calculateIdentity(res.a1, res.a2);
        int gapsSeq1 = countGaps(res.a1);
        int gapsSeq2 = countGaps(res.a2);
        int GapSum = gapsSeq1 + gapsSeq2;
        std::cout << "Gaps: " << GapSum;
        std::cout << "\n-----------------------------------\n";
        std::cout << "Process completed.\n";
    }
};
class SmithWaterman {
    private:
        struct AlignmentResult {
            std::string a1;
            std::string a2;
            int score;
        };
        bool isValidBase(char c) {
            c = std::toupper(c);
            return (c == 'A' || c == 'C' || c == 'G' ||
                    c == 'T' || c == 'N');
        }
        std::string sanitizeSequence(const std::string& line) {
            std::string clean;
            for (char c : line) {
                if (std::isspace(c)) continue;
                if (!isValidBase(c)) {
                    std::cerr << "Error: Invalid character '" << c << "' in FASTA sequence\n";
                    exit(1);
                }
                clean += std::toupper(c);
            }
            return clean;
        }

        int match, mismatch, gap;

        int score(char a, char b) {
            return (a == b) ? match : mismatch;
        }

        AlignmentResult align(const std::string& s1, const std::string& s2) {
            int n = s1.size();
            int m = s2.size();

            std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

            int maxScore = 0;
            int maxI = 0, maxJ = 0;

            // Fill matrix
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    int diag = dp[i - 1][j - 1] + score(s1[i - 1], s2[j - 1]);
                    int up   = dp[i - 1][j] + gap;
                    int left = dp[i][j - 1] + gap;

                    dp[i][j] = std::max({0, diag, up, left});

                    if (dp[i][j] > maxScore) {
                        maxScore = dp[i][j];
                        maxI = i;
                        maxJ = j;
                    }
                }
            }
            // Traceback from max cell
            std::string a1, a2;
            int i = maxI, j = maxJ;

            while (i > 0 && j > 0 && dp[i][j] != 0) {
                if (dp[i][j] == dp[i - 1][j - 1] + score(s1[i - 1], s2[j - 1])) {
                    a1 += s1[i - 1];
                    a2 += s2[j - 1];
                    i--; j--;
                }
                else if (dp[i][j] == dp[i - 1][j] + gap) {
                    a1 += s1[i - 1];
                    a2 += '-';
                    i--;
                }
                else {
                    a1 += '-';
                    a2 += s2[j - 1];
                    j--;
                }
            }

            std::reverse(a1.begin(), a1.end());
            std::reverse(a2.begin(), a2.end());

            AlignmentResult result;
            result.a1 = a1;
            result.a2 = a2;
            result.score = maxScore;

            std::cout << a1 << "\n";
            std::cout << a2 << "\n";
            std::cout << "\n\nAlignment score: " << result.score;

            return result;
        }

        void calculateIdentity(const std::string& a1, const std::string& a2) {
            int identical = 0;
            int aligned = 0;

            for (size_t i = 0; i < a1.size(); i++) {
                if (a1[i] == '-' && a2[i] == '-') continue;
                aligned++;
                if (a1[i] == a2[i])
                    identical++;
            }

            double percentIdentity = (aligned > 0)
                ? (static_cast<double>(identical) / aligned) * 100.0
                : 0.0;

            std::cout << "\nPercentage identity: "
                    << std::fixed << std::setprecision(2)
                    << percentIdentity << "%\n";
        }

        int countGaps(const std::string& alignedSeq) {
            int gaps = 0;
            for (char c : alignedSeq)
                if (c == '-') gaps++;
            return gaps;
        }

    public:
        SmithWaterman(int m = 1, int mm = -1, int g = -2)
            : match(m), mismatch(mm), gap(g) {}

        void FASTA_loader(const std::string& filename) {
            std::ifstream fastaFile(filename);
            if (!fastaFile.is_open()) {
                std::cerr << "Error: Unable to open file: " << filename << "\n";
                exit(1);
            }

            std::vector<std::string> headers;
            std::vector<std::string> sequences;
            std::string line, seq;

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!seq.empty()) {
                        sequences.push_back(seq);
                        seq.clear();
                    }
                    headers.push_back(line.substr(1));
                }
                else {
                    seq += sanitizeSequence(line);
                }
            }

            if (!seq.empty())
                sequences.push_back(seq);

            fastaFile.close();

            if (sequences.size() != 2) {
                std::cerr << "Error: FASTA file must contain exactly two sequences\n";
                exit(1);
            }

            std::cout << "\nSmith-Waterman Pairwise Sequence Alignment\n";
            std::cout << "-----------------------------------\n";
            std::cout << "Sequence 1: " << headers[0] << "\n";
            std::cout << "Sequence 2: " << headers[1] << "\n";
            std::cout << "-----------------------------------\n\n";

            AlignmentResult res = align(sequences[0], sequences[1]);
            calculateIdentity(res.a1, res.a2);

            int gapsSeq1 = countGaps(res.a1);
            int gapsSeq2 = countGaps(res.a2);

            std::cout << "Gaps: " << gapsSeq1 + gapsSeq2 << "\n";
            std::cout << "-----------------------------------\n";
            std::cout << "Process completed.\n";
        }
};
class AssemblyCoverage{
private:
    int CountBases(const std::string& sequence) const {
        int count = 0;
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            if (upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G' || upper == 'N') {
                count++;
            }
        }
        return count;
    }

public:
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << "\n";
            return;
        }
        size_t GenomeSize;
        std::cout << "Enter the genome size: \n";
        std::cin >> GenomeSize;
        std::string line, header, sequence;
        std::cout << "\n----- Coverage(x) -----\n";
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    int bpCount = CountBases(sequence);
                    size_t coveragex = bpCount / GenomeSize;
                    std::cout << ">" << header << "\nCoverage: " << coveragex << "x\n";
                    std::cout << "-----------------------------------\n";
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty()) {
            int bpCount = CountBases(sequence);
            size_t coveragex = bpCount / GenomeSize;
            std::cout << ">" << header << "\nCoverage: " << coveragex << "x\n";
            std::cout << "-----------------------------------\n";
        }
        std::cout << "Process completed.\n";
        fastaFile.close();
    }
};
// Custom Pipelines bellow:
class Nucleostats {
     /*The old Pipeline 2. This is a pipeline that contains the following classes and functions:
    Classes:PurinePyrimidineRatioAnalyzer, MeltingTempCalculator2, GCCalc, BasePairCounter.*/
    private:
    // multithreaded fix later
    struct SeqResult {
        std::string header;
        int purines;
        int pyrimidines;
        int bpCount;
        double tm;
        double GCContent;
    };
    void ppRatio(const std::string& header, const std::string& sequence) const{
        int purines = 0, pyrimidines = 0;
       
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            switch (upper) {
                case 'A':
                case 'G':
                    purines++;
                    break;
                case 'C':
                case 'T':
                    pyrimidines++;
                    break;
                default:
                    break; // Skip N
            }
        }
        std::cout << "\n------------NucleoStats-----------\n";
        std::cout << ">" << header << "\n";
        std::cout << "-----------------------------------\n";
        std::cout << "Purines: " << purines << "\n";
        std::cout << "Pyrimidines: " << pyrimidines << "\n";
        if (pyrimidines == 0) {
            std::cout << "Purine/Pyrimidine Ratio: Undefined (pyrimidines = 0)\n";
        } else {
            double ratio = static_cast<double>(purines) / pyrimidines;
            std::cout << "Purine/Pyrimidine Ratio: " << std::fixed << std::setprecision(3) << ratio << "\n";
        }
    }

    int CountBases(const std::string& sequence) const {
        int count = 0;
        for (char base : sequence) {
            char upper = std::toupper(static_cast<unsigned char>(base));
            if (upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G' || upper == 'N') {
                count++;
            }
        }
        return count;
    }
        struct ThermoParams {
            double dH; // kcal/mol
            double dS; // cal/(mol*K)
        };
        // SantaLucia 1998 parameters
        const std::unordered_map<std::string, ThermoParams> nnParams = {
            {"AA", {-7.9, -22.2}}, {"TT", {-7.9, -22.2}},
            {"AT", {-7.2, -20.4}}, {"TA", {-7.2, -21.3}},
            {"CA", {-8.5, -22.7}}, {"TG", {-8.5, -22.7}},
            {"GT", {-8.4, -22.4}}, {"AC", {-8.4, -22.4}},
            {"CT", {-7.8, -21.0}}, {"AG", {-7.8, -21.0}},
            {"GA", {-8.2, -22.2}}, {"TC", {-8.2, -22.2}},
            {"CG", {-10.6, -27.2}},{"GC", {-9.8, -24.4}},
            {"GG", {-8.0, -19.9}}, {"CC", {-8.0, -19.9}}
        };

        double R = 1.987; // cal/(K*mol)

        double calculateTmNN(const std::string& sequence, double strandConc = 5e-7, double NaConc = 0.05) const {
            if (sequence.size() < 2) return 0.0;

            double dH = 0.0; // kcal/mol
            double dS = 0.0; // cal/(mol*K)

            // initiation correction
            dH += 0.2;
            dS += -5.7;

            // nearest-neighbor summation
            for (size_t i = 0; i < sequence.size() - 1; i++) {
                std::string pair;
                pair += std::toupper(sequence[i]);
                pair += std::toupper(sequence[i+1]);

                auto it = nnParams.find(pair);
                if (it != nnParams.end()) {
                    dH += it->second.dH;
                    dS += it->second.dS;
                }
            }
            // symmetry correction if self-complementary
            bool symmetric = true;
            for (size_t i = 0; i < sequence.size() / 2; i++) {
                if (std::toupper(sequence[i]) != complement(sequence[sequence.size() - 1 - i])) {
                    symmetric = false;
                    break;
                }
            }
            if (symmetric) dS -= 1.4;

            // salt correction
            dS += 0.368 * (sequence.size() - 1) * std::log(NaConc);

            // calculate Tm
            double tm = (1000 * dH) / (dS + R * std::log(strandConc/4.0)) - 273.15;
            return tm;
        }

        static char complement(char base) {
            switch (std::toupper(static_cast<unsigned char>(base))) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'G': return 'C';
                case 'C': return 'G';
                default: return 'N';
            }
        }

        void TmCalc(const std::string& /*header*/, const std::string& sequence) const {
            double tm = calculateTmNN(sequence);

            int a=0,t=0,g=0,c=0;
            for (char base : sequence) {
                switch (std::toupper(static_cast<unsigned char>(base))) {
                    case 'A': ++a; break;
                    case 'T': ++t; break;
                    case 'G': ++g; break;
                    case 'C': ++c; break;
                }
            }
            std::cout << "A: " << a << ", T: " << t << ", G: " << g << ", C: " << c << "\n";
            std::cout << "Melting Temperature (NN model): " << tm << " °C\n";
        }
    
    double GCCon(const std::string& sequence) const {
            int gcCount = 0;
            int validBases = 0;

            for (char base : sequence) {
                char upperBase = std::toupper(base);
                if (upperBase == 'G' || upperBase == 'C') {
                    gcCount++;
                    validBases++;
                } else if (upperBase == 'A' || upperBase == 'T') {
                    validBases++;
                }   
            }
            if (validBases == 0) return 0.0;
            return (static_cast<double>(gcCount) / validBases) * 100.0;
    }
    public:
    void FASTA_loader(const std::string& filename, size_t numThreads = 4) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }
        std::cout << "\n------------ NucleoStats -----------\n";
        std::vector<std::pair<std::string, std::string>> sequences;
        std::string line, header, sequence;
        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!sequence.empty()) {
                    sequences.emplace_back(header, sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty())
            sequences.emplace_back(header, sequence);
        fastaFile.close();
        std::vector<SeqResult> results(sequences.size());
        std::mutex indexMutex;
        size_t index = 0;
        auto worker = [&]() {
            while (true) {
                size_t i;
                {
                    std::lock_guard<std::mutex> lock(indexMutex);
                    if (index >= sequences.size()) return;
                    i = index++;
                }
                const auto& [hdr, seq] = sequences[i];
                SeqResult r;
                r.header = hdr;
                for (char base : seq) {
                    switch (std::toupper(static_cast<unsigned char>(base))) {
                        case 'A':
                        case 'G': r.purines++; break;
                        case 'C':
                        case 'T': r.pyrimidines++; break;
                    }
                }
                r.bpCount = CountBases(seq);
                r.tm = calculateTmNN(seq);
                r.GCContent = GCCon(seq);
                results[i] = r;
            }
        };
        std::vector<std::thread> threads;
        numThreads = std::min(numThreads, sequences.size());
        for (size_t t = 0; t < numThreads; ++t)
            threads.emplace_back(worker);
        for (auto& t : threads)
            t.join();
        for (const auto& r : results) {
            std::cout << r.header << "\n";
            std::cout << "-----------------------------------\n";
            std::cout << "Base pairs: " << r.bpCount << "\n";
            std::cout << "GC Content = "
                    << std::fixed << std::setprecision(2)
                    << r.GCContent << "%\n";
            std::cout << "Melting Temperature (NN model): "
                    << std::fixed << std::setprecision(2)
                    << r.tm << " °C\n";
            std::cout << "Purines: " << r.purines << "\n";
            std::cout << "Pyrimidines: " << r.pyrimidines << "\n";
            if (r.pyrimidines == 0)
                std::cout << "Purine/Pyrimidine Ratio: Undefined\n";
            else
                std::cout << "Purine/Pyrimidine Ratio: "
                        << std::fixed << std::setprecision(3)
                        << static_cast<double>(r.purines) / r.pyrimidines << "\n";   
            std::cout << "-----------------------------------\n";
        }
        std::cout << "Process completed.\n";
}

};
class Proteostats{
    /*The old Pipeline 3. This is a pipeline for structual analysis. It contains the following classes and functions:
    ProteinIsoelectricPoint, ProteinExtinctionCoefficient, MolecularWeightCalculator.
    */
private:
    // aprox pKa values
    const double pKa_N_term = 9.6;
    const double pKa_C_term = 2.4;
    const double pKa_K = 10.5;
    const double pKa_R = 12.5;
    const double pKa_H = 6.0;
    const double pKa_D = 3.9;
    const double pKa_E = 4.1;
    const double pKa_C = 8.3;
    const double pKa_Y = 10.1;

    struct SeqResult {
            std::string header;
            std::string protein;
            double mw;
            double mw_calibrated;
            double pI;
            double epsilon;
        };
    const std::unordered_map<char, double> aaMasses = {
            {'A', 71.0788}, {'R', 156.1875}, {'N', 114.1039}, {'D', 115.0886},
            {'C', 103.1388}, {'E', 129.1155}, {'Q', 128.1307}, {'G', 57.0519},
            {'H', 137.1411}, {'I', 113.1594}, {'L', 113.1594}, {'K', 128.1741},
            {'M', 131.1986}, {'F', 147.1766}, {'P', 97.1167}, {'S', 87.0782},
            {'T', 101.1051}, {'W', 186.2132}, {'Y', 163.1760}, {'V', 99.1326},
            {'X', 110.0} 
        };
        std::string translateToAminoAcids(const std::string& sequence) {
            std::string protein;
            for (size_t i = 0; i + 2 < sequence.size(); i += 3) {
                std::string codon = sequence.substr(i, 3);
                for (char& c : codon) c = std::toupper(c);
                if (codonTable.count(codon)) {
                    protein += codonTable.at(codon);
                } else {
                    protein += 'X';
                }
            }
            return protein;
        }
        double calculateMolecularWeight(const std::string& protein) {
            double totalMass = 0.0;
            for (char aa : protein) {
                char upperAA = std::toupper(aa);
                auto it = aaMasses.find(upperAA);
                totalMass += (it != aaMasses.end()) ? it->second : 110.0;
            }
            return totalMass;
        }
        struct Counts {
        int nterm_len = 0; 
        int cterm_len = 0;
        int D = 0, E = 0, C = 0, Y = 0, H = 0, K = 0, R = 0;
    };
    Counts countIonizable(const std::string& protein) const {
        Counts c;
        c.nterm_len = c.cterm_len = protein.empty() ? 0 : 1;
        for (char aa : protein) {
            switch (std::toupper(static_cast<unsigned char>(aa))) {
                case 'D': c.D++; break;
                case 'E': c.E++; break;
                case 'C': c.C++; break;
                case 'Y': c.Y++; break;
                case 'H': c.H++; break;
                case 'K': c.K++; break;
                case 'R': c.R++; break;
                default: break;
            }
        }
        return c;
    }
    double netChargeAtPH(const Counts& c, double pH) const {
        const double ten = 10.0;
        auto pos = [&](double pKa, int n) {
            if (n == 0) return 0.0;
            double term = 1.0 / (1.0 + std::pow(ten, pH - pKa));
            return n * term;  // +1 
        };
        auto neg = [&](double pKa, int n) {
            if (n == 0) return 0.0;
            double term = 1.0 / (1.0 + std::pow(ten, pKa - pH));
            return -n * term; // -1 
        };
        double charge = 0.0;
        // Termini
        if (c.nterm_len > 0)
            charge += pos(pKa_N_term, 1);
        if (c.cterm_len > 0)
            charge += neg(pKa_C_term, 1);
        // Side chains
        charge += pos(pKa_K, c.K);
        charge += pos(pKa_R, c.R);
        charge += pos(pKa_H, c.H);
        charge += neg(pKa_D, c.D);
        charge += neg(pKa_E, c.E);
        charge += neg(pKa_C, c.C);
        charge += neg(pKa_Y, c.Y);
        return charge;
    }
    double computePI(const std::string& protein) const {
        if (protein.empty()) return 0.0;
        Counts c = countIonizable(protein);
        double low = 0.0, high = 14.0;
        double mid = 7.0;
        // binary search
        for (int iter = 0; iter < 50; ++iter) { 
            mid = (low + high) / 2.0;
            double q = netChargeAtPH(c, mid);
            if (q > 0.0) {
                low = mid;
            } else {
                high = mid;
            }
        }
        return (low + high) / 2.0;
    }
    double calculateExtinction(const std::string& protein) const {
            int C = 0, W = 0, Y = 0;
            for (char aa : protein) {
                char upper = std::toupper(static_cast<unsigned char>(aa));
                if (upper == 'C') C++;
                else if (upper == 'W') W++;
                else if (upper == 'Y') Y++;
            }
            // Gill & von Hippel method: singles + all pairwise interactions
            double epsilon = 0.0;           
            // Single residues
            epsilon += C * 120.0;
            epsilon += W * 5500.0;
            epsilon += Y * 1490.0;           
            // WW pairs
            epsilon += W * (W - 1) * 11000.0 / 2.0;
            // WY + YW pairs  
            epsilon += W * Y * 6990.0;
            // WC + CW pairs
            epsilon += W * C * 5620.0;
            // YY pairs
            epsilon += Y * (Y - 1) * 2980.0 / 2.0;
            // YC + CY pairs
            epsilon += Y * C * 2410.0;
            // CC pairs
            epsilon += C * (C - 1) * 120.0 / 2.0;
            return epsilon;
        }        
public:
    void FASTA_loader(const std::string& filename, size_t numThreads = 4) {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        std::cout << "\n--------- Proteostats ---------\n";
        std::vector<std::pair<std::string, std::string>> sequences;
        std::string line, header, sequence;

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    sequences.emplace_back(header, sequence);
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }
        if (!sequence.empty())
            sequences.emplace_back(header, sequence);

        fastaFile.close();

        std::vector<SeqResult> results(sequences.size());
        std::mutex indexMutex;
        size_t index = 0;

        auto worker = [&]() {
            while (true) {
                size_t i;
                {
                    std::lock_guard<std::mutex> lock(indexMutex);
                    if (index >= sequences.size()) return;
                    i = index++;
                }
                const auto& [hdr, dna] = sequences[i];
                std::string protein = translateToAminoAcids(dna);

                double mw = calculateMolecularWeight(protein);
                double mw_calibrated = (mw - 100.0) / 1000.0;
                double pI = computePI(protein);
                double epsilon = calculateExtinction(protein);
                results[i] = {hdr, protein, mw, mw_calibrated, pI, epsilon};
            }
        };
        std::vector<std::thread> threads;
        numThreads = std::min(numThreads, sequences.size());
        for (size_t t = 0; t < numThreads; ++t)
            threads.emplace_back(worker);
        for (auto& t : threads)
            t.join();
        for (const auto& r : results) {
            std::cout << r.header << "\n";
            std::cout << "-----------------------------------\n";
            std::cout << "Protein: " << r.protein.size() << " AA\n";
            std::cout << "Molecular Weight: " << std::fixed << std::setprecision(3)
                    << r.mw_calibrated << " kDa\n";
            std::cout << "Isoelectric point (pI): " << std::fixed
                    << std::setprecision(2) << r.pI << "\n";
            std::cout << "Extinction Coefficient(ε280): "
                    << std::fixed << std::setprecision(0)
                    << r.epsilon << " M^-1 cm^-1\n";
            std::cout << "-----------------------------------\n";
        }
        std::cout << "Process completed.\n";
    }    
};
class AssemblyStats{
    private:
        struct AssemblyResult {
            size_t numContigs = 0;
            size_t totalLength = 0;
            double gcPercent = 0.0;
            size_t n50 = 0;
            size_t l50 = 0;
            size_t longestContig = 0;
            double meanLength = 0.0;
            int Coverage = 0.0;
        };
        std::string cleanSequence(const std::string& seq) const {
            std::string cleaned;
            cleaned.reserve(seq.size());
            for (char c : seq) {
                if (!std::isspace(static_cast<unsigned char>(c)))
                    cleaned += std::toupper(static_cast<unsigned char>(c));
            }
            return cleaned;
        }
        size_t countGC(const std::string& seq) const {
            size_t gc = 0;
            for (char c : seq)
                if (c == 'G' || c == 'C') gc++;
            return gc;
        }
        void computeN50(std::vector<size_t>& lengths, size_t& n50, size_t& l50) const {
            std::sort(lengths.begin(), lengths.end(), std::greater<size_t>());
            size_t total = 0;
            for (size_t len : lengths)
                total += len;
            size_t cumulative = 0;
            l50 = 0;
            n50 = 0;
            for (size_t len : lengths) {
                cumulative += len;
                l50++;
                if (cumulative >= total / 2) {
                    n50 = len;
                    break;
                }
            }
        }
    public:
        void FASTA_loader(const std::string& filename, size_t numThreads = 4) {
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Error: Cannot open file " << filename << std::endl;
                return;
            }

            std::vector<std::string> sequences;
            std::string line, sequence;

            while (std::getline(file, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        sequences.emplace_back(cleanSequence(sequence));
                        sequence.clear();
                    }
                } else {
                    sequence += line;
                }
            }
            if (!sequence.empty())
                sequences.emplace_back(cleanSequence(sequence));
            file.close();
            if (sequences.empty()) {
                std::cout << "No contigs found.\n";
                return;
            }
            std::mutex mtx;
            AssemblyResult result;
            std::vector<size_t> lengths;
            size_t totalGC = 0;
            size_t index = 0;
            auto worker = [&]() {
                while (true) {
                    size_t i;

                    {
                        std::lock_guard<std::mutex> lock(mtx);
                        if (index >= sequences.size()) return;
                        i = index++;
                    }

                    const std::string& seq = sequences[i];
                    size_t len = seq.size();
                    size_t gc = countGC(seq);

                    {
                        std::lock_guard<std::mutex> lock(mtx);
                        lengths.push_back(len);
                        result.totalLength += len;
                        totalGC += gc;
                        result.longestContig =
                            std::max(result.longestContig, len);
                    }
                }
    };

        numThreads = std::min(numThreads, sequences.size());

        std::vector<std::thread> threads;
        for (size_t t = 0; t < numThreads; ++t)
            threads.emplace_back(worker);

        for (auto& t : threads)
            t.join();

        result.numContigs = sequences.size();
        result.meanLength =
            static_cast<double>(result.totalLength) / result.numContigs;

        result.gcPercent =
            (result.totalLength > 0)
                ? (100.0 * static_cast<double>(totalGC) / result.totalLength)
                : 0.0;

        computeN50(lengths, result.n50, result.l50);
        size_t GenomeSize;
        std::cout << "Enter the genome size: \n";
        std::cin >> GenomeSize;
        size_t Coverage = result.totalLength / GenomeSize;

        //Results
        std::cout << "\n--------- Assembly Stats ---------\n";
        std::cout << "Number of Contigs: " << result.numContigs << "\n";
        std::cout << "Total Length: " << result.totalLength << " bp (" << std::fixed << std::setprecision(2) << static_cast<double>(result.totalLength) / 1e6 << " Mb)\n";
        std::cout << "GC%: " << result.gcPercent << "% (AT%: " << 100.0 - result.gcPercent << "%)\n";
        std::cout << "N50: " << result.n50 << " bp (" << static_cast<double>(result.n50) / 1e6 << " Mb)\n";
        std::cout << "L50: " << result.l50 << "\n";
        std::cout << "Coverage(x): " << Coverage << "x\n";                        ;
        std::cout << "Longest Contig: " << result.longestContig << " bp (" << static_cast<double>(result.longestContig) / 1e6 << " Mb)\n";
        std::cout << "Mean Contig Length: " << result.meanLength << " bp (" << result.meanLength / 1e6 << " Mb)\n";
        std::cout << "----------------------------------\n";
        std::cout << "Process completed.\n\n";        
    }
};
// Main Function 
int main(int argc, char* argv[]){
    if (argc != 3){
        title();
        message();
        return 1;
    }
    title();
    std::string filename = argv[2];
    std::string function = argv[1]; 
    // File verifier
    FileVerifier(filename);
    // Start timer
    auto start_timer = std::chrono::high_resolution_clock::now();
    // Main dispatch map
    std::unordered_map<std::string, std::function<void()>> dispatch {
        {"-gc", [&]() {
            GCCalc GCcalculator;
            GCcalculator.FASTA_loader(filename);
        }},
        {"-nc", [&]() {
            CodonNumber codoncounter;
            codoncounter.FASTA_loader(filename);
        }},
        {"-c", [&]() {
            DNAcomplementary DNAcomp;
            DNAcomp.FASTA_loader(filename);
        }},
        {"-rc", [&]() {
            ReverseComplementDNA revDNA;
            revDNA.FASTA_loader(filename);
        }},
        {"-t", [&]() {
            Transcription x;
            x.FASTA_loader(filename);
        }},
        {"-help", [&]() { helpme(); }},
        {"-p", [&]() {
            ProteinChain protein;
            protein.FASTA_loader(filename);
        }},
        {"-ss", [&]() {
            FASTAChromosomeSeparator splitter;
            splitter.FASTA_loader(filename);
        }},
        {"-sh", [&]() {
            FASTASequenceHeader header;
            header.FASTA_loader(filename);
        }},
        {"-tr", [&]() {
            DNATrimmer trim;
            int start_position, end_position;
            std::cout << "Enter the starting position:";
            std::cin >> start_position;
            std::cout << "Enter the end position:";
            std::cin >> end_position;
            trim.FASTA_loader(filename, start_position, end_position);
        }},
        {"-pp", [&]() {
            PurinePyrimidineRatioAnalyzer ppanalyzer;
            ppanalyzer.FASTA_loader(filename);
        }},
        {"-mt1", [&]() {
            MeltingTempCalculator1 mt1;
            mt1.FASTA_loader(filename);
        }},
        {"-mt2", [&]() {
            MeltingTempCalculator2 mt2;
            mt2.FASTA_loader(filename);
        }},
        {"-cc", [&]() {
            cDNA_colour dnacolour;
            dnacolour.FASTA_loader(filename);
        }},
        {"-nucleo", [&]() {
            //The old pip2
            int num_threads;
            Nucleostats pipeline2;
            
            std::cout << "Select number of threads(default is 2): ";
            std::cin >> num_threads;

            if(num_threads < 2){
                std::cout << "\nWarning!This is a multithreaded function. Cannot be used with less than 2 threads.\n";
                std::cout << "Using 2 threads. . .\n\n\n";
                num_threads = 2;
            }
            pipeline2.FASTA_loader(filename, num_threads);
        }},
        {"-orf", [&]() {
            ORFFinder orffinder;
            orffinder.FASTA_loader(filename);
        }},
        {"-cw", [&]() {
            std::string outputFile;
            std::cout << "Enter output filename: ";
            std::cin >> outputFile;
            DNAcompToFile cdna_fasta;
            cdna_fasta.FASTA_writer(filename, outputFile);
        }},
        {"-rcw", [&]() {
            std::string outputFile;
            std::cout << "Enter output filename: ";
            std::cin >> outputFile;
            ReverseComplementDNAToFile rev_cdna_fasta;
            rev_cdna_fasta.FASTA_writer(filename, outputFile);
        }},
        {"-tw", [&]() {
            std::string outputFile;
            std::cout << "Enter output filename: ";
            std::cin >> outputFile;
            TranscriptionToFile writeRNA;
            writeRNA.FASTA_writer(filename, outputFile);
        }},
        {"-mw", [&]() {
            MolecularWeightCalculator mwcalc;
            mwcalc.FASTA_loader(filename);
        }},
        {"-pi", [&]() {
            ProteinIsoelectricPoint picalc;
            picalc.FASTA_loader(filename);
        }},
        {"-ec", [&]() {
            ProteinExtinctionCoefficient eccalc;
            eccalc.FASTA_loader(filename);
        }},
        {"-cub", [&]() {
            CodonUsageBias cub;
            cub.FASTA_loader(filename);
        }},
        {"-wcub", [&]() {
            CodonUsageBiasCSV cubcsv;
            cubcsv.FASTA_loader(filename);
        }},
        {"-bp", [&]() {
            BasePairCounter bpcounter;
            bpcounter.FASTA_loader(filename);
        }},
        {"-sc", [&]() {
            Seq_colour seqcolour;
            seqcolour.FASTA_loader(filename);
        }},
        {"-mf", [&]() {
            MotifFinder mfinder;
            std::string motif;
            std::cout << "Enter motif: ";
            std::cin >> motif;
            mfinder.FASTAloader(filename, motif);
        }},
        {"-prot", [&](){
            //the old pip3
            int num_threads;
            Proteostats pipeline3;
            
            std::cout << "Select number of threads(default is 2): ";
            std::cin >> num_threads;

            if(num_threads < 2){
                std::cout << "\nWarning!This is a multithreaded function. Cannot be used with less than 2 threads.\n";
                std::cout << "Using 2 threads. . .\n\n\n";
                num_threads = 2;
            }
            pipeline3.FASTA_loader(filename, num_threads);
        }},
        {"-hb", [&](){
            HydrogenBondsCalc hb;
            hb.FASTA_loader(filename);
        }},
        {"-pc", [&](){
            ProteinColor pc;
            pc.FASTA_loader(filename);
        }},
        {"-amb", [&](){
            DNAambiguousStats ambig;
            ambig.FASTA_loader(filename);
        }},
        {"-nw", [&](){
            NeedlemanWunsch nw;
            nw.FASTA_loader(filename);
        }},
        {"-pca", [&](){
            ColorForAminoAcids caa;
            caa.FASTA_loader(filename);
        }},
        {"-asmbl", [&](){
            int numThreads;
            AssemblyStats assstats;
            std::cout << "Select number of threads(default is 2): ";
            std::cin >> numThreads;
            if(numThreads < 2){
                std::cout << "\nWarning!This is a multithreaded function. Cannot be used with less than 2 threads.\n";
                std::cout << "Using 2 threads. . .\n\n\n";
                numThreads = 2;
            }
            assstats.FASTA_loader(filename, numThreads);
        }},
        {"-sw", [&](){
            SmithWaterman swaterman;
            swaterman.FASTA_loader(filename);
        }},
        {"-cx", [&](){
            AssemblyCoverage coverage;
            coverage.FASTA_loader(filename);
        }},
    };
    // Dispatch
    if (dispatch.find(function) != dispatch.end()) {
        dispatch[function]();
    } else {
        message();
        return 1;
    }
    //end clock and print timer results
    auto end_timer = std::chrono::high_resolution_clock::now();
    auto elapsed_timer_s = std::chrono::duration_cast<std::chrono::seconds>(end_timer - start_timer);
    auto elapsed_timer_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_timer - start_timer);

    if (elapsed_timer_s >= std::chrono::seconds(1)){
        auto elapsed_time = elapsed_timer_s.count();
        if (elapsed_time >= 60){
            auto elapsed_time_corrected = elapsed_time / 60;
            std::cout << "Time elapsed: " << elapsed_time_corrected << " m\n";
        }else {
            std::cout << "Time elapsed: " << elapsed_time << " s\n";
        }
    } else {
        std::cout << "Time elapsed: " << elapsed_timer_ms.count() << " ms\n";
    }
    return 0;
}