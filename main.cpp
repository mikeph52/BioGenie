// BioGenie by mikeph_ 2025
#include <iostream>
#include <string>
#include <cctype>
#include <iomanip> 
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>

// Public Functions 
void title(){
    std::cout << "-----------------------\n";
    std::cout << "BioGenie 0.8.0 \nby mikeph_ 2025\n\n";
    //std::cout << "-----------------------------------\n\n";
   
}

void helpme(){
    std::cout << "-----------------------DOCUMENTATION-----------------------\n";
    std::cout << "BioGenie uses functions to execute different tools for different applications.\n";
    std::cout << "Get the complement DNA sequence --> '-c'.\n";
    std::cout << "Get the reverse complement DNA sequence --> '-rc'.\n";
    std::cout << "Get the codon number --> '-nc'.\n";
    std::cout << "Get the mRNA --> '-t'.\n";
    std::cout << "GC percentage calculation --> '-gc'.\n";
    std::cout << "Generate the aminoacids(Protein chain) ---> '-p'.\n";
    std::cout << "More functions will be added in the future.\n";
    std::cout << "-----------------------------------------------------------\n";
}

//Genetic code
const std::unordered_map<std::string, char> codonTable = {
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

// Arg Classes
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

            std::string line;
            std::string sequence;
            std::string Header;

            std::cout << "\n-----------------------------------\n";

            while (std::getline(fastaFile, line)) {
                if (line.empty()) continue;

                if (line[0] == '>') {
                    if (!sequence.empty()) {
                        double gcContent = GCContent(sequence);
                        std::cout << Header << ":\nGC Content = " << std::fixed << std::setprecision(2) << gcContent << "%\n";
                        std::cout << "\n-----------------------------------\n";
                        sequence.clear();
                    }
                    Header = line.substr(1);
                } else {
                    sequence += line;
                }
            }

            if (!sequence.empty()) {
                double gc = GCContent(sequence);
                std::cout << Header << ":\nGC Content = " << std::fixed << std::setprecision(2) << gc << "%\n";
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
                std::cout << ">" << header << "RNA:\n" << complement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";
            fastaFile.close();
        }
};

// Main Function 
int main(int argc, char* argv[]){
    if (argc != 3){
        std::cout << "-----------------------\n";
        std::cout << "BioGenie 0.8.0 \nby mikeph_ 2025\n\n";
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence]\n";
        std::cerr << "[-nc codon number][-t mRNA][-gc GC percentage calculator][-p protein chain]\n";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n ";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
        return 1;
    }
    
    title();
    std::string filename = argv[2];
    std::string function = argv[1];

    if (function == "-gc"){
        GCCalc GCcalculator;
        GCcalculator.FASTA_loader(filename);
    } else if (function == "-nc"){
        CodonNumber codoncounter;
        codoncounter.FASTA_loader(filename);
    } else if (function == "-c"){
        DNAcomplementary DNAcomp;
        DNAcomp.FASTA_loader(filename);
    } else if (function == "-rc"){
        ReverseComplementDNA revDNA;
        revDNA.FASTA_loader(filename);
    } else if (function == "-t"){
        Transcription transciptedRNA;
        transciptedRNA.FASTA_loader(filename);
    } else if(function == "-help"){
        helpme();
    } else if(function == "-p"){
        ProteinChain protein;
        protein.FASTA_loader(filename);

    }else {
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence]\n";
        std::cerr << "[-nc codon number][-t mRNA][-gc GC percentage calculator][-p protein chain]\n";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n ";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
        return 1;
    }


    return 0;
}