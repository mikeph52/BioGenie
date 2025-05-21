#include <iostream>
#include <string>
#include <cctype>
#include <iomanip> 
#include <fstream>
#include <cstdlib>
#include <algorithm>

// Public Functions 
void title(){
    std::cout << "-----------------------\n";
    std::cout << "BioGenie 0.4.0 for macOS\nby mikeph_ 2025\n";
    //std::cout << "-----------------------------------\n\n";
   
}

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
    
class DNAcomplimentary{
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

// Main Function 
int main(int argc, char* argv[]){
    if (argc != 2){
        std::cout << "-----------------------\n";
        std::cout << "BioGenie 0.4.0 for macOS\nby mikeph_ 2025\n";
        std::cerr << "Usage: gcgenie <FASTA_file_path>\n";
        return 1;
    }
    
    title();
    std::string filename = argv[1];

   /* GCCalc GCcalculator;
    GCcalculator.FASTA_loader(filename);

    DNAcomplimentary DNAcomp;
    DNAcomp.FASTA_loader(filename);

    Transcription transciptedRNA;
    transciptedRNA.FASTA_loader(filename); 

    ReverseComplementDNA revDNA;
    revDNA.FASTA_loader(filename);*/

    CodonNumber codoncounter;
    codoncounter.FASTA_loader(filename);



    return 0;
}