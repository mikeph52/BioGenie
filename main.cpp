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
    std::cout << "BioGenie 0.13.0\nby mikeph_ 2025\n\n";
    //std::cout << "-----------------------------------\n\n";
    
}

void helpme(){
    std::cout << "-----------------------DOCUMENTATION-----------------------\n";
    std::cout << "BioGenie uses functions to execute different tools for different applications.\n\n";
    std::cout << "Get the complement DNA sequence --> '-c'.\n";
    std::cout << "Get the reverse complement DNA sequence --> '-rc'.\n";
    std::cout << "Get the codon number --> '-nc'.\n";
    std::cout << "Get the mRNA --> '-t'.\n";
    std::cout << "GC percentage calculation --> '-gc'.\n";
    std::cout << "Generate the aminoacids(Protein chain) ---> '-p'.\n";
    std::cout << "Separate different sequencies in a FASTA file ---> '-ss'\n";
    std::cout << "Print the different sequence headers from a FASTA file ---> '-sh'\n";
    std::cout << "Trim DNA ---> 'tr'. It uses 0-based indexing (start = 0 is the first base).\n";
    std::cout << "Get the purine/pyrimidine ratio --> '-pp'.\n";
    std::cout << "Calculate melting temperature (Tm) of DNA sequences --> '-mt'.\n";
    std::cout << "Preset pipeline 1 ---> -pip1. Returns the codon number and GC%.\n";
    std::cout << "Preset pipeline 2 ---> -pip2. Returns the purine/pyrimidine ratio, GC% and Melting temperature.\nIdeal for Primer design.\n\n";
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
                std::cout << ">" << header << "Protein:\n" << complement << "\n";
            }

            std::cout << "-----------------------------------\n\n\n";
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
        std::cout << "\n--------- Sequence Headers ---------\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                std::cout << line << "\n";
            }
        }

        std::cout << "-----------------------------------\n\n\n";
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

class MeltingTempCalculator {
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

class Pipeline1 {
    /*This is a pipeline that contains the following classes and functions:
    Classes:GCCalc, CodonNumber.*/
    private:
    // GC content function
        double GCContent1(const std::string& sequence) const {
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
    // codon count function
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
                exit(1);
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
                        std::cout << ">" << header << "\nCodon count:" << codons << "\n";
                        double gcContent = GCContent1(sequence);
                        std::cout << "GC Content = " << std::fixed << std::setprecision(2) << gcContent << "%\n";
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
                std::cout << ">" << header << "\nCodon count:" << codons << "\n";
                double gc = GCContent1(sequence);
                std::cout << "GC Content = " << std::fixed << std::setprecision(2) << gc << "%\n";
            }

            std::cout << "-----------------------------------\n\n\n";
            std::cout << "Process completed.\n";

            fastaFile.close();
        }
};

class Pipeline2 {
     /*This is a pipeline that contains the following classes and functions:
    Classes:PurinePyrimidineRatioAnalyzer, MeltingTempCalculator, GCCalc.*/
    private:
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
                    break; // Skip 'N' or unknown bases
            }
        }

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
        //std::cout << "\n-----------------------------------\n";
    }
    void TmCalc(const std::string& header, const std::string& sequence) const {
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

        //std::cout << ">" << header << "\n";
        std::cout << "A: " << a << ", T: " << t << ", G: " << g << ", C: " << c << "\n";
        std::cout << "Melting Temperature (Tm): " << tm << "°C\n";
    }
    double GCCon(const std::string& sequence)const {
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
    void FASTA_loader(const std::string& filename) const {
        std::ifstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            std::cerr << "Error: Unable to open file: " << filename << "\n";
            return;
        }

        std::string line, header, sequence;
        std::cout << "\n-----"<< filename<<"------\n";
        //std::cout << "\n-----------------------\n";

        while (std::getline(fastaFile, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!sequence.empty()) {
                    ppRatio(header, sequence);
                    TmCalc(header, sequence);
                    double GCContent = GCCon(sequence);
                    std::cout << "GC Content = " << std::fixed << std::setprecision(2) << GCContent << "%\n";
                    std::cout << "-----------------------------------\n";
                    sequence.clear();
                }
                header = line.substr(1);
            } else {
                sequence += line;
            }
        }

        if (!sequence.empty()) {
            ppRatio(header, sequence);
            TmCalc(header, sequence);
            double GCContent = GCCon(sequence);
            std::cout << "GC Content = " << std::fixed << std::setprecision(2) << GCContent << "%\n";
            std::cout << "-----------------------------------\n";
        }

        std::cout << "Process completed.\n";
        fastaFile.close();
    }

};

// Main Function 
int main(int argc, char* argv[]){
    if (argc != 3){
        title();
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence]\n";
        std::cerr << "[-nc codon number][-t mRNA][-gc GC percentage calculator][-p protein chain]\n";
        std::cerr << "[-ss FASTA sequencies separator][-sh FASTA sequencies headers][-tr DNA Trimmer]\n";
        std::cerr << "[-pp purine/pyrimidine ratio][-mt melting temp. calculator][-pip1 Preset pipeline 1]\n";
        std::cerr << "[-pip2 Preset pipeline 2]";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n ";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
        return 1;
    }
    
    title();
    std::string filename = argv[2];
    std::string function = argv[1];
    
    //FASTA verifier
    FastaVerifier verifier(filename);
    if (verifier.verify()) {
        std::cout << "FASTA file status  [OK]\n";
    } else {
        std::cerr << "FASTA file status  [FAULT]\n";
    }
    // Main if body
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

    }else if(function == "-ss"){
        FASTAChromosomeSeparator splitter;
        splitter.FASTA_loader(filename);

    }else if(function == "-sh"){
        FASTASequenceHeader headers;
        headers.FASTA_loader(filename);

    }else if(function == "-tr"){
        DNATrimmer trim;
        int start_position, end_position;
        std::cout << "Enter the starting position:";
        std::cin >> start_position;
        std::cout << "Enter the end position:";
        std::cin >> end_position;
        trim.FASTA_loader(filename, start_position, end_position);

    }else if(function == "-pip1"){
        // THIS IS A TESTING FUNCTION
        Pipeline1 pipeline1;
        pipeline1.FASTA_loader(filename);
        

    }else if(function == "-pp"){
        PurinePyrimidineRatioAnalyzer ppanalyzer;
        ppanalyzer.FASTA_loader(filename);

    }else if(function == "-mt"){
        MeltingTempCalculator mtcalc;
        mtcalc.FASTA_loader(filename);

    }else if(function == "-pip2"){
        Pipeline2 pipeline2;
        pipeline2.FASTA_loader(filename);
    }
    else {
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence]\n";
        std::cerr << "[-nc codon number][-t mRNA][-gc GC percentage calculator][-p protein chain]\n";
        std::cerr << "[-ss FASTA sequencies separator][-sh FASTA sequencies headers][-tr DNA Trimmer]\n";
        std::cerr << "[-pp purine/pyrimidine ratio][-mt melting temp. calculator][-pip1 Preset pipeline 1]\n";
        std::cerr << "[-pip2 Preset pipeline 2]";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n ";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
        return 1;
    }


    return 0;
}