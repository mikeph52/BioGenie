// BioGenie by mikeph_ 2025
// Current version 0.23 
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
// Public Functions 
void title(){
    std::cout << "-----------------------\n";
    std::cout << "BioGenie 0.24.0 pre-release \nby mikeph_ 2025\n\n";
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
    std::cout << "Get the purine/pyrimidine ratio ---> '-pp'.\n";
    std::cout << "Calculate melting temperature (Tm) of DNA sequences using the Wallace Rule(only valid for oligos <20bp) --> '-mt1'.\n";
    std::cout << "Calculate melting temperature (Tm) of DNA sequences using the SantaLucia 1998 nearest-neighbor method --> '-mt2'.\n";
    std::cout << "Get the DNA sequence with colour(EXPERIMENTAL) --> '-sc'.\n";
    std::cout << "Get the complement DNA sequence with colour(EXPERIMENTAL) --> '-cc'.\n";
    std::cout << "Get the Open Reading Frame(ORF) ---> '-orf'.\n";
    std::cout << "Generate cDNA sequence FASTA ---> '-cw'.\n";
    std::cout << "Generate Reverse cDNA sequence FASTA ---> '-rcw'.\n";
    std::cout << "Generate mRNA sequence FASTA ---> '-tw'.\n";
    std::cout << "Calculate Codon Usage Bias(CUB) ---> '-cub'.\n";
    std::cout << "Export Codon Usage Bias(CUB) to CSV file ---> '-wcub'.\n";
    std::cout << "Calculate the Number of Base Pairs(bp) ---> '-bp'.\n";
    std::cout << "Preset pipeline 1 ---> '-pip1'. Returns the codon number and GC%.\n";
    std::cout << "Preset pipeline 2 ---> '-pip2'. Ideal for Primer design.\n\n";
    std::cout << "For more info visit the github page: https://github.com/mikeph52/BioGenie/blob/main/documentation.md\n";
    std::cout << "More functions will be added in the future.\n\n";
    std::cout << "-----------------------------------------------------------\n";
}

void message(){
        std::cerr << "Usage: biogenie <function> <FASTA_file_path>\n\n";
        std::cerr << "[-c complement DNA sequence][-rc reverse complement DNA sequence]\n";
        std::cerr << "[-nc codon number][-t mRNA][-gc GC percentage calculator][-p protein chain]\n";
        std::cerr << "[-ss FASTA sequencies separator][-sh FASTA sequencies headers][-tr DNA Trimmer][-bp Base Pairs]\n";
        std::cerr << "[-pp purine/pyrimidine ratio][-mt1 melting temp.(Wallace rule)][-mt2 melting temp.(Nearest-neighbour)]\n";
        std::cerr << "[-cc cDNA coloured][-orf ORF Finder][-cw generate cDNA fasta][-rcw Reverse cDNA fasta][-tw mRNA fasta]\n";
        std::cerr << "[-cub Codon Usage Bias][-wcub Codon Usage Bias to CSV][-sc colour sequence]\n";
        std::cerr << "[-pip1 Preset pipeline 1][-pip2 Preset pipeline 2]\n";
        std::cerr << "[Use '-help me' for documentation.]\n\n\n ";
        std::cerr << "For more info visit the github page:\nhttps://github.com/mikeph52/BioGenie\n\n";
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

                if (codon == "ATG") { // start codon found
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
                    // move to next orf after stop codon
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
            default:  return 'N'; // Unknown base
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

        // Build AA -> Codons map
        for (auto &kv : codonTable) {
            if (kv.second != "*") // skip stops
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

        // Build AA -> Codons map (skip stops)
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

        // Save CUB as CSV
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
    Classes:PurinePyrimidineRatioAnalyzer, MeltingTempCalculator2, GCCalc, BasePairCounter.*/
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

            //std::cout << ">" << header << "\n";   <--- fix
            std::cout << "A: " << a << ", T: " << t << ", G: " << g << ", C: " << c << "\n";
            std::cout << "Melting Temperature (NN model): " << tm << " °C\n";
            // std::cout << "-----------------------------------\n";
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
                    int bpCount = CountBases(sequence);
                    std::cout << "Base pairs: " << bpCount << "\n";
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
            int bpCount = CountBases(sequence);
            std::cout << "Base pairs: " << bpCount << "\n";
            TmCalc(header, sequence);
            double GCContent = GCCon(sequence);
            std::cout << "GC Content = " << std::fixed << std::setprecision(2) << GCContent << "%\n";
            std::cout << "-----------------------------------\n";
        }

        std::cout << "Process completed.\n";
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
        Pipeline1 pipeline1;
        pipeline1.FASTA_loader(filename);
        
    }else if(function == "-pp"){
        PurinePyrimidineRatioAnalyzer ppanalyzer;
        ppanalyzer.FASTA_loader(filename);

    }else if(function == "-mt1"){
        MeltingTempCalculator1 mtcalc1;
        mtcalc1.FASTA_loader(filename);

    }else if(function == "-cc"){
        cDNA_colour dnacolour;
        dnacolour.FASTA_loader(filename);

    }    else if(function == "-pip2"){
        Pipeline2 pipeline2;
        pipeline2.FASTA_loader(filename);
    }
    else if(function == "-mt2"){
        MeltingTempCalculator2 mtcalc2;
        mtcalc2.FASTA_loader(filename);

    }else if(function == "-orf"){
        ORFFinder orffinder;
        orffinder.FASTA_loader(filename);

    }else if (function == "-cw" && !filename.empty()) {
        std::string outputFile;
        std::cout << "Enter output filename: ";
        std::cin >> outputFile;
        DNAcompToFile writecomplimentary;
        writecomplimentary.FASTA_writer(filename, outputFile);

    }else if(function == "-rcw" && !filename.empty()){
        std::string outputFile;
        std::cout << "Enter output filename: ";
        std::cin >> outputFile;
        ReverseComplementDNAToFile writereverse;
        writereverse.FASTA_writer(filename, outputFile);

    }else if(function == "-tw" && !filename.empty()){
        std::string outputFile;
        std::cout << "Enter output filename: ";
        std::cin >> outputFile;
        TranscriptionToFile writeRNA;
        writeRNA.FASTA_writer(filename, outputFile);

    }else if(function == "-cub"){
        CodonUsageBias cub;
        cub.FASTA_loader(filename);

    }else if(function == "-wcub"){
        CodonUsageBiasCSV cubcsv;
        cubcsv.FASTA_loader(filename);

    }else if(function == "-bp"){
        BasePairCounter bpcounter;
        bpcounter.FASTA_loader(filename);

    }else if(function == "-sc"){
        Seq_colour seqcolour;
        seqcolour.FASTA_loader(filename);

    } else {
        message();
        return 1;
    }
    return 0;
}