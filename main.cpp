#include <iostream>
#include <string>
#include <cctype>
#include <iomanip> 
#include <fstream>
#include <cstdlib>

void title(){
    std::cout << "-----------------------\n";
    std::cout << "BioGenie 0.1 for macOS\nby mikeph_ 2025\n";
    //std::cout << "-----------------------------------\n\n";
   
}

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
    

int main(int argc, char* argv[]){
    if (argc != 2){
        std::cerr << "Usage: gcgenie <FASTA_file_path>\n";
        return 1;
    }
    
    title();
    std::string filename = argv[1];
    GCCalc GCcalculator;
    GCcalculator.FASTA_loader(filename);
    return 0;
}