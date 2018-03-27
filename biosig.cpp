#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>

// DEFAULTS -----------------------
static uint KMER_LEN = 5;
static uint SIGNATURE_WIDTH = 1024;
static uint SIGNATURE_DENSITY = 19;
// --------------------------------

using namespace std;

vector<int> hashKmer(string kmer) {
    srand(hash<string>{}(kmer));

    vector<int> kmerHash(SIGNATURE_WIDTH);

    uint pos = 0;
    uint set_limit = SIGNATURE_WIDTH / SIGNATURE_DENSITY / 2;

    // Fill half with 1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_WIDTH;

        if (!kmerHash[pos]) {
            kmerHash[pos] = 1;
            set++;
        }
    }

    // Fill other half with -1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_WIDTH;

        if (!kmerHash[pos]) {
            kmerHash[pos] = -1;
            set++;
        }
    }

    return kmerHash;
}

void generateSignature(string filename) {
    ifstream file;
    file.open(filename);

    if(!file.is_open()) {
        cout << "Error opening file: " << filename << endl;
        return;
    }

    vector<int> signature(SIGNATURE_WIDTH);

    string kmerBuffer[KMER_LEN];
    uint maxBufferIndex = 1;

    char ch;
    while((ch = file.get()) != EOF) {
        if (ch == '>') {
            while ((ch = file.get()) != '\n');
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            for (uint i = 0; i < maxBufferIndex; i++) {
                kmerBuffer[i] += ch;

                if (kmerBuffer[i].length() == KMER_LEN) {
                    vector<int> kmerHash = hashKmer(kmerBuffer[i]);

                    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
                        signature[i] += kmerHash[i];
                    }
                    
                    kmerBuffer[i].clear();
                }
            }

            if (maxBufferIndex < KMER_LEN) maxBufferIndex++;
        }
    }

    // Output flattened signature
    uint popcount = 0;
    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
        uint bit = (signature[i] > 0 ? 1 : 0);
        popcount += bit;
        cout << bit;
    }
    cout << endl;
}

void readDirectory(string directory) {

}

int main(int argc, char * argv[]) {
    if (argc < 2) {
        cout << "Please specify filename." << endl;
        return 1;
    }

    vector<string> inputFiles;

    // Argument parsing
    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);

        if (arg[0] == '-') {
            // Configuration
            string setting = arg.substr(1, arg.length());

            if (setting == "kmerlen") {
                KMER_LEN = stoi(argv[++i]);
                continue;
            }

            if (setting == "sigwidth") {
                SIGNATURE_WIDTH = stoi(argv[++i]);
                cerr << SIGNATURE_WIDTH << endl;
                continue;
            }

            if (setting == "sigdensity") {
                SIGNATURE_DENSITY = stoi(argv[++i]);
                continue;
            }
        } else {
            // Input file
            inputFiles.push_back(arg);
        }
    }

    for (string file : inputFiles) {
        cerr << file << endl;
        cerr << "\tIndexing...";
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        generateSignature(file);

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        cerr << time_span.count() << 's' << endl;
    }

    return 0;
}