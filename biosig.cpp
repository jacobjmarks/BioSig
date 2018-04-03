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

enum Use {
    GENERAL, INDEX, SEARCH
};

string Usage(Use use) {
    switch(use) {
        case GENERAL:
            return
            "Usage: ./biosig [OPTION]\n"
            "OPTION\n"
            "    index     Generates a signature file from the given sequence files.\n"
            "    search    Searches the documents in a given signature file.\n";
        case INDEX:
            return
            "Usage: ./biosig index [OPTIONS] inputfile1 [inputfile2 [...]] > outputfile.bsig\n"
            "OPTIONS\n"
            "    -kmerlen       Kmer length to hash.\n"
            "                   DEFAULT: 5\n"
            "    -sigwidth      Signature size in bits.\n"
            "                   DEFAULT: 1024\n"
            "    -sigdensity    Signature density (1/x).\n"
            "                   DEFAULT: 19\n";
        case SEARCH:
        default:
            return NULL;
            break;
    }
}

vector<int> HashKmer(string kmer) {
    srand(hash<string>{}(kmer));

    vector<int> kmer_hash(SIGNATURE_WIDTH);

    uint pos = 0;
    uint set_limit = SIGNATURE_WIDTH / SIGNATURE_DENSITY / 2;

    // Fill half with 1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_WIDTH;

        if (!kmer_hash[pos]) {
            kmer_hash[pos] = 1;
            set++;
        }
    }

    // Fill other half with -1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_WIDTH;

        if (!kmer_hash[pos]) {
            kmer_hash[pos] = -1;
            set++;
        }
    }

    return kmer_hash;
}

void GenerateSignature(string filename) {
    ifstream file;
    file.open(filename);

    if(!file.is_open()) {
        cerr << "Error opening file" << endl;
        return;
    }

    cout << '>' << filename.substr(filename.find_last_of('/') + 1, filename.length()) << endl;

    vector<int> signature(SIGNATURE_WIDTH);

    string kmer_buffer[KMER_LEN];
    uint max_buffer_index = 1;

    char ch;
    while((ch = file.get()) != EOF) {
        if (ch == '>') {
            while ((ch = file.get()) != '\n');
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            for (uint i = 0; i < max_buffer_index; i++) {
                kmer_buffer[i] += ch;

                if (kmer_buffer[i].length() == KMER_LEN) {
                    vector<int> kmer_hash = HashKmer(kmer_buffer[i]);

                    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
                        signature[i] += kmer_hash[i];
                    }
                    
                    kmer_buffer[i].clear();
                }
            }

            if (max_buffer_index < KMER_LEN) max_buffer_index++;
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

void ReadDirectory(string directory) {

}

int main(int argc, char * argv[]) {
    if (argc < 2) {
        cerr << Usage(GENERAL) << endl;
        return 1;
    }

    string use = argv[1];

    if (use == "index") {
        if (argc < 3) {
            cerr << Usage(INDEX) << endl;
            return 1;
        }

        vector<string> input_files;

        // Argument parsing
        for (int i = 2; i < argc; i++) {
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
                
                cerr << "Unknown param: -" << setting << endl;
                cerr << Usage(INDEX) << endl;
                return 1;
            } else {
                // Input file
                input_files.push_back(arg);
            }
        }

        if (!input_files.size()) {
            cerr << "Please specify input file(s)." << endl;
            cerr << Usage(INDEX) << endl;
            return 1;
        }

        // Initial config metadata
        cout << KMER_LEN << ',' << SIGNATURE_WIDTH << ',' << SIGNATURE_DENSITY << endl;

        for (string file : input_files) {
            cerr << file << endl;
            cerr << "\tIndexing...";
            chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

            GenerateSignature(file);

            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            cerr << time_span.count() << 's' << endl;
        }

        return 0;
    }

    if (use == "search") {
        cerr << "NOT IMPLEMENTED" << endl;
        return 0;
    }

    cerr << Usage(GENERAL) << endl;
    return 1;
}