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
using namespace chrono;

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
            return
            "Usage: ./biosig search sigfile.bsig queryFile1 [queryFile2 [...]] > results.txt";
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

string GenerateSignature(string filename) {
    ifstream file;
    file.open(filename);

    if(!file.is_open()) {
        cerr << endl;
        __throw_runtime_error(("Error opening file: " + filename).c_str());
    }

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

    file.close();

    // Flatten and return
    string flatsig;
    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
        flatsig += (signature[i] > 0 ? '1' : '0');
    }

    return flatsig;
}

void ReadDirectory(string directory) {

}

int main(int argc, char * argv[]) {
    if (argc < 2) {
        cerr << Usage(GENERAL) << endl;
        return 1;
    }

    string use = argv[1];

    // INDEX --------------------------------------------------------------------------------------
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
            high_resolution_clock::time_point t1 = high_resolution_clock::now();

            cout << '>' << file.substr(file.find_last_of('/') + 1, file.length()) << endl;

            string signature = GenerateSignature(file);

            for (char bit : signature) {
                cout << bit;
            }
            cout << endl;

            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
            cerr << time_span.count() << 's' << endl;
        }

        return 0;
    }
    // --------------------------------------------------------------------------------------------

    // SEARCH -------------------------------------------------------------------------------------
    if (use == "search") {
        if (argc < 4) {
            cerr << "Please specify target signature file and query sequence file(s)." << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        string signature_filepath;
        vector<string> query_files;

        // Argument parsing
        for (int i = 2; i < argc; i++) {
            string arg = string(argv[i]);

            if (arg[0] == '-') {
                // Configuration
                string setting = arg.substr(1, arg.length());                
                cerr << "Unknown param: -" << setting << endl;
                cerr << Usage(SEARCH) << endl;
                return 1;
            } else {
                if (signature_filepath.empty()) {
                    signature_filepath = arg;
                } else {
                    query_files.push_back(arg);
                }
            }
        }

        ifstream signature_file;
        signature_file.open(signature_filepath);

        if (!signature_file.is_open()) {
            cerr << "Error opening signature file: " << signature_filepath << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        string line;

        // Read and match signature generation metadata
        getline(signature_file, line, ',');
        KMER_LEN = stoi(line);

        getline(signature_file, line, ',');
        SIGNATURE_WIDTH = stoi(line);

        getline(signature_file, line);
        SIGNATURE_DENSITY = stoi(line);

        streampos start = signature_file.tellg();

        for (string query_file : query_files) {
            string query_signature = GenerateSignature(query_file);

            cout << query_file << endl;

            while (getline(signature_file, line)) {
                if (line[0] == '>') {
                    // Metadata
                    cout << '\t' << line.substr(1, line.length()) << '\t';
                } else {
                    // Signature
                    uint hamming_dist = 0;

                    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
                        if (query_signature[i] != line[i]) {
                            hamming_dist++;
                        }
                    }

                    cout << abs((double)hamming_dist - SIGNATURE_WIDTH) / SIGNATURE_WIDTH << endl;
                }
            }

            signature_file.clear();
            signature_file.seekg(start);
        }

        signature_file.close();
        return 0;
    }
    // --------------------------------------------------------------------------------------------

    cerr << Usage(GENERAL) << endl;
    return 1;
}