#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <omp.h>

// DEFAULTS -----------------------
static uint KMER_LEN = 5;
static uint SIGNATURE_WIDTH = 1024;
static uint SIGNATURE_DENSITY = 19;
// --------------------------------
static uint SIGNATURE_BYTE_COUNT = SIGNATURE_WIDTH / 8;

using namespace std;
using namespace chrono;

ofstream OUTFILE;

struct Sequence {
    string id;
    string value;
};

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
            "Usage: ./biosig index [OPTIONS] inputfile1 [inputfile2 [...]] -o outfile.bsig\n"
            "OPTIONS\n"
            "    -kmerlen       Kmer length to hash.\n"
            "                   DEFAULT: 5\n"
            "    -sigwidth      Signature size in bits.\n"
            "                   DEFAULT: 1024\n"
            "    -sigdensity    Signature density (1/x).\n"
            "                   DEFAULT: 19\n";
        case SEARCH:
            return
            "Usage: ./biosig search sigfile.bsig queryFile1 [queryFile2 [...]] -o outfile.tsv";
        default:
            return NULL;
            break;
    }
}

void GenerateSignature(const string &sequence, string &result) {
    default_random_engine random_generator;
    uniform_int_distribution<int> int_distribution(0, SIGNATURE_WIDTH - 1);

    vector<int> signature(SIGNATURE_WIDTH);

    for (uint i = 0; i < sequence.length() - KMER_LEN + 1; i++) {
        seed_seq seed(sequence.begin() + i, sequence.begin() + i + KMER_LEN);
        random_generator.seed(seed);

        uint pos = 0;
        uint set_limit = SIGNATURE_WIDTH / SIGNATURE_DENSITY / 2;
        vector<bool> bits_set(SIGNATURE_WIDTH, false);

        for (uint set = 0; set < set_limit;) {
            pos = int_distribution(random_generator);

            if (!bits_set[pos]) {
                signature[pos] += 1;
                bits_set[pos] = true;
                set++;
            }
        }

        for (uint set = 0; set < set_limit;) {
            pos = int_distribution(random_generator);

            if (!bits_set[pos]) {
                signature[pos] -= 1;
                bits_set[pos] = true;
                set++;
            }
        }
    }

    result.reserve(SIGNATURE_WIDTH);

    for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
        result += (signature[i] > 0 ? '1' : '0');
    }
}

void ForEachSequence(const string filepath, function<void(const Sequence &sequence)> func) {
    ifstream file;
    file.open(filepath);

    Sequence sequence;

    char ch;
    while((ch = file.get()) != EOF) {
        if (ch == '>') {
            if (!sequence.value.empty()) {
                func(sequence);
            }
            getline(file, sequence.id);
            sequence.value.clear();
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            sequence.value += ch;
        }
    }

    if (!sequence.id.empty()) {
        func(sequence);
    }
}

int main(int argc, char * argv[]) {
    if (argc < 2) {
        cerr << Usage(GENERAL) << endl;
        return 1;
    }

    string use = argv[1];

    // INDEX --------------------------------------------------------------------------------------
    if (use == "index") {
        if (argc < 5) {
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

                if (setting == "o") {
                    OUTFILE.open(argv[++i]);
                    if (!OUTFILE.is_open()) {
                        cerr << "Cannot open file for writing: " << argv[i] << endl;
                        cerr << Usage(SEARCH) << endl;
                        return 1;
                    }
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

        if (!OUTFILE.is_open()) {
            cerr << "Please specifiy output file with '-o'" << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        // Initial config metadata
        OUTFILE << KMER_LEN << ',' << SIGNATURE_WIDTH << ',' << SIGNATURE_DENSITY << endl;

        #pragma omp parallel
        {
            #pragma omp single
            {
                for (string filepath : input_files) {
                    cerr << filepath << endl;
                    cerr << "\tIndexing...";
                    high_resolution_clock::time_point t1 = high_resolution_clock::now();

                    ForEachSequence(filepath, [](const Sequence &sequence){
                        #pragma omp task
                        {
                            string signature;
                            GenerateSignature(sequence.value, signature);
                            #pragma omp critical(write_out)
                            {
                                OUTFILE << '>' << sequence.id << endl;
                                OUTFILE << signature << endl;
                            }
                        }
                    });

                    #pragma omp taskwait

                    high_resolution_clock::time_point t2 = high_resolution_clock::now();
                    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
                    cerr << time_span.count() << 's' << endl;
                }
            }
        }

        return 0;
    }
    // --------------------------------------------------------------------------------------------

    // SEARCH -------------------------------------------------------------------------------------
    // if (use == "search") {
    //     if (argc < 6) {
    //         cerr << Usage(SEARCH) << endl;
    //         return 1;
    //     }

    //     string signature_filepath;
    //     vector<string> query_files;

    //     // Argument parsing
    //     for (int i = 2; i < argc; i++) {
    //         string arg = string(argv[i]);

    //         if (arg[0] == '-') {
    //             // Configuration
    //             string setting = arg.substr(1, arg.length());     

    //             if (setting == "o") {
    //                 OUTFILE.open(argv[++i]);
    //                 if (!OUTFILE.is_open()) {
    //                     cerr << "Cannot open file for writing: " << argv[i] << endl;
    //                     cerr << Usage(SEARCH) << endl;
    //                     return 1;
    //                 }
    //                 continue;
    //             }

    //             cerr << "Unknown param: -" << setting << endl;
    //             cerr << Usage(SEARCH) << endl;
    //             return 1;
    //         } else {
    //             if (signature_filepath.empty()) {
    //                 signature_filepath = arg;
    //             } else {
    //                 query_files.push_back(arg);
    //             }
    //         }
    //     }

    //     if (!OUTFILE.is_open()) {
    //         cerr << "Please specifiy output file with '-o'" << endl;
    //         cerr << Usage(SEARCH) << endl;
    //         return 1;
    //     }

    //     ifstream signature_file;
    //     signature_file.open(signature_filepath);

    //     if (!signature_file.is_open()) {
    //         cerr << "Error opening signature file: " << signature_filepath << endl;
    //         cerr << Usage(SEARCH) << endl;
    //         return 1;
    //     }

    //     string line;

    //     // Read and match signature generation metadata
    //     getline(signature_file, line, ',');
    //     KMER_LEN = stoi(line);

    //     getline(signature_file, line, ',');
    //     SIGNATURE_WIDTH = stoi(line);

    //     getline(signature_file, line);
    //     SIGNATURE_DENSITY = stoi(line);

    //     streampos start = signature_file.tellg();

    //     vector<double> results[query_files.size()];
    //     vector<string> signature_targets;

    //     for (uint i = 0; i < query_files.size(); i++) {
    //         string query_signature = GenerateSignature(query_files[i]);

    //         while (getline(signature_file, line)) {
    //             if (line[0] == '>') {
    //                 if (i == 0) signature_targets.push_back(line.substr(1, line.length()));
    //             } else {
    //                 // Signature
    //                 uint hamming_dist = 0;

    //                 for (uint i = 0; i < SIGNATURE_WIDTH; i++) {
    //                     if (query_signature[i] != line[i]) {
    //                         hamming_dist++;
    //                     }
    //                 }

    //                 results[i].push_back(abs((double)hamming_dist - SIGNATURE_WIDTH) / SIGNATURE_WIDTH);
    //                 // results[i].push_back(hamming_dist);
    //             }
    //         }

    //         signature_file.clear();
    //         signature_file.seekg(start);
    //     }

    //     signature_file.close();

    //     // Output result tsv
    //     for (string query_file : query_files) {
    //         OUTFILE << '\t' << query_file;
    //     }

    //     OUTFILE << endl;

    //     for (uint i = 0; i < signature_targets.size(); i++) {
    //         OUTFILE << signature_targets[i];

    //         for (vector<double> result : results) {
    //             OUTFILE << '\t' << to_string(result[i]);
    //         }

    //         OUTFILE << endl;
    //     }

    //     OUTFILE.close();
    //     return 0;
    // }
    // --------------------------------------------------------------------------------------------

    cerr << Usage(GENERAL) << endl;
    return 1;
}