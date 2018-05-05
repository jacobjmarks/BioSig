#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <omp.h>
#include <algorithm>

// DEFAULTS -----------------------
static uint KMER_LEN = 5;
static uint SIGNATURE_WIDTH = 1024;
static uint SIGNATURE_DENSITY = 19;
static double DIST_THRESHOLD = 0.5;
// --------------------------------

using namespace std;
using namespace chrono;

ofstream OUTFILE;

typedef struct {
    string id;
    string value;
} Sequence, Signature;

struct SearchResult {
    string target;
    uint hamming_dist;
    double normalised_dist;

    SearchResult(string target, uint hamming_dist, double normalised_dist) {
        SearchResult::target = target;
        SearchResult::hamming_dist = hamming_dist;
        SearchResult::normalised_dist = normalised_dist;
    }
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

    result.reserve(SIGNATURE_WIDTH / 8);

    char bit = 0;
    for (uint i = 0; i < SIGNATURE_WIDTH / 8; i++) {
        for (uint j = 0; j < 8; j++) {
            bit ^= (-(signature[(i * 8) + j] > 0 ? 1 : 0) ^ bit) & (1UL << (-j + 7));
        }
        result += bit;
    }
}

template<typename FUNC>
void ForEachSequence(const string filepath, const FUNC &func) {
    ifstream file;
    file.open(filepath);

    #pragma omp critical(exit)
    if (!file.is_open()) {
        cerr << endl;
        cerr << "Error opening file: " << filepath << endl;
        exit(1);
    }

    Sequence sequence;

    char ch;
    while((ch = file.get()) != EOF) {
        if (ch == '>') {
            if (!sequence.value.empty()) {
                func(sequence);
            }
            sequence.id.clear();
            sequence.value.clear();
            while((ch = file.get()) != '\n' && ch != '\r') sequence.id += ch;
        } else if (ch == '\n' || ch == '\r') {
            // Skip
        } else {
            sequence.value += ch;
        }
    }

    if (!sequence.id.empty()) {
        func(sequence);
    }

    file.close();
}

template<typename FUNC>
void ForEachSignature(const string filepath, const FUNC &func) {
    ifstream file;
    file.open(filepath);

    #pragma omp critical(exit)
    if (!file.is_open()) {
        cerr << endl;
        cerr << "Error opening file: " << filepath << endl;
        exit(1);
    }

    string line;

    // Skip meta
    getline(file, line);

    Signature signature;

    char ch;
    while((ch = file.get()) != EOF) {
        if (ch == '>') {
            getline(file, signature.id);
            signature.value.clear();
            for (uint i = 0; i < SIGNATURE_WIDTH / 8; i++) {
                signature.value += (ch = file.get());
            }
            func(signature);
        }
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
        #pragma omp single
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

        return 0;
    }
    // --------------------------------------------------------------------------------------------

    // SEARCH -------------------------------------------------------------------------------------
    if (use == "search") {
        if (argc < 6) {
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        string target_filepath;
        vector<string> query_files;

        // Argument parsing
        for (int i = 2; i < argc; i++) {
            string arg = string(argv[i]);

            if (arg[0] == '-') {
                // Configuration
                string setting = arg.substr(1, arg.length());     

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
                cerr << Usage(SEARCH) << endl;
                return 1;
            } else {
                if (target_filepath.empty()) {
                    target_filepath = arg;
                } else {
                    query_files.push_back(arg);
                }
            }
        }

        if (!OUTFILE.is_open()) {
            cerr << "Please specifiy output file with '-o'" << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        ifstream target_file;
        target_file.open(target_filepath);

        if (!target_file.is_open()) {
            cerr << "Error opening signature file: " << target_filepath << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        string line;

        // Read and match signature generation metadata
        getline(target_file, line, ',');
        KMER_LEN = stoi(line);

        getline(target_file, line, ',');
        SIGNATURE_WIDTH = stoi(line);

        getline(target_file, line);
        SIGNATURE_DENSITY = stoi(line);

        target_file.close();

        cerr << "Searching...";
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        #pragma omp parallel
        #pragma omp single
        for (string query_filepath : query_files) {
            ForEachSignature(query_filepath, [&](const Signature &query_signature) {
                #pragma omp task
                {
                    vector<SearchResult> these_results;

                    ForEachSignature(target_filepath, [&](const Signature &target_signature) {
                        uint hamming_dist = 0;

                        for (uint i = 0; i < SIGNATURE_WIDTH / 8; i++) {
                            char xOR = query_signature.value[i] ^ target_signature.value[i];

                            // Brian Kernighan pop count
                            while (xOR != 0) {
                                xOR &= (xOR - 1);
                                hamming_dist++;
                            }
                        }

                        double normalised_dist =
                            abs((double)hamming_dist - SIGNATURE_WIDTH) / SIGNATURE_WIDTH;

                        if (normalised_dist >= DIST_THRESHOLD) {
                            these_results.push_back(
                                SearchResult(target_signature.id, hamming_dist, normalised_dist)
                            );
                        }
                    });

                    sort(these_results.begin(), these_results.end(),
                        [](const SearchResult &a, const SearchResult &b) {
                            return a.hamming_dist < b.hamming_dist;
                        }
                    );

                    #pragma omp critical(write_out)
                    for (const SearchResult &result : these_results) {
                        OUTFILE
                            << query_signature.id
                            << '\t'
                            << result.target
                            << '\t'
                            << result.hamming_dist
                            << '\t'
                            << result.normalised_dist
                            << endl;
                    }
                }
            });
        }

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
        cerr << time_span.count() << 's' << endl;

        return 0;
    }
    // --------------------------------------------------------------------------------------------

    cerr << Usage(GENERAL) << endl;
    return 1;
}