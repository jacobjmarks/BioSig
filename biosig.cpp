#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <omp.h>
#include <algorithm>
#include <queue>
#include <regex>

// DEFAULTS -----------------------
static uint KMER_LEN = 5;
static uint SIGNATURE_WIDTH = 1024;
static uint SIGNATURE_DENSITY = 19;
static double DIST_THRESHOLD = 0.0;
static uint RESULT_LIMIT = 0;
// --------------------------------

using namespace std;
using namespace chrono;

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
            "    index     Generates a signature file from the given sequence file/s.\n"
            "    search    Searches the documents in a given signature file.\n";
        case INDEX:
            return
            "Usage: ./biosig index [OPTIONS] sequenceFile [sequenceFile2 [...]] -o outfile.bsig\n"
            "OPTIONS\n"
            "    -kmerlen       Kmer length to hash.\n"
            "                     DEFAULT: 5\n"
            "    -sigwidth      Signature size in bits.\n"
            "                     DEFAULT: 1024\n"
            "    -sigdensity    Signature density (1/x).\n"
            "                     DEFAULT: 19\n";
        case SEARCH:
            return
            "Usage: ./biosig search [OPTIONS] targetSig.bsig querySig.bsig [querySig2.bsig [...]] -o resultfile\n"
            "OPTIONS\n"
            "    -threshold    Filter results lower than the given threshold (0-1).\n"
            "                    DEFAULT: 0.0\n"
            "    -top          Retain only the top k results.\n"
            "                    DEFAULT: 0 (disabled)\n"
            "    -format       Result output format. (tsv | csv | trec)\n"
            "                    DEFAULT: tsv\n";
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
void ForEachSignature(const string sigfilepath, const FUNC &func) {
    ifstream sigfile;
    sigfile.open(sigfilepath, ifstream::binary);

    ifstream headfile;
    headfile.open(sigfilepath + ".head");

    if (!sigfile.is_open()) {
        cerr << endl;
        cerr << "Error opening signature file: " << sigfilepath << endl;
        #pragma omp critical(exit)
        exit(1);
    }

    if (!headfile.is_open()) {
        cerr << endl;
        cerr << "Error opening header file: " << sigfilepath + ".head" << endl;
        #pragma omp critical(exit)
        exit(1);
    }

    string line;

    // Skip meta
    getline(headfile, line);

    Signature signature;

    char ch;
    while (getline(headfile, line)) {
        signature.id = line;
        for (uint i = 0; i < SIGNATURE_WIDTH / 8; i++) {            
            signature.value += (ch = sigfile.get());
        }
        func(signature);
        signature.value.clear();
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

        static vector<string> input_files;
        static ofstream headfile;
        static ofstream sigfile;

        static bool use_regex = false;
        static regex seqid_regex;

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

                if (setting == "match") {
                    use_regex = true;
                    seqid_regex.assign(argv[++i]);
                    continue;
                }

                if (setting == "o") {
                    sigfile.open(argv[++i]);
                    headfile.open(argv[i] + (string)".head");

                    if (!sigfile.is_open()) {
                        cerr << "Cannot open file for writing: " << argv[i] << endl;
                        return 1;
                    }

                    if (!headfile.is_open()) {
                        cerr << "Cannot open file for writing: " << argv[i] << ".head" << endl;
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

        if (!headfile.is_open() || !sigfile.is_open()) {
            cerr << "Please specifiy output file with '-o'" << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        // Prepass; Get sequence count
        cerr << "Prepass...";

        static unsigned long long int sequence_count = 0;

        for (string filepath : input_files) {
            ForEachSequence(filepath, [](Sequence &sequence) {
                sequence_count++;
            });
        }

        cerr << "DONE" << endl;

        // Write initial metadata
        headfile << KMER_LEN << ',' << SIGNATURE_WIDTH << ',' << SIGNATURE_DENSITY << ',' << sequence_count << endl;

        cerr << "Indexing " << input_files.size() << " files with " << sequence_count << " sequences..." << endl;

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        static unsigned long long int sequence_index = 0;
        static uint percent_complete = 0;

        cerr << "0%";

        #pragma omp parallel
        #pragma omp single
        for (string filepath : input_files) {
            ForEachSequence(filepath, [](Sequence &sequence){
                // Delegate signature generation
                #pragma omp task
                {
                    string signature;
                    GenerateSignature(sequence.value, signature);

                    if (use_regex) {
                        smatch match;
                        regex_search(sequence.id, match, seqid_regex);
                        if (match.empty()) {
                            #pragma omp critical(exit)
                            {
                                cerr << "ABORTED: Regular expression resulted in empty match given the following sequence ID..." << endl;
                                cerr << sequence.id << endl;
                                exit(1);
                            }
                        }
                        sequence.id = string((match.size() > 1 ? match[1].str() : match[0].str()));
                    }

                    #pragma omp critical(write_out)
                    {
                        headfile << sequence.id << endl;
                        sigfile << signature;
                    }

                    #pragma omp critical(update_progress)
                    if (int((float)sequence_index++ / sequence_count * 100.0) != percent_complete) {
                        cerr << '\r' << ++percent_complete << '%';
                    }
                }
            });
        }
        cerr << "\r100%" << endl;

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
        cerr << time_span.count() << 's' << endl;

        return 0;
    }
    // --------------------------------------------------------------------------------------------

    // SEARCH -------------------------------------------------------------------------------------
    if (use == "search") {
        if (argc < 6) {
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        static string target_filepath;
        static vector<string> query_files;
        static ifstream target_sigfile;
        static ifstream target_headfile;
        static ofstream resultfile;

        static bool unique = false;
        static string result_format = "tsv";

        static bool use_regex = false;
        static regex seqid_regex;

        // Argument parsing
        for (int i = 2; i < argc; i++) {
            string arg = string(argv[i]);

            if (arg[0] == '-') {
                // Configuration
                string setting = arg.substr(1, arg.length());

                if (setting == "threshold") {
                    DIST_THRESHOLD = stod(argv[++i]);
                    continue;
                }

                if (setting == "top") {
                    RESULT_LIMIT = stoi(argv[++i]);
                    continue;
                }

                if (setting == "format") {
                    string format = argv[++i];

                    if (format != "tsv" && format != "csv" && format != "trec") {
                        cerr << "Unknown format: " << format << endl;
                        return 1;
                    }

                    result_format = format;
                    continue;
                }

                if (setting == "unique") {
                    unique = true;
                    continue;
                }

                if (setting == "match") {
                    use_regex = true;
                    seqid_regex.assign(argv[++i]);
                    continue;
                }

                if (setting == "o") {
                    resultfile.open(argv[++i]);
                    
                    if (!resultfile.is_open()) {
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

        if (!resultfile.is_open()) {
            cerr << "Please specifiy output with '-o'" << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        target_sigfile.open(target_filepath);

        if (!target_sigfile.is_open()) {
            cerr << "Error opening signature file: " << target_filepath << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        target_sigfile.close();

        target_headfile.open(target_filepath + ".head");

        if (!target_headfile.is_open()) {
            cerr << "Error opening header file: " << target_filepath + ".head" << endl;
            cerr << Usage(SEARCH) << endl;
            return 1;
        }

        string line;

        // Read and match signature generation metadata
        getline(target_headfile, line, ',');
        KMER_LEN = stoi(line);

        getline(target_headfile, line, ',');
        SIGNATURE_WIDTH = stoi(line);

        getline(target_headfile, line);
        SIGNATURE_DENSITY = stoi(line);

        target_headfile.close();

        cerr << "Searching...";
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        auto result_comparator = [](const SearchResult &a, const SearchResult &b) {
            return a.hamming_dist < b.hamming_dist;
        };

        #pragma omp parallel
        #pragma omp single
        for (string query_filepath : query_files) {
            ForEachSignature(query_filepath, [&](Signature &query_signature) {
                if (use_regex) {
                    smatch match;
                    regex_search(query_signature.id, match, seqid_regex);
                    if (match.empty()) {
                        #pragma omp critical(exit)
                        {
                            cerr << "ABORTED: Regular expression resulted in empty match given the following sequence ID..." << endl;
                            cerr << query_signature.id << endl;
                            exit(1);
                        }
                    }
                    query_signature.id = string((match.size() > 1 ? match[1].str() : match[0].str()));
                }

                #pragma omp task
                {
                    priority_queue<SearchResult, vector<SearchResult>, decltype(result_comparator)>
                        these_results(result_comparator);

                    ForEachSignature(target_filepath, [&](Signature &target_signature) {
                        if (use_regex) {
                            smatch match;
                            regex_search(target_signature.id, match, seqid_regex);
                            if (match.empty()) {
                                #pragma omp critical(exit)
                                {
                                    cerr << "ABORTED: Regular expression resulted in empty match given the following sequence ID..." << endl;
                                    cerr << target_signature.id << endl;
                                    exit(1);
                                }
                            }
                            target_signature.id = string((match.size() > 1 ? match[1].str() : match[0].str()));
                        }

                        if (unique && query_signature.id == target_signature.id) return;

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
                            if ((!RESULT_LIMIT || these_results.size() < RESULT_LIMIT)
                                    || hamming_dist < these_results.top().hamming_dist) {
                                these_results.emplace(
                                    target_signature.id,
                                    hamming_dist,
                                    normalised_dist
                                );
                                
                                if (RESULT_LIMIT && these_results.size() > RESULT_LIMIT) {
                                    these_results.pop();
                                }
                            }
                        }
                    });

                    deque<SearchResult> sorted_results;

                    while(!these_results.empty()) {
                        sorted_results.push_front(these_results.top());
                        these_results.pop();
                    }

                    // Output results based on format ...

                    if (result_format == "tsv") {
                        #pragma omp critical(write_out)
                        for(SearchResult &result : sorted_results) {
                            resultfile
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

                    if (result_format == "csv") {
                        #pragma omp critical(write_out)
                        for(SearchResult &result : sorted_results) {
                            resultfile
                                << query_signature.id
                                << ','
                                << result.target
                                << ','
                                << result.hamming_dist
                                << ','
                                << result.normalised_dist
                                << endl;
                        }
                    }

                    if (result_format == "trec") {
                        uint index = 1;
                        #pragma omp critical(write_out)
                        for(SearchResult &result : sorted_results) {
                            resultfile
                                << query_signature.id
                                << " Q0 "
                                << result.target
                                << ' ' << index++ << ' '
                                << result.normalised_dist
                                << " biosig" << endl;
                        }
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