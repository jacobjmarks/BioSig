#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>

// DEFAULTS ---------------------
static uint KMER_LEN = 5;
static uint SIGNATURE_LEN = 1024;
static uint SIGNATURE_DENSITY = 19;
// ------------------------------

using namespace std;

vector<int> hashKmer(string kmer) {
    srand(hash<string>{}(kmer));

    vector<int> kmerHash(SIGNATURE_LEN);

    uint pos = 0;
    uint set_limit = SIGNATURE_LEN / SIGNATURE_DENSITY / 2;

    // Fill half with 1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_LEN;

        if (!kmerHash[pos]) {
            kmerHash[pos] = 1;
            set++;
        }
    }

    // Fill other half with -1
    for (uint set = 0; set < set_limit;) {
        pos = rand() % SIGNATURE_LEN;

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

    vector<int> signature(SIGNATURE_LEN);

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

                    for (uint i = 0; i < SIGNATURE_LEN; i++) {
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
    for (uint i = 0; i < SIGNATURE_LEN; i++) {
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

    cerr << "Generating signature...";
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    generateSignature(argv[1]);

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cerr << time_span.count() << 's' << endl;

    return 0;
}