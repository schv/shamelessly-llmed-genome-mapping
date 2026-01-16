#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "lib/bio.hpp"

using namespace std;

// Parse FASTA file - concatenate all sequences
string loadFasta(const string& filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Error: Cannot open " << filename << endl;
        exit(1);
    }
    
    string genome, line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        genome += line;
    }
    return genome;
}

// FASTQ read structure
struct Read {
    string id;
    string seq;
    string qual;
};

// Parse single read from FASTQ (4 lines)
bool readFastq(ifstream& in, Read& read) {
    string plus;
    if (!getline(in, read.id)) return false;
    if (!getline(in, read.seq)) return false;
    if (!getline(in, plus)) return false;
    if (!getline(in, read.qual)) return false;
    if (!read.id.empty() && read.id[0] == '@') {
        read.id = read.id.substr(1);
    }
    return true;
}

// Mapping result
enum class MapStatus { Unmapped, Unique, Multi };

struct MappingResult {
    MapStatus status;
    int position;
    int edit_dist;
};

// Map single read using seed-and-extend
MappingResult mapRead(const string& genome, const vector<int>& sa,
                      const string& read, int seed_len, int max_errors) {
    MappingResult result{MapStatus::Unmapped, -1, -1};
    
    // Skip reads starting with N (common Illumina artifact)
    if (!read.empty() && read[0] == 'N') return result;
    
    // Try exact match first (fast path)
    int lo = bio::suffixArrayLowerBound(genome, sa, read);
    int hi = bio::suffixArrayUpperBound(genome, sa, read);
    
    if (hi > lo) {
        if (hi - lo == 1) {
            result.status = MapStatus::Unique;
            result.position = sa[lo];
            result.edit_dist = 0;
        } else {
            result.status = MapStatus::Multi;
            result.position = sa[lo];
            result.edit_dist = 0;
        }
        return result;
    }
    
    // Seed-and-extend: try multiple seeds
    vector<int> candidates;
    int num_seeds = 3;
    int step = (read.size() - seed_len) / max(1, num_seeds - 1);
    
    for (int i = 0; i < num_seeds && i * step + seed_len <= (int)read.size(); i++) {
        string seed = read.substr(i * step, seed_len);
        
        // Skip seeds with N
        if (seed.find('N') != string::npos) continue;
        
        int slo = bio::suffixArrayLowerBound(genome, sa, seed);
        int shi = bio::suffixArrayUpperBound(genome, sa, seed);
        
        // Limit candidates per seed to avoid explosion
        int max_hits = 100;
        for (int j = slo; j < shi && j < slo + max_hits; j++) {
            int genome_start = sa[j] - i * step;
            if (genome_start >= 0 && genome_start + (int)read.size() <= (int)genome.size()) {
                candidates.push_back(genome_start);
            }
        }
    }
    
    if (candidates.empty()) return result;
    
    // Remove duplicates
    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    
    // Verify candidates with edit distance
    int best_dist = max_errors + 1;
    int best_pos = -1;
    int best_count = 0;
    
    for (int cand : candidates) {
        string ref_seg = genome.substr(cand, read.size());
        int dist = bio::editDistance<10>(ref_seg, read);
        
        if (dist < best_dist) {
            best_dist = dist;
            best_pos = cand;
            best_count = 1;
        } else if (dist == best_dist && cand != best_pos) {
            best_count++;
        }
    }
    
    if (best_dist <= max_errors) {
        if (best_count == 1) {
            result.status = MapStatus::Unique;
        } else {
            result.status = MapStatus::Multi;
        }
        result.position = best_pos;
        result.edit_dist = best_dist;
    }
    
    return result;
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    string genome_file = "data/GCF_000005845.2_ASM584v2_genomic.fna";
    string reads_file = "data/ERR022075_1.fastq";
    int max_reads = -1;  // -1 = all reads
    int seed_len = 20;
    int max_errors = 3;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-g" && i + 1 < argc) genome_file = argv[++i];
        else if (arg == "-r" && i + 1 < argc) reads_file = argv[++i];
        else if (arg == "-n" && i + 1 < argc) max_reads = stoi(argv[++i]);
        else if (arg == "-s" && i + 1 < argc) seed_len = stoi(argv[++i]);
        else if (arg == "-e" && i + 1 < argc) max_errors = stoi(argv[++i]);
        else if (arg == "-h") {
            cerr << "Usage: " << argv[0] << " [options]\n"
                 << "  -g <file>  Reference genome (FASTA)\n"
                 << "  -r <file>  Reads file (FASTQ)\n"
                 << "  -n <num>   Max reads to process (-1 = all)\n"
                 << "  -s <len>   Seed length (default: 20)\n"
                 << "  -e <num>   Max errors allowed (default: 3)\n";
            return 0;
        }
    }
    
    auto start_time = chrono::high_resolution_clock::now();
    
    // Load reference genome
    cerr << "Loading reference genome..." << endl;
    string genome = loadFasta(genome_file);
    cerr << "Genome size: " << genome.size() << " bp" << endl;
    
    // Build suffix array
    cerr << "Building suffix array..." << endl;
    auto sa_start = chrono::high_resolution_clock::now();
    vector<int> sa = bio::buildSuffixArray(genome);
    auto sa_end = chrono::high_resolution_clock::now();
    cerr << "Suffix array built in " 
         << chrono::duration_cast<chrono::milliseconds>(sa_end - sa_start).count() 
         << " ms" << endl;
    
    // Open reads file
    ifstream reads_in(reads_file);
    if (!reads_in) {
        cerr << "Error: Cannot open " << reads_file << endl;
        return 1;
    }
    
    // Mapping statistics
    long long total_reads = 0;
    long long mapped_reads = 0;
    long long unique_mapped = 0;
    long long multi_mapped = 0;
    long long total_edit_dist = 0;
    vector<int> coverage(genome.size(), 0);
    
    cerr << "Mapping reads..." << endl;
    Read read;
    int progress_interval = 100000;
    
    while (readFastq(reads_in, read)) {
        if (max_reads >= 0 && total_reads >= max_reads) break;
        total_reads++;
        
        MappingResult result = mapRead(genome, sa, read.seq, seed_len, max_errors);
        
        if (result.status != MapStatus::Unmapped) {
            mapped_reads++;
            total_edit_dist += result.edit_dist;
            
            if (result.status == MapStatus::Unique) {
                unique_mapped++;
                // Update coverage
                int start = result.position;
                int end = min(start + (int)read.seq.size(), (int)genome.size());
                for (int i = start; i < end; i++) {
                    coverage[i]++;
                }
            } else {
                multi_mapped++;
            }
        }
        
        if (total_reads % progress_interval == 0) {
            cerr << "\rProcessed " << total_reads << " reads... " 
                 << (100.0 * mapped_reads / total_reads) << "% mapped" << flush;
        }
    }
    cerr << endl;
    
    auto end_time = chrono::high_resolution_clock::now();
    double total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    
    // Calculate coverage statistics
    long long covered_bases = 0;
    long long total_coverage = 0;
    for (int c : coverage) {
        if (c > 0) covered_bases++;
        total_coverage += c;
    }
    
    // Output report
    cout << "=== Genome Mapping Report ===" << endl;
    cout << endl;
    cout << "Algorithms used:" << endl;
    cout << "  - Suffix array O(n log^2 n) construction" << endl;
    cout << "  - Seed-and-extend with " << seed_len << "-mer seeds" << endl;
    cout << "  - Band-limited edit distance (max " << max_errors << " errors)" << endl;
    cout << endl;
    cout << "Reference: " << genome_file << endl;
    cout << "Genome size: " << genome.size() << " bp" << endl;
    cout << endl;
    cout << "Reads file: " << reads_file << endl;
    cout << "Total reads processed: " << total_reads << endl;
    cout << endl;
    cout << "Mapping statistics:" << endl;
    cout << "  Mapped reads: " << mapped_reads 
         << " (" << fixed << setprecision(2) << (100.0 * mapped_reads / total_reads) << "%)" << endl;
    cout << "  Unmapped reads: " << (total_reads - mapped_reads)
         << " (" << fixed << setprecision(2) << (100.0 * (total_reads - mapped_reads) / total_reads) << "%)" << endl;
    cout << endl;
    cout << "  Uniquely mapped: " << unique_mapped 
         << " (" << fixed << setprecision(2) << (100.0 * unique_mapped / total_reads) << "%)" << endl;
    cout << "  Multi-mapped: " << multi_mapped
         << " (" << fixed << setprecision(2) << (100.0 * multi_mapped / total_reads) << "%)" << endl;
    cout << endl;
    cout << "Alignment quality:" << endl;
    cout << "  Average edit distance: " << fixed << setprecision(2) 
         << (mapped_reads > 0 ? (double)total_edit_dist / mapped_reads : 0) << endl;
    cout << endl;
    cout << "Genome coverage (from uniquely mapped reads):" << endl;
    cout << "  Covered bases: " << covered_bases 
         << " (" << fixed << setprecision(2) << (100.0 * covered_bases / genome.size()) << "%)" << endl;
    cout << "  Average depth: " << fixed << setprecision(2) 
         << (double)total_coverage / genome.size() << "x" << endl;
    cout << endl;
    cout << "Total runtime: " << fixed << setprecision(1) << total_time << " seconds" << endl;
    
    return 0;
}
