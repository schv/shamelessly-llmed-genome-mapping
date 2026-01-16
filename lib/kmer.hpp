#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

namespace bio {

namespace detail {
    inline unsigned long long dnaCharValue(char c) {
        switch(c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
        }
        return 0;
    }
    
    constexpr unsigned long long KMER_HASH_BASE = 5;
    constexpr unsigned long long KMER_HASH_MOD = 1000000000000000003ULL;
}

// Compute rolling hash for a k-mer
inline unsigned long long computeKmerHash(const std::string& kmer) {
    unsigned long long hash = 0;
    for (char c : kmer) {
        hash = (hash * detail::KMER_HASH_BASE + detail::dnaCharValue(c)) % detail::KMER_HASH_MOD;
    }
    return hash;
}

// Count frequency of all k-mers in a string
// Returns map from k-mer hash to (frequency, first_position)
inline std::unordered_map<unsigned long long, std::pair<int, int>> countKmers(const std::string& s, int k) {
    std::unordered_map<unsigned long long, std::pair<int, int>> result;
    
    int n = s.size();
    if (k > n) return result;
    
    unsigned long long base_pow = 1;
    for (int i = 0; i < k - 1; i++) {
        base_pow = (base_pow * detail::KMER_HASH_BASE) % detail::KMER_HASH_MOD;
    }
    
    // Initial hash
    unsigned long long hash = 0;
    for (int i = 0; i < k; i++) {
        hash = (hash * detail::KMER_HASH_BASE + detail::dnaCharValue(s[i])) % detail::KMER_HASH_MOD;
    }
    result[hash] = {1, 0};
    
    // Rolling hash
    for (int i = 1; i <= n - k; i++) {
        hash = (hash + detail::KMER_HASH_MOD - (detail::dnaCharValue(s[i-1]) * base_pow) % detail::KMER_HASH_MOD) % detail::KMER_HASH_MOD;
        hash = (hash * detail::KMER_HASH_BASE + detail::dnaCharValue(s[i+k-1])) % detail::KMER_HASH_MOD;
        
        auto it = result.find(hash);
        if (it == result.end()) {
            result[hash] = {1, i};
        } else {
            it->second.first++;
        }
    }
    
    return result;
}

// Find the most frequent k-mer in a string
// Returns pair of (k-mer string, frequency)
inline std::pair<std::string, int> findMostFrequentKmer(const std::string& s, int k) {
    int n = s.size();
    if (k > n) return {s, 1};
    
    unsigned long long base_pow = 1;
    for (int i = 0; i < k - 1; i++) {
        base_pow = (base_pow * detail::KMER_HASH_BASE) % detail::KMER_HASH_MOD;
    }
    
    unsigned long long hash = 0;
    for (int i = 0; i < k; i++) {
        hash = (hash * detail::KMER_HASH_BASE + detail::dnaCharValue(s[i])) % detail::KMER_HASH_MOD;
    }
    
    std::unordered_map<unsigned long long, int> freq;
    std::unordered_map<unsigned long long, int> first_pos;
    
    freq[hash] = 1;
    first_pos[hash] = 0;
    
    int max_freq = 1;
    unsigned long long best_hash = hash;
    int best_pos = 0;
    
    for (int i = 1; i <= n - k; i++) {
        hash = (hash + detail::KMER_HASH_MOD - (detail::dnaCharValue(s[i-1]) * base_pow) % detail::KMER_HASH_MOD) % detail::KMER_HASH_MOD;
        hash = (hash * detail::KMER_HASH_BASE + detail::dnaCharValue(s[i+k-1])) % detail::KMER_HASH_MOD;
        freq[hash]++;
        
        if (first_pos.find(hash) == first_pos.end()) {
            first_pos[hash] = i;
        }
        
        if (freq[hash] > max_freq) {
            max_freq = freq[hash];
            best_hash = hash;
            best_pos = first_pos[hash];
        }
    }
    
    return {s.substr(best_pos, k), max_freq};
}

// Get all k-mers from a string as vector of strings
inline std::vector<std::string> extractKmers(const std::string& s, int k) {
    std::vector<std::string> result;
    int n = s.size();
    if (k > n) return result;
    
    result.reserve(n - k + 1);
    for (int i = 0; i <= n - k; i++) {
        result.push_back(s.substr(i, k));
    }
    return result;
}

} // namespace bio
