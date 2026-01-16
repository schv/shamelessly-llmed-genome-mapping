#pragma once

#include <vector>
#include <string>
#include <algorithm>

namespace bio {

// Compute Burrows-Wheeler Transform of a string
// Appends '$' as sentinel and returns BWT
inline std::string computeBWT(const std::string& input) {
    std::string s = input + '$';
    int n = s.size();
    
    std::vector<int> sa(n), rk(n), tmp(n), cnt(std::max(256, n));
    
    for (int i = 0; i < n; i++) rk[i] = (unsigned char)s[i];
    
    for (int k = 1; k < n; k *= 2) {
        std::fill(cnt.begin(), cnt.end(), 0);
        for (int i = 0; i < n; i++) cnt[(i + k < n) ? rk[i + k] + 1 : 0]++;
        for (int i = 1; i < (int)cnt.size(); i++) cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--) tmp[--cnt[(i + k < n) ? rk[i + k] + 1 : 0]] = i;
        
        std::fill(cnt.begin(), cnt.end(), 0);
        for (int i = 0; i < n; i++) cnt[rk[i] + 1]++;
        for (int i = 1; i < (int)cnt.size(); i++) cnt[i] += cnt[i - 1];
        for (int i = 0; i < n; i++) sa[cnt[rk[tmp[i]]]++] = tmp[i];
        
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            int a = sa[i - 1], b = sa[i];
            bool eq = rk[a] == rk[b] && ((a + k < n ? rk[a + k] : -1) == (b + k < n ? rk[b + k] : -1));
            tmp[sa[i]] = tmp[sa[i - 1]] + !eq;
        }
        rk = tmp;
        if (rk[sa[n - 1]] == n - 1) break;
    }
    
    std::string bwt(n, ' ');
    for (int i = 0; i < n; i++) bwt[i] = s[(sa[i] + n - 1) % n];
    
    return bwt;
}

namespace detail {
    inline int dnaCharToIdx(char c) {
        if (c == '$') return 0;
        if (c == 'A') return 1;
        if (c == 'C') return 2;
        if (c == 'G') return 3;
        return 4; // T
    }
}

// Inverse Burrows-Wheeler Transform
// Input: BWT string containing $, A, C, G, T
// Output: original DNA string (without $)
inline std::string inverseBWT(const std::string& bwt) {
    int n = bwt.size();
    
    // Count occurrences of each character
    std::vector<int> cnt(5, 0);
    for (char c : bwt) cnt[detail::dnaCharToIdx(c)]++;
    
    // Cumulative counts (C array)
    std::vector<int> C(5);
    C[0] = 0;
    for (int i = 1; i < 5; i++) C[i] = C[i-1] + cnt[i-1];
    
    // Build rank array
    std::vector<int> rank(n);
    std::vector<int> seen(5, 0);
    for (int i = 0; i < n; i++) {
        int idx = detail::dnaCharToIdx(bwt[i]);
        rank[i] = seen[idx]++;
    }
    
    // Find position of $
    int pos = 0;
    for (int i = 0; i < n; i++) {
        if (bwt[i] == '$') { pos = i; break; }
    }
    
    // Reconstruct original string
    std::string result;
    result.reserve(n - 1);
    
    for (int i = 0; i < n - 1; i++) {
        pos = C[detail::dnaCharToIdx(bwt[pos])] + rank[pos];
        result += bwt[pos];
    }
    
    std::reverse(result.begin(), result.end());
    return result;
}

// FM-Index helper: build occurrence table for BWT
// Returns occ[c][i] = count of character c in bwt[0..i-1]
inline std::vector<std::vector<int>> buildOccurrenceTable(const std::string& bwt) {
    int n = bwt.size();
    std::vector<std::vector<int>> occ(5, std::vector<int>(n + 1, 0));
    
    for (int i = 0; i < n; i++) {
        for (int c = 0; c < 5; c++) {
            occ[c][i + 1] = occ[c][i];
        }
        occ[detail::dnaCharToIdx(bwt[i])][i + 1]++;
    }
    return occ;
}

// Get cumulative character counts (C array) for FM-index
inline std::vector<int> buildCumulativeCounts(const std::string& bwt) {
    std::vector<int> cnt(5, 0);
    for (char c : bwt) cnt[detail::dnaCharToIdx(c)]++;
    
    std::vector<int> C(5);
    C[0] = 0;
    for (int i = 1; i < 5; i++) C[i] = C[i-1] + cnt[i-1];
    return C;
}

} // namespace bio
