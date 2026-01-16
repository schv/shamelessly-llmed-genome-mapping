#pragma once

#include <vector>
#include <string>
#include <algorithm>

namespace bio {

inline std::vector<int> buildSuffixArray(const std::string& s) {
    int n = s.size();
    std::vector<int> sa(n), rk(n), tmp(n);
    
    for (int i = 0; i < n; i++) {
        sa[i] = i;
        rk[i] = (unsigned char)s[i];
    }
    
    for (int k = 1; k < n; k *= 2) {
        auto cmp = [&](int a, int b) {
            if (rk[a] != rk[b]) return rk[a] < rk[b];
            int ra = a + k < n ? rk[a + k] : -1;
            int rb = b + k < n ? rk[b + k] : -1;
            return ra < rb;
        };
        std::sort(sa.begin(), sa.end(), cmp);
        
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            tmp[sa[i]] = tmp[sa[i-1]] + (cmp(sa[i-1], sa[i]) ? 1 : 0);
        }
        rk = tmp;
        if (rk[sa[n-1]] == n - 1) break;
    }
    return sa;
}

// Find lower bound: first suffix >= pattern
inline int suffixArrayLowerBound(const std::string& s, const std::vector<int>& sa, const std::string& pat) {
    int lo = 0, hi = sa.size();
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (s.compare(sa[mid], pat.size(), pat) < 0)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

// Find upper bound: first suffix > pattern
inline int suffixArrayUpperBound(const std::string& s, const std::vector<int>& sa, const std::string& pat) {
    int lo = 0, hi = sa.size();
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (s.compare(sa[mid], pat.size(), pat) <= 0)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

// Find all occurrences of pattern in text using suffix array
inline std::vector<int> findAllOccurrences(const std::string& text, const std::vector<int>& sa, const std::string& pattern) {
    int lo = suffixArrayLowerBound(text, sa, pattern);
    int hi = suffixArrayUpperBound(text, sa, pattern);
    std::vector<int> result;
    for (int i = lo; i < hi; i++) {
        result.push_back(sa[i]);
    }
    return result;
}

// Check if pattern has unique occurrence
inline bool hasUniqueMatch(const std::string& text, const std::vector<int>& sa, const std::string& pattern) {
    int lo = suffixArrayLowerBound(text, sa, pattern);
    int hi = suffixArrayUpperBound(text, sa, pattern);
    return (hi - lo) == 1;
}

// Get unique match position, returns -1 if not unique
inline int getUniqueMatchPosition(const std::string& text, const std::vector<int>& sa, const std::string& pattern) {
    int lo = suffixArrayLowerBound(text, sa, pattern);
    int hi = suffixArrayUpperBound(text, sa, pattern);
    if (hi - lo == 1) return sa[lo];
    return -1;
}

} // namespace bio
