#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace bio {

// Band-limited edit distance calculation
// Optimized for cases where edit distance is guaranteed to be small (≤ maxDist)
// Returns edit distance, or maxDist+1 if distance exceeds maxDist
template<int maxDist = 100>
inline int editDistance(const std::string& s, const std::string& t) {
    int n = s.size();
    int m = t.size();
    
    constexpr int W = 2 * maxDist + 1;
    constexpr int INF = maxDist + 1;
    
    int dp[2][W];
    
    for (int j = 0; j < W; ++j) {
        dp[0][j] = INF;
        dp[1][j] = INF;
    }
    dp[0][maxDist] = 0;
    
    for (int i = 0; i <= n; ++i) {
        int curr = i & 1;
        int prev = 1 - curr;
        
        for (int j = 0; j < W; ++j) dp[curr][j] = INF;
        
        int j_min = std::max(0, i - maxDist);
        int j_max = std::min(m, i + maxDist);
        
        for (int j = j_min; j <= j_max; ++j) {
            int d = j - i;
            int idx = d + maxDist;
            
            if (i == 0 && j == 0) {
                dp[curr][idx] = 0;
            } else if (i == 0) {
                dp[curr][idx] = j;
            } else if (j == 0) {
                dp[curr][idx] = i;
            } else {
                // Substitution or match
                dp[curr][idx] = dp[prev][idx] + (s[i-1] != t[j-1] ? 1 : 0);
                
                // Insertion
                if (idx > 0 && j > 0) {
                    dp[curr][idx] = std::min(dp[curr][idx], dp[curr][idx-1] + 1);
                }
                // Deletion
                if (idx < W - 1) {
                    dp[curr][idx] = std::min(dp[curr][idx], dp[prev][idx+1] + 1);
                }
            }
        }
    }
    
    return dp[n & 1][m - n + maxDist];
}

// Standard edit distance without band limitation
// Use for general cases where distance may be large
inline int editDistanceFull(const std::string& s, const std::string& t) {
    int n = s.size();
    int m = t.size();
    
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));
    
    for (int i = 0; i <= n; i++) dp[i][0] = i;
    for (int j = 0; j <= m; j++) dp[0][j] = j;
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            if (s[i-1] == t[j-1]) {
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = 1 + std::min({dp[i-1][j], dp[i][j-1], dp[i-1][j-1]});
            }
        }
    }
    
    return dp[n][m];
}

// Check if two strings are within a given edit distance
template<int maxDist = 100>
inline bool withinEditDistance(const std::string& s, const std::string& t, int threshold) {
    if (std::abs((int)s.size() - (int)t.size()) > threshold) return false;
    return editDistance<maxDist>(s, t) <= threshold;
}

} // namespace bio
