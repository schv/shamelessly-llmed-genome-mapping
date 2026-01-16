#pragma once

// Genome Mapping Library
// Algorithms for DNA sequence analysis

#include "suffix_array.hpp"
#include "bwt.hpp"
#include "kmer.hpp"
#include "edit_distance.hpp"

// Library namespace: bio
//
// Available functions:
//
// suffix_array.hpp:
//   - buildSuffixArray(s)              : O(n log² n) suffix array construction
//   - suffixArrayLowerBound(s, sa, p)  : binary search lower bound
//   - suffixArrayUpperBound(s, sa, p)  : binary search upper bound  
//   - findAllOccurrences(s, sa, p)     : find all pattern occurrences
//   - hasUniqueMatch(s, sa, p)         : check for unique match
//   - getUniqueMatchPosition(s, sa, p) : get position of unique match
//
// bwt.hpp:
//   - computeBWT(s)                    : Burrows-Wheeler Transform
//   - inverseBWT(bwt)                  : inverse BWT
//   - buildOccurrenceTable(bwt)        : FM-index occurrence table
//   - buildCumulativeCounts(bwt)       : FM-index C array
//
// kmer.hpp:
//   - computeKmerHash(kmer)            : rolling hash for k-mer
//   - countKmers(s, k)                 : count all k-mer frequencies
//   - findMostFrequentKmer(s, k)       : find most frequent k-mer
//   - extractKmers(s, k)               : extract all k-mers as strings
//
// edit_distance.hpp:
//   - editDistance<maxDist>(s, t)      : band-limited edit distance
//   - editDistanceFull(s, t)           : standard edit distance
//   - withinEditDistance<maxDist>(s, t, threshold) : check distance threshold
