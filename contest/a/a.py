import sys


def edit_distance_banded(s: str, t: str, k: int = 100) -> int:
    """
    Compute edit distance between s and t using banded DP.
    Only computes values within k diagonals of the main diagonal.
    Time complexity: O((n + m) * k)
    """
    n: int = len(s)
    m: int = len(t)

    # If length difference exceeds k, distance > k (but problem guarantees <= 100)
    if abs(n - m) > k:
        return k + 1

    # Pre-allocate both arrays at maximum band size to avoid repeated allocations
    max_band_size: int = 2 * k + 1
    prev: list[int] = [0] * max_band_size
    curr: list[int] = [0] * max_band_size

    # Initialize for row 0: dp[0][j] = j for j in range(0, min(m, k) + 1)
    for j in range(min(m, k) + 1):
        prev[j] = j

    for i in range(1, n + 1):
        # Column range for row i
        j_min: int = max(0, i - k)
        j_max: int = min(m, i + k)

        for j in range(j_min, j_max + 1):
            idx: int = j - j_min  # Index in curr array

            if j == 0:
                curr[idx] = i
            else:
                # Previous row index for column j
                prev_j_min: int = max(0, (i - 1) - k)
                prev_j_max: int = min(m, (i - 1) + k)

                cost: int = 0 if s[i - 1] == t[j - 1] else 1

                candidates: list[int] = []

                # Substitution/match: dp[i-1][j-1] + cost
                if prev_j_min <= j - 1 <= prev_j_max:
                    candidates.append(prev[(j - 1) - prev_j_min] + cost)

                # Deletion from s: dp[i-1][j] + 1
                if prev_j_min <= j <= prev_j_max:
                    candidates.append(prev[j - prev_j_min] + 1)

                # Insertion to s: dp[i][j-1] + 1
                if j - 1 >= j_min:
                    candidates.append(curr[idx - 1] + 1)

                curr[idx] = min(candidates) if candidates else k + 1

        # Swap arrays by swapping references
        prev, curr = curr, prev

    # Result is dp[n][m]
    final_j_min: int = max(0, n - k)
    if final_j_min <= m <= min(m, n + k):
        return prev[m - final_j_min]
    return k + 1


def main() -> None:
    input_data: list[str] = sys.stdin.read().split()
    s: str = input_data[0]
    t: str = input_data[1]
    print(edit_distance_banded(s, t, 100))


if __name__ == "__main__":
    main()
