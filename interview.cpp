// Author: Lucas
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <deque>
#include <iomanip>



typedef long long ll;
ll MOD=1000000007;
ll power(ll base, ll exp) {
	ll res = 1;
	base %= MOD;
	while (exp>0) {
		if (exp%2==1) {
			res *= base;
			res %= MOD;
		}

		base = base * base;
		base %= MOD;
		exp = exp/2;
	}
	return res;
}

ll mod_inv(ll base) {
	return power(base, MOD-2);
}




static const long long MOD = 1e9 + 7;

long long modpow(long long a, long long b) {
  long long res = 1;
  while (b) {
      if (b & 1) res = (res * a) % MOD;
      a = (a * a) % MOD;
      b >>= 1;
  }
  return res;
}

long long comb(int n, int k, vector<long long>& fact, vector<long long>& invfact) {
  if (n < k) return 0;
  return fact[n] * invfact[k] % MOD * invfact[n - k] % MOD;
}

vector<long long> fact(N+1, 1), invfact(N+1, 1);

for (int i = 1; i <= N; ++i) {
   fact[i] = fact[i - 1] * i % MOD;
}


invfact[N] = modpow(fact[N], MOD - 2);
for (int i = N-1; i >= 0; --i) {
   invfact[i] = invfact[i + 1] * (i + 1) % MOD;
}


class UnionFind {
public:
    int count;
    vector<int> parent;
    vector<int> rank;   

    UnionFind(int N) {
        count = N;
        parent.resize(N);
        rank.assign(N, 1);

        for (int i = 0; i < N; ++i) {
            parent[i] = i;
        }
    }

    int find(int p) {
        if (p != parent[p]) {
            parent[p] = find(parent[p]);  
        }
        return parent[p];
    }

    bool unite(int p, int q) {
        int prt = find(p);
        int qrt = find(q);

        if (prt == qrt) return false;

      
        if (rank[prt] >= rank[qrt]) {
            parent[qrt] = prt;
            rank[prt] += rank[qrt];
        } else {
            parent[prt] = qrt;
            rank[qrt] += rank[prt];
        }

        count--;
        return true;
    }
};

set<int> primes;
vector<int> primes_list;


void buildPrimes() {
    const int N = 200000;
    vector<bool> isPrime(N + 1, true);
    isPrime[0] = isPrime[1] = false;

    for (int i = 2; i * i <= N; i++) {
        if (isPrime[i]) {
            for (int j = i * i; j <= N; j += i) {
                isPrime[j] = false;
            }
        }
    }

    for (int i = 2; i <= N; i++) {
        if (isPrime[i]) {
            primes.insert(i);
            primes_list.push_back(i);
        }
    }
}

#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds; 

template <typename num_t>
using ordered_set = tree<num_t, null_type, less<num_t>, rb_tree_tag, tree_order_statistics_node_update>;

template <typename num_t>
struct ordered_multiset {
	ordered_set<pair<num_t, int> > vals;
	set<pair<num_t, int> > best; /* start at -1 */
	
	/* helper, find the lowest value that represents the element */
	int findbest(num_t val) {
		return (*best.upper_bound(make_pair(val - 1, 0))).second;
	}
	
	/* is element in set */
	bool contains(num_t val) {
		return vals.find(make_pair(val, -1)) != vals.end();
	}
	
	void insert(num_t val) {
		if (contains(val)) { /* already in, update lowest value and insert a new one */
			int loc = findbest(val);
			best.erase(make_pair(val, loc));
			best.insert(make_pair(val, loc - 1));
			vals.insert(make_pair(val, loc - 1));
		} else { /* make lowest value -1 and insert it */
			best.insert(make_pair(val, -1));
			vals.insert(make_pair(val, -1));
		}
	}
	
	void erase(num_t val) { /* erases one */
		if (!contains(val)) return; /* not in */
		num_t loc = findbest(val);
		
		/* remove the element and its best */
		best.erase(make_pair(val, loc));
		vals.erase(make_pair(val, loc));
		if (loc != -1) best.insert(make_pair(val, loc + 1)); /* more elements in set, update best */
	}
	
	/* unmodified functions */
	num_t find_by_order(int k) { return (*vals.find_by_order(k)).first; }
	int order_of_key(num_t k) { return vals.order_of_key(make_pair(k - 1, 0)); }
	auto begin() { return vals.begin(); }
	auto end() { return vals.end(); }
	auto rbegin() { return vals.rbegin(); }
	auto rend() { return vals.rend(); }
	int size() { return vals.size(); }
	void clear() { vals.clear(); best.clear(); }
	int count(num_t k) { return vals.order_of_key({k, 0}) - vals.order_of_key({k - 1, 0}); }
	auto lower_bound(num_t k) { return vals.lower_bound(make_pair(k - 1, 0)); }
	auto upper_bound(num_t k) { return vals.upper_bound(make_pair(k, 0)); }
};





// Matrix Multiply AB,  % MOD
    vector<vector<long long>> multiply(const vector<vector<long long>>& A, const vector<vector<long long>>& B) {
        int sz = A.size();
        vector<vector<long long>> C(sz, vector<long long>(sz, 0));
        for (int i = 0; i < sz; ++i) {
            for (int k = 0; k < sz; ++k) { 
                if (A[i][k] == 0) continue; 
                for (int j = 0; j < sz; ++j) {
                    C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % MOD;
                }
            }
        }
        return C;
    }

// Binary Exponentiation: M^p, log(sz))
vector<vector<long long>> power(vector<vector<long long>> M, long long p) {
	int sz = M.size();
	vector<vector<long long>> res(sz, vector<long long>(sz, 0));
	for (int i = 0; i < sz; ++i) res[i][i] = 1; // Identity matrix
	while (p > 0) {
		if (p & 1) res = multiply(res, M);
		M = multiply(M, M);
		p >>= 1;
	}
	return res;
}




class MaxSegmentTree {
private:
    int n;
    vector<long long> tree;

public:
    MaxSegmentTree(int size) {
        n = 1;
        while (n < size) n <<= 1;
        tree.assign(2 * n, 0LL); 
    }

    /**
     * Updates value at index 'idx' with 'val'
     */
    void update(int idx, long long val) {
        idx += n; 
        if (tree[idx] >= val) return; 

        tree[idx] = val;
        for (idx >>= 1; idx >= 1; idx >>= 1) {
            tree[idx] = max(tree[2 * idx], tree[2 * idx + 1]);
        }
    }

    /**
     * Queries the maximum value in the range [l, r] inclusive
     */
    long long query_max(int l, int r) {
        if (l > r || l < 0) return 0LL;
        
        long long res = 0;
        l += n;
        r += n;

        while (l <= r) {
            if (l % 2 == 1) res = max(res, tree[l++]);
            if (r % 2 == 0) res = max(res, tree[r--]);
            l /= 2;
            r /= 2;
        }
        return res;
    }
};
