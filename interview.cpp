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

