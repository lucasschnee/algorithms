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






struct Treap{ /// hash = 96814
    int len;
    const int ADD = 1000010;
    const int MAXVAL = 1000000010;
    tr1::unordered_map <long long, int> mp; /// Change to int if only int in treap
    tree<long long, null_type, less<long long>, rb_tree_tag, tree_order_statistics_node_update> T;

    Treap(){
        len = 0;
        T.clear(), mp.clear();
    }

    inline void clear(){
        len = 0;
        T.clear(), mp.clear();
    }

    inline void insert(long long x){
        len++, x += MAXVAL;
        int c = mp[x]++;
        T.insert((x * ADD) + c);
    }

    inline void erase(long long x){
        x += MAXVAL;
        int c = mp[x];
        if (c){
            c--, mp[x]--, len--;
            T.erase((x * ADD) + c);
        }
    }

    /// 1-based index, returns the K'th element in the treap, -1 if none exists
    inline long long kth(int k){
        if (k < 1 || k > len) return -1;
        auto it = T.find_by_order(--k);
        return ((*it) / ADD) - MAXVAL;
    }

    /// Count of value < x in treap
    inline int count(long long x){
        x += MAXVAL;
        int c = mp[--x];
        return (T.order_of_key((x * ADD) + c));
    }

    /// Number of elements in treap
    inline int size(){
        return len;
    }
};

