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

