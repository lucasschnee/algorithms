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

// PBDS for SortedList equivalent
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

#define ll long long
typedef long double ld;

const ll INF = 1e18;
const int MOD = 1e9 + 7;

// PBDS Ordered Set (Equivalent to Python's SortedList)
// Supports: find_by_order(k) and order_of_key(x)
template<class T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<class T> using ordered_multiset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

// Directions
vector<pair<int, int>> directions = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}};

/* =========================================
   MATH & COMBINATORICS
========================================= */
bool is_prime(ll n) {
    if (n <= 1) return false;
    for (ll i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

vector<ll> get_divisors(ll n) {
    vector<ll> divs;
    for (ll i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            divs.push_back(i);
            if (i * i != n) divs.push_back(n / i);
        }
    }
    return divs;
}

// Combinatorics Precomputation
const int MX = 200005;
ll fac[MX], ifac[MX];

ll power(ll base, ll exp) {
    ll res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

void precompute_comb() {
    fac[0] = 1;
    ifac[0] = 1;
    for (int i = 1; i < MX; i++) {
        fac[i] = (fac[i - 1] * i) % MOD;
    }
    ifac[MX - 1] = power(fac[MX - 1], MOD - 2);
    for (int i = MX - 2; i >= 1; i--) {
        ifac[i] = (ifac[i + 1] * (i + 1)) % MOD;
    }
}

ll n_choose_k(int n, int k) {
    if (k < 0 || k > n) return 0;
    return fac[n] * ifac[n - k] % MOD * ifac[k] % MOD;
}

// Matrix Exponentiation
typedef vector<vector<ll>> matrix;
matrix mat_mult(matrix& A, matrix& B) {
    int size = A.size();
    matrix result(size, vector<ll>(size, 0));
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                result[i][j] = (result[i][j] + A[i][k] * B[k][j]) % MOD;
    return result;
}

matrix mat_pow(matrix mat, ll exp) {
    int size = mat.size();
    matrix result(size, vector<ll>(size, 0));
    for (int i = 0; i < size; i++) result[i][i] = 1;
    while (exp > 0) {
        if (exp % 2 == 1) result = mat_mult(result, mat);
        mat = mat_mult(mat, mat);
        exp /= 2;
    }
    return result;
}

/* =========================================
   STRINGS
========================================= */
class RabinKarp {
private:
    ll m = 1111111111111111111LL;
    ll p = 1000000007LL;
    vector<ll> pows, roll;
public:
    RabinKarp(const string& s) {
        int n = s.length();
        pows.assign(n + 1, 1);
        roll.assign(n + 1, 0);
        for (int i = 0; i < n; i++) {
            pows[i + 1] = (__int128(pows[i]) * p) % m;
            roll[i + 1] = (__int128(roll[i]) * p + s[i]) % m;
        }
    }
    ll query(int i, int j) {
        return (roll[j] - (__int128(roll[i]) * pows[j - i]) % m + m) % m;
    }
};

vector<int> kmp(const string& pattern, const string& text) {
    int m = pattern.length();
    vector<int> lps(m, 0);
    for (int i = 1, k = 0; i < m; i++) {
        while (k > 0 && pattern[k] != pattern[i]) k = lps[k - 1];
        if (pattern[k] == pattern[i]) k++;
        lps[i] = k;
    }
    vector<int> ans;
    for (int i = 0, k = 0; i < text.length(); i++) {
        while (k > 0 && (k == m || pattern[k] != text[i])) k = lps[k - 1];
        if (pattern[k] == text[i]) k++;
        ans.push_back(k);
    }
    return ans;
}

vector<int> z_function(string s) {
    int n = s.length();
    vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++z[i];
        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    }
    return z;
}

string manacher(string s) {
    string ss = "#";
    for (char c : s) { ss += c; ss += "#"; }
    int n = ss.length();
    vector<int> hlen(n, 0);
    int center = 0, right = 0, maxLen = 0, centerIndex = 0;
    for (int i = 0; i < n; i++) {
        if (i < right) hlen[i] = min(right - i, hlen[2 * center - i]);
        while (i - 1 - hlen[i] >= 0 && i + 1 + hlen[i] < n && ss[i - 1 - hlen[i]] == ss[i + 1 + hlen[i]]) {
            hlen[i]++;
        }
        if (i + hlen[i] > right) {
            center = i;
            right = i + hlen[i];
        }
        if (hlen[i] > maxLen) {
            maxLen = hlen[i];
            centerIndex = i;
        }
    }
    return s.substr((centerIndex - maxLen) / 2, maxLen);
}

/* =========================================
   DATA STRUCTURES (Trees & Tries)
========================================= */
class TrieNode {
public:
    bool is_word = false;
    unordered_map<char, TrieNode*> children;
};

class Trie {
    TrieNode* root;
public:
    Trie() { root = new TrieNode(); }
    void insert(string word) {
        TrieNode* node = root;
        for (char c : word) {
            if (!node->children.count(c)) node->children[c] = new TrieNode();
            node = node->children[c];
        }
        node->is_word = true;
    }
    bool search(string word) {
        TrieNode* node = root;
        for (char c : word) {
            if (!node->children.count(c)) return false;
            node = node->children[c];
        }
        return node->is_word;
    }
};

class XORTrieNode {
public:
    int count = 0;
    XORTrieNode* child[2] = {nullptr, nullptr};
};

class XORTrie {
    XORTrieNode* root;
public:
    XORTrie() { root = new XORTrieNode(); }
    void add(int num) {
        XORTrieNode* node = root;
        for (int i = 16; i >= 0; i--) {
            int val = (num >> i) & 1;
            if (!node->child[val]) node->child[val] = new XORTrieNode();
            node = node->child[val];
            node->count++;
        }
    }
};

class UnionFind {
    vector<int> parent, rank;
public:
    int count;
    UnionFind(int n) : count(n), parent(n), rank(n, 1) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    int find(int p) {
        if (p != parent[p]) parent[p] = find(parent[p]);
        return parent[p];
    }
    bool unite(int p, int q) {
        int prt = find(p), qrt = find(q);
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

/* =========================================
   BINARY INDEXED TREES & SEGMENT TREES
========================================= */
class BIT {
    int n;
    vector<ll> tree;
public:
    BIT(int n) : n(n), tree(n + 1, 0) {}
    void update(int i, ll val) {
        for (++i; i <= n; i += i & -i) tree[i] += val;
    }
    ll query(int i) {
        ll res = 0;
        for (++i; i > 0; i -= i & -i) res += tree[i];
        return res;
    }
};

class RBIT {
    int n;
    vector<ll> B1, B2;
    void add(vector<ll>& b, int idx, ll val) {
        for (; idx <= n; idx += idx & -idx) b[idx] += val;
    }
    ll sum(vector<ll>& b, int idx) {
        ll total = 0;
        for (; idx > 0; idx -= idx & -idx) total += b[idx];
        return total;
    }
    ll prefix_sum(int idx) {
        return sum(B1, idx) * idx - sum(B2, idx);
    }
public:
    RBIT(int n) : n(n), B1(n + 1, 0), B2(n + 1, 0) {}
    void range_add(int l, int r, ll x) {
        l++; r++;
        add(B1, l, x); add(B1, r + 1, -x);
        add(B2, l, x * (l - 1)); add(B2, r + 1, -x * r);
    }
    ll range_query(int l, int r) {
        return prefix_sum(r + 1) - prefix_sum(l);
    }
};

class LazySegTree {
    int n;
    vector<ll> tree, lazy;
    void push(int idx, int start, int end) {
        if (lazy[idx] != 0) {
            tree[idx] += (end - start + 1) * lazy[idx];
            if (start != end) {
                lazy[2 * idx + 1] += lazy[idx];
                lazy[2 * idx + 2] += lazy[idx];
            }
            lazy[idx] = 0;
        }
    }
public:
    LazySegTree(int n) : n(n), tree(4 * n, 0), lazy(4 * n, 0) {}
    void update_range(int l, int r, ll val, int idx = 0, int start = 0, int end = -1) {
        if (end == -1) end = n - 1;
        push(idx, start, end);
        if (start > r || end < l) return;
        if (l <= start && end <= r) {
            lazy[idx] += val;
            push(idx, start, end);
            return;
        }
        int mid = start + (end - start) / 2;
        update_range(l, r, val, 2 * idx + 1, start, mid);
        update_range(l, r, val, 2 * idx + 2, mid + 1, end);
        tree[idx] = tree[2 * idx + 1] + tree[2 * idx + 2];
    }
    ll query_range(int l, int r, int idx = 0, int start = 0, int end = -1) {
        if (end == -1) end = n - 1;
        push(idx, start, end);
        if (start > r || end < l) return 0;
        if (l <= start && end <= r) return tree[idx];
        int mid = start + (end - start) / 2;
        return query_range(l, r, 2 * idx + 1, start, mid) + query_range(l, r, 2 * idx + 2, mid + 1, end);
    }
};

class DynamicSegTreeNode {
public:
    DynamicSegTreeNode *left = nullptr, *right = nullptr;
    ll marked = 0;
    bool lazy = false;
};

class DynamicSegTree {
    DynamicSegTreeNode* root;
    ll L, R;
    void push(DynamicSegTreeNode* node, ll l, ll r) {
        if (node->lazy) {
            ll mid = l + (r - l) / 2;
            if (!node->left) node->left = new DynamicSegTreeNode();
            if (!node->right) node->right = new DynamicSegTreeNode();
            node->left->marked = mid - l + 1;
            node->right->marked = r - mid;
            node->left->lazy = true;
            node->right->lazy = true;
            node->lazy = false;
        }
    }
    void add(DynamicSegTreeNode* node, ll l, ll r, ll ql, ll qr) {
        if (ql > r || qr < l) return;
        if (ql <= l && r <= qr) {
            node->marked = r - l + 1;
            node->lazy = true;
            return;
        }
        push(node, l, r);
        ll mid = l + (r - l) / 2;
        if (!node->left) node->left = new DynamicSegTreeNode();
        if (!node->right) node->right = new DynamicSegTreeNode();
        add(node->left, l, mid, ql, qr);
        add(node->right, mid + 1, r, ql, qr);
        node->marked = (node->left ? node->left->marked : 0) + (node->right ? node->right->marked : 0);
    }
public:
    DynamicSegTree(ll L = 0, ll R = 1e9) : L(L), R(R) { root = new DynamicSegTreeNode(); }
    void add(ll l, ll r) { add(root, L, R, l, r); }
    ll query() { return root ? root->marked : 0; }
};

/* =========================================
   GRAPHS & TREES
========================================= */
vector<ll> dijkstra(int n, int src, const vector<vector<pair<int, ll>>>& adj) {
    vector<ll> dist(n, INF);
    dist[src] = 0;
    priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<pair<ll, int>>> pq;
    pq.push({0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        for (auto& edge : adj[u]) {
            int v = edge.first;
            ll w = edge.second;
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist;
}

vector<int> topo_sort(int n, const vector<vector<int>>& adj) {
    vector<int> indegree(n, 0), res;
    for (int u = 0; u < n; u++) {
        for (int v : adj[u]) indegree[v]++;
    }
    queue<int> q;
    for (int i = 0; i < n; i++) {
        if (indegree[i] == 0) q.push(i);
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        res.push_back(u);
        for (int v : adj[u]) {
            if (--indegree[v] == 0) q.push(v);
        }
    }
    return res;
}

class BinaryLifting {
    int n, LOG;
    vector<vector<int>> parent;
    vector<int> depth;
    vector<vector<pair<int, ll>>> adj;

    void dfs(int u, int p) {
        parent[u][0] = p;
        for (auto edge : adj[u]) {
            int v = edge.first;
            if (v != p) {
                depth[v] = depth[u] + 1;
                dfs(v, u);
            }
        }
    }

public:
    BinaryLifting(int n, const vector<vector<pair<int, ll>>>& adj) : n(n), adj(adj) {
        LOG = ceil(log2(n)) + 1;
        parent.assign(n, vector<int>(LOG, -1));
        depth.assign(n, 0);
        dfs(0, -1);
        for (int j = 1; j < LOG; j++) {
            for (int i = 0; i < n; i++) {
                if (parent[i][j - 1] != -1) {
                    parent[i][j] = parent[parent[i][j - 1]][j - 1];
                }
            }
        }
    }

    int lca(int u, int v) {
        if (depth[u] < depth[v]) swap(u, v);
        for (int i = LOG - 1; i >= 0; i--) {
            if (parent[u][i] != -1 && depth[parent[u][i]] >= depth[v]) {
                u = parent[u][i];
            }
        }
        if (u == v) return u;
        for (int i = LOG - 1; i >= 0; i--) {
            if (parent[u][i] != parent[v][i]) {
                u = parent[u][i];
                v = parent[v][i];
            }
        }
        return parent[u][0];
    }
};

/* =========================================
   GEOMETRY
========================================= */
ll cross_product(pair<ll, ll> o, pair<ll, ll> a, pair<ll, ll> b) {
    return (a.first - o.first) * (b.second - o.second) - (a.second - o.second) * (b.first - o.first);
}

vector<pair<ll, ll>> convex_hull(vector<pair<ll, ll>>& points) {
    int n = points.size(), k = 0;
    if (n <= 2) return points;
    vector<pair<ll, ll>> hull(2 * n);
    sort(points.begin(), points.end());
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross_product(hull[k - 2], hull[k - 1], points[i]) <= 0) k--;
        hull[k++] = points[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross_product(hull[k - 2], hull[k - 1], points[i]) <= 0) k--;
        hull[k++] = points[i];
    }
    hull.resize(k - 1);
    return hull;
}
