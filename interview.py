fmin = lambda x, y: x if x < y else y
fmax = lambda x, y: x if x > y else y

directions = [(0, 1), (1, 0), (-1, 0), (0, -1)]

for dx, dy in directions:
    nx, ny = i + dx, j + dy

    if not (0 <= nx < M) or not (0 <= ny < N):
	continue

def rabin_karp(pattern, text):
    """Rabin Karp string search algorithm (Monte Carlo)"""
    R, Q = 128, 997
    hp = ht = 0 
    rm = 1/R
    for p, t in zip(pattern, text):
        hp = (hp * R + ord(p)) % Q
        ht = (ht * R + ord(t)) % Q 
        rm = rm * R % Q 
    if hp == ht: return True 

    m, n = len(pattern), len(text)
    for i in range(m, n):
        ht = (ht + Q - ord(text[i-m])) % Q
        ht = (ht + ord(text[i])) % Q 
        if ht == hp: return True 

    return False 

class RabinKarp: 

    def __init__(self, s): 
        """Calculate rolling hash of s"""
        self.m = 1_111_111_111_111_111_111
        self.pow = [1]
        self.roll = [0] # rolling hash 
        self.N = len(s)

        p = 1_000_000_007
        for x in s: 
            self.pow.append(self.pow[-1] * p % self.m)
            self.roll.append((self.roll[-1] * p + x) % self.m)

    def query(self, i, j): 
        """Return rolling hash of s[i:j]"""
        return (self.roll[j] - self.roll[i] * self.pow[j-i]) % self.m

def kmp(pattern, text):
    '''
    Given a word pattern and target text, returns an array ans
    where ans[i] means the length of the longest prefix of the pattern that is a suffix of text[:i+1].
    In other words, ans[i] is the length of the longest prefix of the pattern that matches a substring ending at position i in the text.
    '''
    k = 0
    lps = [0] 
    for i in range(1, len(pattern)):
	while k and pattern[k] != pattern[i]: k = lps[k-1]
	if pattern[k] == pattern[i]: k += 1
	lps.append(k)
    k = 0
    ans = []
    for i, ch in enumerate(text): 
	while k and (k == len(pattern) or pattern[k] != ch): k = lps[k-1]
	if pattern[k] == ch: k += 1
	ans.append(k)
    return ans
        
class BIT:
    def __init__(self, n):
        self.n = n
        self.tree = [0] * (n + 1)

    def query(self, i):
        ans = 0
        i += 1
        while i > 0:
            ans += self.tree[i]
            i -= (i & (-i))
        return ans

    def update(self, i, value):
        i += 1
        while i <= self.n:
            self.tree[i] += value
            i += (i & (-i))

# Range Binary Indexed Tree
class RBIT:
    def __init__(self, N: int):
        self.b1 = [0] * (N + 1)
        self.b2 = [0] * (N + 1)

    def incr(self, arr, i, x):
        while i <= len(arr):
            arr[i] += x
            i += i & -i

    def range_add(self, l, r, x):
        self.incr(self.b1, l, x)
        self.incr(self.b1, r + 1, -x)
        self.incr(self.b2, l, x * (l - 1))
        self.incr(self.b2, r + 1, -x * r)

    def total(self, arr, i):
        s = 0
        while i > 0:
            s += arr[i]
            i -= i & -i
        return s

    def rtotal(self, i):
        return self.total(self.b1, i) * i - self.total(self.b2, i)

	

class MinSegTree: 
    """A segment tree, aka a statistic tree, is a tree data structure used for 
    storing information about intervals. It allows querying which of the stored 
    segments contain a given point."""

    def __init__(self, arr: List[int]): 
        """Build the setmentation tree."""
        self.n = len(arr)
        self.tree = [0]*(4*self.n)
        self._build(arr, 0, 0, self.n)

    def _build(self, arr: List[int], k: int, lo: int, hi: int) -> None: 
        """Build segment tree from array."""
        if lo+1 == hi: 
            self.tree[k] = arr[lo]
            return 
        mid = lo + hi >> 1
        self._build(arr, 2*k+1, lo, mid)
        self._build(arr, 2*k+2, mid, hi)
        self.tree[k] = min(self.tree[2*k+1], self.tree[2*k+2])

    def update(self, idx: int, val: int, k: int = 0, lo: int = 0, hi: int = 0) -> None:
        """Update segment tree when an array value is changed."""
        if not hi: hi = self.n
        if lo+1 == hi: 
            self.tree[k] = val 
            return 
        mid = lo + hi >> 1
        if idx < mid: self.update(idx, val, 2*k+1, lo, mid) 
        else: self.update(idx, val, 2*k+2, mid, hi)
        self.tree[k] = min(self.tree[2*k+1], self.tree[2*k+2])

    def query(self, qlo: int, qhi: int, k: int = 0, lo: int = 0, hi: int = 0) -> int: 
        """Query value from qlo (inclusive) and qhi (exclusive)."""
        if not hi: hi = self.n
        if qlo <= lo and hi <= qhi: return self.tree[k] # total overlap 
        if qhi <= lo or  hi <= qlo: return inf # no overlap 
        mid = lo + hi >> 1 # partial overlap 
        return min(self.query(qlo, qhi, 2*k+1, lo, mid), self.query(qlo, qhi, 2*k+2, mid, hi))



INF = 10 ** 30
class SegTree:
    def __init__(self, nums):

        self.ar = [0] * (4 * len(nums) + 1)
        self.build(0, 0, len(nums) - 1, nums)
    

    def build(self, i, l, r, nums):
        if l == r:
           
            self.ar[i] = nums[l]
        else:
            m = (l+r) // 2
            self.build(2*i+1, l, m, nums)
            self.build(2*i+2, m+1, r, nums)
            self.ar[i] = self.ar[2*i+1] + self.ar[2*i+2]
        
    def update_v(self, i, l, r, idx, v):
        # self.ST.update_v(0, 0, self.N - 1, index, val)
        if l == r and l == idx:
            self.ar[i] = v

        else:
            m = (l+r)//2
            if m < idx:
                self.update_v(2*i+2, m+1, r, idx, v)
            else:
                self.update_v(2*i+1, l, m, idx, v)

            self.ar[i] = self.ar[2*i+1] + self.ar[2*i+2]


    def query(self, qlo, qhi, i, l, r) -> int:    
        # self.ST.query(left, right, 0, 0, self.N - 1)
        if qlo <= l and r <= qhi:  
            return self.ar[i]
        if qhi < l or r < qlo:  
            return 0
        m = (l + r) // 2
        return self.query(qlo, qhi, 2 * i + 1, l, m) + self.query(qlo, qhi, 2 * i + 2, m + 1, r)


class RangeSegmentTree:
    def __init__(self, size):
        self.n = size
        self.tree = [0] * (4 * size)
        self.lazy = [0] * (4 * size)

    def build(self, start, end, idx, data):
        if start == end:
            self.tree[idx] = data[start]
        else:
            mid = (start + end) // 2
            self.build(start, mid, 2 * idx + 1, data)
            self.build(mid + 1, end, 2 * idx + 2, data)
            self.tree[idx] = self.tree[2 * idx + 1] + self.tree[2 * idx + 2]

    def update_range(self, l, r, val, idx=0, start=0, end=None):
        # ST.update_range(left_bound, right_bound, val)
        # ST.update_range(0, 10, 1)
        if end is None:
            end = self.n - 1

        if self.lazy[idx] != 0:
            self.tree[idx] += (end - start + 1) * self.lazy[idx]
            if start != end: 
                self.lazy[2 * idx + 1] += self.lazy[idx]
                self.lazy[2 * idx + 2] += self.lazy[idx]
            self.lazy[idx] = 0

        if start > r or end < l:
            return

        if l <= start and end <= r:
            self.tree[idx] += (end - start + 1) * val
            if start != end:
                self.lazy[2 * idx + 1] += val
                self.lazy[2 * idx + 2] += val
            return

        mid = (start + end) // 2
        self.update_range(l, r, val, 2 * idx + 1, start, mid)
        self.update_range(l, r, val, 2 * idx + 2, mid + 1, end)
        self.tree[idx] = self.tree[2 * idx + 1] + self.tree[2 * idx + 2]

    def query_range(self, l, r, idx=0, start=0, end=None):
        # ST.query_range(left_bound, right_bound)
        # ST.query_range(0, 10)
        if end is None:
            end = self.n - 1

        if self.lazy[idx] != 0:
            self.tree[idx] += (end - start + 1) * self.lazy[idx]
            if start != end:
                self.lazy[2 * idx + 1] += self.lazy[idx]
                self.lazy[2 * idx + 2] += self.lazy[idx]
            self.lazy[idx] = 0

        if start > r or end < l:
            return 0

        if l <= start and end <= r:
            return self.tree[idx]

        mid = (start + end) // 2
        left_query = self.query_range(l, r, 2 * idx + 1, start, mid)
        right_query = self.query_range(l, r, 2 * idx + 2, mid + 1, end)
        return left_query + right_query


'''
segment tree for flipping 1s to 0s and 0s to 1s
'''
class RangeSegmentTree2:
    def __init__(self, size):
        self.n = size
        self.tree = [0] * (4 * size)
        self.lazy = [0] * (4 * size)

    def build(self, start, end, idx, data):
        '''
        0, N - 1, 0, nums
        '''
        if start == end:
            self.tree[idx] = data[start]
        else:
            mid = (start + end) // 2
            self.build(start, mid, 2 * idx + 1, data)
            self.build(mid + 1, end, 2 * idx + 2, data)
            self.tree[idx] = self.tree[2 * idx + 1] + self.tree[2 * idx + 2]

    def propagate(self, start, end, idx):
        if self.lazy[idx] == 1:  # Flip bits if lazy flag is set
            self.tree[idx] = (end - start + 1) - self.tree[idx]  # Flip values in the range
            if start != end:  # Propagate to children
                self.lazy[2 * idx + 1] ^= 1  # Toggle lazy flags for children
                self.lazy[2 * idx + 2] ^= 1
            self.lazy[idx] = 0  # Clear the lazy flag for the current node

    def update_range(self, l, r, idx=0, start=0, end=None):
        if end is None:
            end = self.n - 1

        self.propagate(start, end, idx)  # Ensure the current node is up-to-date

        if start > r or end < l:
            return  # No overlap

        if l <= start and end <= r:
            # Fully within range: flip this range
            self.tree[idx] = (end - start + 1) - self.tree[idx]
            if start != end:  # Mark children for lazy propagation
                self.lazy[2 * idx + 1] ^= 1
                self.lazy[2 * idx + 2] ^= 1
            return

        # Partial overlap
        mid = (start + end) // 2
        self.update_range(l, r, 2 * idx + 1, start, mid)
        self.update_range(l, r, 2 * idx + 2, mid + 1, end)
        self.tree[idx] = self.tree[2 * idx + 1] + self.tree[2 * idx + 2]

    def query_range(self, l, r, idx=0, start=0, end=None):
        if end is None:
            end = self.n - 1

        self.propagate(start, end, idx)  # Ensure the current node is up-to-date

        if start > r or end < l:
            return 0  # No overlap

        if l <= start and end <= r:
            return self.tree[idx]  # Fully within range

        # Partial overlap
        mid = (start + end) // 2
        left_query = self.query_range(l, r, 2 * idx + 1, start, mid)
        right_query = self.query_range(l, r, 2 * idx + 2, mid + 1, end)
        return left_query + right_query



class RangeZeroSegmentTree:
    def __init__(self, size):
        self.size = size
        self.tree = [{'cnt_add': 0, 'total': 0} for _ in range(4 * size)]
        self.init()

    def init(self, x=1, l=0, r=None):
        if r is None:
            r = self.size

        self.tree[x] = {'cnt_add': 0, 'total': 0}
        if r - l < 2:
            return

        m = (l + r) // 2
        self.init(2 * x, l, m)
        self.init(2 * x + 1, m, r)

    def update_range(self, L, R, v, x=1, l=0, r=None):
        # ST.update_range(0, 9, 1)  # update range [0, 9) with val
        if r is None:
            r = self.size

        if r <= L or R <= l:
            return
        if L <= l and r <= R:
            self.tree[x]['cnt_add'] += v
            if self.tree[x]['cnt_add']:
                self.tree[x]['total'] = r - l
            else:
                if r - l > 1:
                    self.tree[x]['total'] = self.tree[2 * x]['total'] + self.tree[2 * x + 1]['total']
                else:
                    self.tree[x]['total'] = 0
            return

        m = (l + r) // 2
        self.update_range(L, R, v, 2 * x, l, m)
        self.update_range(L, R, v, 2 * x + 1, m, r)

        if self.tree[x]['cnt_add']:
            self.tree[x]['total'] = r - l
        else:
            self.tree[x]['total'] = self.tree[2 * x]['total'] + self.tree[2 * x + 1]['total']

    def get_zeroes(self):
        # returns number of 0s in the array
        return self.size - self.tree[1]['total']



INF = 10 ** 30
class SpecialSegTree:
    def __init__(self, nums):

        self.ar = [[None] * 4 for _ in range(4 * len(nums) + 1)]
        self.build(0, 0, len(nums) - 1, nums)
    
    def combine(self, t1, t2):
        nn1, ny1, yn1, yy1 = t1
        nn2, ny2, yn2, yy2 = t2
        nn = max(nn1 + yn2, nn1 + nn2, ny1 + nn2)
        ny = max(nn1 + yy2, nn1 + ny2, ny1 + ny2)
        yn = max(yn1 + yn2, yn1 + nn2, yy1 + nn2)
        yy = max(yn1 + yy2, yn1 + ny2, yy1 + ny2)
        return [nn, ny, yn, yy]
    
    def build(self, i, l, r, nums):
        
        if l == r:
            # One element, YN and NY doesnt exist
            self.ar[i] = [0, -INF, -INF, nums[l]]
        else:
            m = (l+r) // 2
            self.build(2*i+1, l, m, nums)
            self.build(2*i+2, m+1, r, nums)
            self.ar[i] = self.combine(self.ar[2*i+1], self.ar[2*i+2])
        
    def update_v(self, i, l, r, idx, v):
        if l == r and l == idx:

            self.ar[i] = [0, -INF, -INF, v]

        else:
            m = (l+r)//2
            if m < idx:
                self.update_v(2*i+2, m+1, r, idx, v)
            else:
                self.update_v(2*i+1, l, m, idx, v)
            self.ar[i] = self.combine(self.ar[2*i+1], self.ar[2*i+2])




def manacher(s: str) -> str:               
    """Return longest palindromic substring via Manacher's algo."""
    ss = "#" + "#".join(s) + "#" # augmented string (even-length palindromes)
    n = len(ss)
    hlen = [0] * n # half-length
    center = right = 0
    for i in range(n):
        if i < right: hlen[i] = min(right-i, hlen[2*center-i])
        while 0 <= i-1-hlen[i] and i+1+hlen[i] < len(ss) and ss[i-1-hlen[i]] == ss[i+1+hlen[i]]: hlen[i] += 1
        if right < i+hlen[i]: center, right = i, i+hlen[i]
    xx, ii = max((x, i) for i, x in enumerate(hlen))
    return s[(ii-xx)//2 : (ii+xx)//2]


class CompressedSegmentTree:
    def __init__(self, xs):
        xs = list(sorted(set(xs)))
        self.xs = xs
        xl = {}
        for i, x in enumerate(xs):
            xl[x] = i

        self.xl = xl
        self.N = len(xl)
        # how many segments cover this node
        self.count = [0] * (4 * self.N + 1)
        # how much of this node is covered
        self.total = [0] * (4 * self.N + 1)

    def _update_range(self, node, l, r, ql, qr, val):
        if ql >= qr:
            return
        
        if ql == self.xs[l] and qr == self.xs[r]:
            self.count[node] += val
        else:
            mid = (l + r) // 2
            if ql < self.xs[mid]:
                self._update_range(2 * node, l, mid, ql, min(qr, self.xs[mid]), val)
            if qr > self.xs[mid]:
                self._update_range(2 * node + 1, mid, r, max(ql, self.xs[mid]), qr, val)

        if self.count[node] > 0:
            self.total[node] = self.xs[r] - self.xs[l]
        else:
            self.total[node] = self.total[node * 2] + self.total[node * 2 + 1] if node * 2 + 1 < len(self.total) else 0
            
    def update_range(self, ql, qr, val):
        """ Adds 'val' to all indices in range [ql, qr]. """
        self._update_range(1, 0, self.N - 1, ql, qr, val)

    # get total number of marked nodes
    def query(self):
        return self.total[1]


# marks vals as 1 (up to 10 ** 9)
class Node:
    def __init__(self):
        self.left = None
        self.right = None
        self.marked = 0  # Count of marked spots in this range
        self.lazy = 0     # Lazy propagation flag

class SegmentTree:
    def __init__(self, L=0, R=10**9):
        self.root = Node()
        self.L = L
        self.R = R

    def _push(self, node, l, r):
        """ Propagate the lazy update if needed. """
        if node.lazy:
            mid = (l + r) // 2
            if not node.left:
                node.left = Node()
            if not node.right:
                node.right = Node()
            node.left.marked = mid - l + 1
            node.right.marked = r - mid
            node.left.lazy = 1
            node.right.lazy = 1
            node.lazy = 0  # Clear lazy flag

    def _add(self, node, l, r, ql, qr):
        """ Marks the range [ql, qr] as 1. """
        if ql > r or qr < l:
            return
        if ql <= l and r <= qr:
            node.marked = r - l + 1
            node.lazy = 1
            return
        self._push(node, l, r)
        mid = (l + r) // 2
        if not node.left:
            node.left = Node()
        if not node.right:
            node.right = Node()
        self._add(node.left, l, mid, ql, qr)
        self._add(node.right, mid + 1, r, ql, qr)
        node.marked = (node.left.marked if node.left else 0) + (node.right.marked if node.right else 0)

    def add(self, l, r):
        self._add(self.root, self.L, self.R, l, r)

    def query(self):
        return self.root.marked if self.root else 0



def z_function(s):
    z, l, r, n = [0] * len(s), 0, 0, len(s)
    for i in range(1, n):
	if i < r:
	    z[i] = min(r - i, z[i - l])
	while i + z[i] < n and s[i + z[i]] == s[z[i]]:
	    z[i] += 1
	if i + z[i] > r:
	    l, r = i, i + z[i]

    return z

def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

class UnionFind:

    def __init__(self, N):
        self.count = N              
        self.parent = [i for i in range(N)]
        self.rank = [1] * N
        
        
    def find(self, p):
        if p != self.parent[p]:
            self.parent[p] = self.find(self.parent[p]) 
        return self.parent[p]

    def union(self, p, q):
        prt, qrt = self.find(p), self.find(q)
        if prt == qrt: return False
        if self.rank[prt] >= self.rank[qrt]: 
            self.parent[qrt] = prt
            self.rank[prt] += self.rank[qrt]
        else:
            self.parent[prt] = qrt
            self.rank[qrt] += self.rank[prt]
            
        self.count -= 1 
        return True 
	    

class BinaryLifting:
    def __init__(self, n, edges):
        self.n = n
        self.LOG = n.bit_length()
        self.C = 27
        self.graph = defaultdict(list)
        self.parent = [[-1] * n for _ in range(self.LOG)]
        self.depth = [0] * n
        self.weight_count = [[0] * self.C for _ in range(n)]
        self._build_graph(edges)
        self._preprocess()

    def _build_graph(self, edges):
        for u, v, w in edges:
            self.graph[u].append((v, w))
            self.graph[v].append((u, w))

    def _dfs(self, node, par):
        self.parent[0][node] = par
        for neighbor, weight in self.graph[node]:
            if neighbor == par:
                continue
            self.depth[neighbor] = self.depth[node] + 1
            self.weight_count[neighbor] = self.weight_count[node].copy()
            self.weight_count[neighbor][weight] += 1
            self._dfs(neighbor, node)

    def _preprocess(self):
        self._dfs(0, -1)
        for i in range(1, self.LOG):
            for j in range(self.n):
                if self.parent[i - 1][j] != -1:
                    self.parent[i][j] = self.parent[i - 1][self.parent[i - 1][j]]

    def get_kth_ancestor(self, node, k):
        for i in range(self.LOG):
            if k & (1 << i):
                node = self.parent[i][node]
                if node == -1:
                    break
        return node

    def lca(self, u, v):
        if self.depth[u] < self.depth[v]:
            u, v = v, u
        diff = self.depth[u] - self.depth[v]
        u = self.get_kth_ancestor(u, diff)
        if u == v:
            return u
        for i in range(self.LOG - 1, -1, -1):
            if self.parent[i][u] != self.parent[i][v]:
                u = self.parent[i][u]
                v = self.parent[i][v]
        return self.parent[0][u]

    def min_operations_queries(self, queries: List[List[int]]) -> List[int]:
        res = []
        for x, y in queries:
            l = self.lca(x, y)
            length = self.depth[x] + self.depth[y] - 2 * self.depth[l]
            max_z = max(self.weight_count[x][z] + self.weight_count[y][z] - 2 * self.weight_count[l][z] for z in range(self.C))
            res.append(length - max_z)
        return res






def merge_sort(arr):
    if len(arr) <= 1:
        return arr

    mid = len(arr) // 2
    left = merge_sort(arr[:mid])
    right = merge_sort(arr[mid:])

    return merge(left, right)

def merge(left, right):
    merged = []
    l, r = 0, 0

    while l < len(left) and r < len(right):
        if left[l] <= right[r]:
            merged.append(left[l])
            l += 1
        else:
            merged.append(right[r])
            r += 1

    merged.extend(left[l:])
    merged.extend(right[r:])

    return merged


class MAX_BIT:
    def __init__(self, N):
        self.INF = 10 ** 20
        self.stree = [-self.INF] * (N + 1)
    
    def update(self, i, x):
	i += 1
        while i < len(self.stree):
            self.stree[i] = max(self.stree[i], x)
            i |= (i + 1)

    def query_max(self, i):
        s = -self.INF
	i += 1

        while i >= 0:
            s = max(s, self.stree[i])
            i &= (i + 1)
            i -= 1

        return s


# Floyd-Warshall
INF = 10 ** 20
e = [[INF] * N for _ in range(N)]

for i in range(N):
    e[i][i] = 0

for u, v, w in edges:
    e[u][v] = min(e[u][v], w)
    e[v][u] = min(e[v][u], w)

for k in range(N):
	for i in range(N):
	    for j in range(N):
		e[i][j] = min(e[i][j], e[i][k] + e[k][j])






def lca(node):

    if not node:
	return None

    if node in nodes: 
	return node

    left = lca(node.left)

    right = lca(node.right)

    if left and right:
	return node

    if left:
	return left

    return right


class TrieNode:
    def __init__(self):
        self.is_word = False
        self.children = {}

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        node = self.root
        for c in word:
            if c not in node.children:
                node.children[c] = TrieNode()
            node = node.children[c]
        node.is_word = True

    def search(self, word: str) -> bool:
        node = self.root
        for c in word:
            if c not in node.children:
                return False
            node = node.children[c]
        return node.is_word

    def starts_with(self, prefix: str) -> bool:
        node = self.root
        for c in prefix:
            if c not in node.children:
                return False
            node = node.children[c]
        return True



class BIT:
    # ranged query and ranged update 
    def __init__(self, size):
        self.n = size
        self.B1 = [0] * (self.n+1)
        self.B2 = [0] * (self.n+1)

    def incr(self, s, e):
        # [s, e]
        self.rangeUpdate(s, e, 1)

    def rangeUpdate(self, s, e, d):
        s += 1
        e += 1
        self._add(self.B1, s, d)
        self._add(self.B1, e+1, -d)
        self._add(self.B2, s, d*(s-1))
        self._add(self.B2, e+1, -d*e)

    def _add(self, B, k, d):
        while k <= self.n:
            B[k] += d
            k += k & -k

    def rangeQuery(self, s, e):
        # [s, e]
        s, e = s+1, e+1
        return self._sum(e) - self._sum(s-1)

    def _sum(self, k):
        return self._sumOf(self.B1, k) * k - self._sumOf(self.B2, k)

    def _sumOf(self, B, k):
        ans = 0
        while k > 0:
            ans += B[k]
            k -= k & -k
        return ans

B = BIT(N + 5)
B.rangeQuery(l, r)
B.incr(l, r)


def get_divisors(n):
    divisors = set()
    for i in range(1, ceil(n**0.5) + 1):
        if n % i == 0:
            divisors.add(i)
            divisors.add(n // i)
            
    return divisors


# primes
prime = [True] * N
prime[0] = prime[1] = False
for i in range(2, N):
    if prime[i]:
	x = i + i
	while x < N:
	    prime[x] = False
	    x += i


# Rotate 90 degrees
matrix = list(zip(*reversed(matrix)))

# Transpose
matrix = list(zip(*matrix))



# tree diameter and tree center
def bfs_farthest_node(start, graph):
    n = len(graph)
    visited = [False] * n
    dist = [0] * n
    queue = deque([start])
    visited[start] = True
    farthest_node = start

    while queue:
        node = queue.popleft()
        for neighbor in graph[node]:
            if not visited[neighbor]:
                visited[neighbor] = True
                dist[neighbor] = dist[node] + 1
                queue.append(neighbor)
                if dist[neighbor] > dist[farthest_node]:
                    farthest_node = neighbor

    return farthest_node, dist[farthest_node]

def find_tree_diameter(edges):
    n = len(edges) + 1
    # Build the adjacency list
    graph = defaultdict(list)
    for u, v in edges:
        graph[u].append(v)
        graph[v].append(u)

    # Start from any node (Node 0 in this case)
    farthest_node, _ = bfs_farthest_node(0, graph)
    # Find the farthest node from the previously found farthest node
    farthest_node, diameter = bfs_farthest_node(farthest_node, graph)
    return diameter

def find_tree_center(edges):
    n = len(edges) + 1
    if n == 1:
        return [0]

    # Build the adjacency list
    graph = defaultdict(list)
    degree = [0] * n
    for u, v in edges:
        graph[u].append(v)
        graph[v].append(u)
        degree[u] += 1
        degree[v] += 1

    # Initialize leaves
    leaves = deque()
    for i in range(n):
        if degree[i] == 1:
            leaves.append(i)

    # Trim leaves until reaching the center
    remaining_nodes = n
    while remaining_nodes > 2:
        leaves_count = len(leaves)
        remaining_nodes -= leaves_count
        for _ in range(leaves_count):
            leaf = leaves.popleft()
            for neighbor in graph[leaf]:
                degree[neighbor] -= 1
                if degree[neighbor] == 1:
                    leaves.append(neighbor)

    return list(leaves)



def convex_hull(points):
    def clockwise(p1,p2,p3):
	x1,y1=p1
	x2,y2=p2
	x3,y3=p3
	
	return ((y3-y2)*(x2-x1)-(y2-y1)*(x3-x2))
    points.sort()
    upper=[]
    lower=[]
    for t in points:
	while len(upper)>1 and clockwise(upper[-2],upper[-1],t)>0:
	    upper.pop()
	while len(lower)>1 and clockwise(lower[-2],lower[-1],t)<0:
	    lower.pop()
	upper.append(tuple(t))
	lower.append(tuple(t))
	
    return list(set(upper+lower))   





class XORTrieNode:
    def __init__(self):
        self.count = 0
        self.children = {}

class XORTrie:
    def __init__(self):
        self.root = TrieNode()

    def add(self, num: int) -> None:
        node = self.root
        for i in range(16, -1, -1):
            val = (num >> i) & 1
            if val not in node.children:
                node.children[val] = TrieNode()
            node = node.children[val]
            node.count += 1



    def count(self, num: int, upper: int) -> int:
        node = self.root
        total = 0

        for i in range(16, -1, -1):
            val = (num >> i) & 1
            val2 = (upper >> i) & 1

            if val == 1 and val2 == 1:
                if 1 in node.children:
                    total += node.children[1].count
            
                if 0 not in node.children:
                    return total
                
                node = node.children[0]

                
            elif val == 1 and val2 == 0:
                if 1 not in node.children:
                    return total
                
                node = node.children[1]

            elif val == 0 and val2 == 1:
                if 0 in node.children:
                    total += node.children[0].count
            
                if 1 not in node.children:
                    return total
                
                node = node.children[1]

            
            elif val == 0 and val2 == 0:
                if 0 not in node.children:
                    return total
                node = node.children[0]
                
                
        total += node.count
        return total





def shortest_path_calc(vertex_weights, lookup):
    INF = 10 ** 20
    min_dist = [INF] * N
    min_dist[0] = vertex_weights[0]

    h = []
    heapq.heappush(h, (vertex_weights[0], 0))

    while h:
        val, u = heapq.heappop(h)

        if min_dist[u] < val:
            continue

        for v, w in lookup[u]:
            vertex_weight = vertex_weights[v]

            total_weight = w + val + vertex_weight

            if total_weight >= min_dist[v]:
                continue

            min_dist[v] = total_weight
            heapq.heappush(h, (total_weight, v))

    return min_dist

  def shortest_path_calc(src):
    min_dist = [INF] * N
    min_dist[src] = 0

    h = []
    heapq.heappush(h, (0, src))

    while h:
	val, u = heapq.heappop(h)

	if min_dist[u] < val:
	    continue

	for v, w in lookup[u]:
	    
	    total_weight = w + val 

	    if total_weight >= min_dist[v]:
		continue

	    min_dist[v] = total_weight
	    heapq.heappush(h, (total_weight, v))

    return min_dist



  def lca(node):
    if not node:
	return None

    if node in [p, q]:
	return node

    left = lca(node.left)
    right = lca(node.right)

    if left and right:
	return node

    if left:
	return left

    return right


# topological sort
edges_to = defaultdict(set)
edges_from = defaultdict(set)

for u, v in edges:
    edges_to[v].add(u)
    edges_from[u].add(v)

def topo_sort(edges_to, edges_from):
    q = deque()
    for v in range(1, k + 1):
	if len(edges_from[v]) == 0:
	    q.append(v)

    arr = []

    while q:
	u = q.popleft()
	arr.append(u)
	for v in edges_to[u]:
	    edges_from[v].remove(u)
	    if not edges_from[v]:
		q.append(v)

    arr.reverse()
    return arr


MOD = 10 ** 9 + 7
MX = N + 30
 
fac = [1] * MX
for i in range(1, MX):
    fac[i] = fac[i-1] * i % MOD
ifac = [pow(fac[MX - 1], MOD-2, MOD)] * MX
for i in range(MX - 1, 0, -1):
    ifac[i-1] = ifac[i] * i % MOD
 
def n_choose_k(N,K):
    return (fac[N]*ifac[N-K]*ifac[K])%MOD











#start SortedList
from bisect import bisect_left
from bisect import bisect_right
from functools import reduce

class FenwickTree:
    def __init__(self, x):
        bit = self.bit = list(x)
        size = self.size = len(bit)
        for i in range(size):
            j = i | (i + 1)
            if j < size:
                bit[j] += bit[i]

    def update(self, idx, x):
        """updates bit[idx] += x"""
        while idx < self.size:
            self.bit[idx] += x
            idx |= idx + 1

    def __call__(self, end):
        """calc sum(bit[:end])"""
        x = 0
        while end:
            x += self.bit[end - 1]
            end &= end - 1
        return x

    def find_kth(self, k):
        """Find largest idx such that sum(bit[:idx]) <= k"""
        idx = -1
        for d in reversed(range(self.size.bit_length())):
            right_idx = idx + (1 << d)
            if right_idx < self.size and self.bit[right_idx] <= k:
                idx = right_idx
                k -= self.bit[idx]
        return idx + 1, k


class SortedList:
    block_size = 700

    def __init__(self, iterable=()):
        self.macro = []
        self.micros = [[]]
        self.micro_size = [0]
        self.fenwick = FenwickTree([0])
        self.size = 0
        for item in iterable:
            self.add(item)

    def add(self, x):
        i = bisect_left(self.macro, x)
        j = bisect_right(self.micros[i], x)
        self.micros[i].insert(j, x)
        self.size += 1
        self.micro_size[i] += 1
        self.fenwick.update(i, 1)
        if len(self.micros[i]) >= self.block_size:
            self.micros[i:i + 1] = self.micros[i][:self.block_size >> 1], self.micros[i][self.block_size >> 1:]
            self.micro_size[i:i + 1] = self.block_size >> 1, self.block_size >> 1
            self.fenwick = FenwickTree(self.micro_size)
            self.macro.insert(i, self.micros[i + 1][0])

    def remove(self,x):
        i = self.bisect_left(x)
        if i >= self.size:
            return
        if self.__getitem__(i) == x:
            self.pop(i)
        return

    def pop(self, k=-1):
        i, j = self._find_kth(k)
        self.size -= 1
        self.micro_size[i] -= 1
        self.fenwick.update(i, -1)
        return self.micros[i].pop(j)

    def __getitem__(self, k):
        if isinstance(k,slice):
            start, stop, step = k.indices(self.size)
            return [self.__getitem__(x) for x in range(start,stop,step)]
        i, j = self._find_kth(k)
        return self.micros[i][j]

    def count(self, x):
        return self.bisect_right(x) - self.bisect_left(x)

    def __contains__(self, x):
        return self.count(x) > 0

    def bisect_left(self, x):
        i = bisect_left(self.macro, x)
        return self.fenwick(i) + bisect_left(self.micros[i], x)

    def bisect_right(self, x):
        i = bisect_right(self.macro, x)
        return self.fenwick(i) + bisect_right(self.micros[i], x)

    def _find_kth(self, k):
        return self.fenwick.find_kth(k + self.size if k < 0 else k)

    def __len__(self):
        return self.size

    def __iter__(self):
        return (x for micro in self.micros for x in micro)

    def __repr__(self):
        return str(list(self))


# linear matrix exponentiation
M = 10 ** 9 + 7

# perform AB matrix multiplication
def mat_mult(A, B):
    size = len(A)
    result = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                result[i][j] = (result[i][j] + A[i][k] * B[k][j]) % M
    return result


# perform mat ** exp in log exp time
def mat_pow(mat, exp):
    size = len(mat)

    # itentity matrix
    result = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    base = mat
    while exp > 0:
        if exp % 2 == 1:
            # return current base
            result = mat_mult(result, base)

        # otherwise, base = base * base
        base = mat_mult(base, base)
        exp //= 2
    return result


# Precompute factorials and inverse factorials for combinations
fact = [1]*(n+1)
invfact = [1]*(n+1)
for i in range(1, n+1):
    fact[i] = (fact[i-1] * i) % MOD
invfact[n] = pow(fact[n], MOD-2, MOD)
for i in reversed(range(n)):
    invfact[i] = (invfact[i+1]*(i+1)) % MOD

def comb(n, r):
    if r < 0 or r > n:
	return 0
    return (fact[n]*((invfact[r]*invfact[n-r]) % MOD)) % MOD

# Build prefix sums of C(x, 0..r) for r up to k
# combPrefix[x][r] = sum of C(x, m) for m = 0..r
combPrefix = [[0]*(k+1) for _ in range(n+1)]
for x in range(n+1):
    combPrefix[x][0] = 1
    for r in range(1, k+1):
	combPrefix[x][r] = (combPrefix[x][r-1] + comb(x, r)) % MOD


'''
triangle calculation

given an interval [l, r] and a point index, coun the number of subintervals within interval that 1. contain index and 2. are at most k
'''
# triangle sum for interval length at most k
interval_length = min(k, left_df + right_df + 1)
t = interval_length * (interval_length + 1) // 2

# cut smaller triangles
if left_df < interval_length - 1:
    q = interval_length - 1 - left_df
    t -= q * (q + 1) // 2

if right_df < interval_length - 1:
    q = interval_length - 1 - right_df
    t -= q * (q + 1) // 2

total += t * key



# has cycle

# has cycle
def has_cycle(adj, num_nodes):
    indegree = [0] * num_nodes
    for u in range(num_nodes):
        for v in adj[u]:
            indegree[v] += 1

    q = deque(i for i in range(num_nodes) if indegree[i] == 0)
    count = 0

    while q:
        u = q.popleft()
        count += 1
        for v in adj[u]:
            indegree[v] -= 1
            if indegree[v] == 0:
                q.append(v)

    return count != num_nodes
	
NEW, VISITING, DONE = range(3)
def has_cycle(adj):
    color = defaultdict(int)
    def dfs(node):
        if color[node] != NEW:
            return color[node] == VISITING
        color[node] = VISITING
        if any(dfs(nei) for nei in adj[node]):
            return True
        color[node] = DONE
        return False

    return any(dfs(u) for u in list(adj.keys()))






# Adding number of ways to choose all possible paths from group of branches
total = 0
for A in branches.values():
    s = 1
    for x in A:
	total += s * x
	s += x
return total




# string pattern thing -> Lex largest substring
N = len(word)
i = 0
j = 1

while j < N:
    k = 0
    while j + k < N and word[i + k] == word[j + k]:
	k += 1
    if j + k < N and word[i + k] < word[j + k]:
	i, j = j, max(j + 1, i + k + 1)
    else:
	j = j + k + 1

# Hungarian Algorithm
# n^3 time n^2 space
# cost matrix
# rows: strength[i], cols: X = j + 1, 1 <= X <= N
mat =[[0] * N for i in range(N)]
for i in range(N):
    for j in range(N):
	# for num := strength[i], and X = j + 1, calc time and strore in matrix
	mat[i][j]= strength[i] // (j + 1) + int(strength[i] % (j + 1) > 0)

# use hungarian algo to get best is and js pairs
A, B = linear_sum_assignment(mat)



# cool binary lifting example
 A = []
for index, x in enumerate(nums):
    A.append((x, index))
A.sort()
lookup = defaultdict(int)
for index, (v, i) in enumerate(A):
    lookup[i] = index

mx = len(bin(n)[2:]) + 1

binary_lift = [[-1] * mx for _ in range(n)]

r = 0

for l in range(n):
    while r < n and A[r][0] - A[l][0] <= maxDiff:
	r += 1
    r -= 1
    binary_lift[l][0] = r

for b in range(1, mx):
    for i in range(n):
	binary_lift[i][b] = binary_lift[binary_lift[i][b-1]][b - 1]


INF = 10 ** 20
def find(u, v, b):
    if u == v:
	return 0

    if binary_lift[u][b] < v:
	return INF

    if binary_lift[u][0] >= v:
	return 1

    
    for bb in range(b, -1, -1):
	if binary_lift[u][bb] < v:
	    return (1 << bb) + find(binary_lift[u][bb], v, bb)



ans = []

for u, v in queries:
    u, v = lookup[u], lookup[v]
    u, v = min(u, v), max(u, v)
    res = find(u, v, mx - 1)
    if res < INF:
	ans.append(res)
    else:
	ans.append(-1)

return ans



class BitTrieNode:
    def __init__(self):
        self.count = 0
        self.child = [None] * 2
 
class BitTrie:
    def __init__(self):
        self.root = BitTrieNode()

    def increase(self, number, d):
        cur = self.root 
        for i in range(17, -1, -1):
            bit = (number >> i) & 1

            if not cur.child[bit]: 
                cur.child[bit] = BitTrieNode()

            cur = cur.child[bit]
            cur.count += d

    def findMax(self, number):

        cur, ans = self.root, 0
        for i in range(17, -1, -1):

            bit = (number >> i) & 1

            if cur.child[1 - bit] and cur.child[1-bit].count > 0:
                cur = cur.child[1 - bit]
                ans |= (1 << i)

            elif cur.child[bit] and cur.child[bit].count > 0:
                cur = cur.child[bit]

            else:
                return ans

        return ans
