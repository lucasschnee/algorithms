n & n-1 # unset last set bit 
n & -n # retain last set bit 

directions = [(0, 1), (1, 0), (-1, 0), (0, -1)]

for dx, dy in directions:
    nx, ny = x + dx, y + dy

    if not (0 <= nx < M) or not (0 <= ny < N):
	continue

def boyer_moore(pattern, text):
    """Boyer-Moore string search algorithm"""
    m, n = len(pattern), len(text)
    jump = {}
    for i, c in enumerate(pattern):
        jump[c] = i 

    offset = 0
    while offset < n-m:
        skip = 0
        for j in reversed(range(m)):
            if text[offset+j] != pattern[j]:
                skip = max(1, j - jump.get(pattern[j], -1))
                break 
        if skip == 0: return True
        offset += skip
    return False #not found 


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


 def KMP(pref, s):
    m, n = len(pref), len(s)
    lps = [0] * m
    j = 0  # length of the previous longest prefix suffix
    
    # Preprocessing the pattern
    for i in range(1, m):
	while j and pref[j] != pref[i]:
	    j = lps[j - 1]
	if pref[j] == pref[i]:
	    j += 1
	lps[i] = j

    # Searching the pattern in the text
    occurrences = []
    j = 0
    for i in range(n):
	while j and pref[j] != s[i]:
	    j = lps[j - 1]
	if pref[j] == s[i]:
	    j += 1
	if j == m:
	    occurrences.append(i - m + 1)
	    j = lps[j - 1]
    return occurrences
        
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
    if n <= 3:
	return True
    if n % 2 == 0 or n % 3 == 0:
	return False
    i = 5
    while i * i <= n:
	if n % i == 0 or n % (i + 2) == 0:
	    return False
	i += 6
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


# divisors, factors
divisors = set()
for i in range(1, ceil(n**0.5) + 1): 
    if n % i == 0: 
	divisors.add(i)   
	divisors.add(n // i)  


# primes
primes = [0] * n

total = 0

for i in range(2, n):
    
    if not primes[i]:
	total += 1

	for j in range(i, n, i):
	    primes[j] = 1


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

