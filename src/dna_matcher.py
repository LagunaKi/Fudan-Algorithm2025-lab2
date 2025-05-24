import collections
from typing import List, Tuple

class Anchor:
    def __init__(self, q_start, r_start, length, strand):
        self.q_start = q_start  # query起始坐标
        self.r_start = r_start  # ref起始坐标
        self.length = length    # 匹配长度
        self.strand = strand    # 1=正向, -1=反向

    def __repr__(self):
        s = '+' if self.strand == 1 else '-'
        return f"Anchor(q={self.q_start}, r={self.r_start}, len={self.length}, strand={s})"

# 反向互补
def reverse_complement(seq: str) -> str:
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

# minimizer稀疏化提取
def minimize(seq: str, k: int, w: int) -> List[Tuple[str, int]]:
    """
    提取minimizer稀疏化的kmer及其在seq中的起始位置。
    每个窗口w内只保留字典序最小的kmer。
    """
    minimizers = []
    n = len(seq)
    if n < k:
        return minimizers
    for i in range(n - w - k + 2):
        window = [seq[j:j+k] for j in range(i, i+w)]
        min_kmer = min(window) # 取长度为w的窗口内字典序最小的kmer
        min_pos = i + window.index(min_kmer) # 记录位置
        if not minimizers or minimizers[-1][1] != min_pos: # 不重复记录
            minimizers.append((min_kmer, min_pos))
    return minimizers

# 合并连续锚点
'''！可模糊合并'''
def merge_anchors(anchors: List[Anchor]) -> List[Anchor]:
    """
    合并query和ref上都连续的锚点为一个更长的锚点。
    """
    if not anchors:
        return []
    anchors.sort(key=lambda a: (a.q_start, a.r_start, a.strand))
    merged = [anchors[0]]
    for a in anchors[1:]:
        last = merged[-1] # 已合并的最后一个锚点
        # 判断是否可合并：方向相同，且query和ref都连续
        if (a.strand == last.strand and
            a.q_start == last.q_start + last.length and
            a.r_start == last.r_start + last.length):
            last.length += a.length
        else:
            merged.append(a)
    return merged

# 提取正向和反向的minimizer锚点，并合并连续锚点
def build_anchors(query: str, ref: str, k: int, w: int) -> List[Anchor]:
    anchors = []
    # 构建ref的kmer哈希表（kmer->起点位置）
    ref_kmers = collections.defaultdict(list)
    for i in range(len(ref) - k + 1):
        kmer = ref[i:i+k]
        ref_kmers[kmer].append(i)
    # 正向minimize
    for kmer, q_pos in minimize(query, k, w): # 遍历query的minimizer
        for r_pos in ref_kmers.get(kmer, []): # 在ref中查找kmer
            anchors.append(Anchor(q_pos, r_pos, k, 1))
    # 反向minimize
    rc_query = reverse_complement(query)
    for kmer, rc_q_pos in minimize(rc_query, k, w):
        for r_pos in ref_kmers.get(kmer, []):
            q_start = len(query) - (rc_q_pos + k)
            anchors.append(Anchor(q_start, r_pos, k, -1))
    # 合并连续锚点
    return merge_anchors(anchors)

# 用二分查找优化的DP，O(A log A)
def dp(anchors: List[Anchor]) -> List[Anchor]:
    """
    用二分查找优化的DP，保证主链无重叠，支持倒位。
    """
    # 按ref坐标排序，便于二分查找
    anchors.sort(key=lambda a: (a.r_start, a.q_start, a.strand))
    n = len(anchors)
    dp = [a.length for a in anchors]  # 以第i个锚点结尾的最长链长度，初始化为锚点自身长度
    prev = [-1]*n # 第i个锚点前驱下标
    q_ends = []  # (query_end, idx)，每个锚点qry区间终点和下标
    # dp_vals = [] # dp值
    for i, a in enumerate(anchors):
        # 二分查找所有query_end < a.q_start且ref_end < a.r_start的最大dp
        max_dp = 0
        max_j = -1
        for j in range(len(q_ends)):
            q_end, idx = q_ends[j]
            ref_end = anchors[idx].r_start + anchors[idx].length - 1
            if q_end < a.q_start:
                if dp[idx] > max_dp:
                    max_dp = dp[idx]
                    max_j = idx
        if max_j != -1:
            dp[i] = max_dp + a.length
            prev[i] = max_j
        # 记录当前锚点的query_end
        q_ends.append((a.q_start + a.length - 1, i))
        # dp_vals.append(dp[i])
    # 回溯最长链
    idx = max(range(n), key=lambda i: dp[i]) if n else -1
    chain = []
    while idx != -1:
        chain.append(anchors[idx])
        idx = prev[idx]
    return chain[::-1]

# Smith-Waterman局部比对
def smith_waterman(q, r, match_score=2, mismatch_score=-1, gap_open=-2):
    m, n = len(q), len(r)
    H = [[0]*(n+1) for _ in range(m+1)]
    max_score = 0
    max_pos = (0, 0)
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = H[i-1][j-1] + (match_score if q[i-1] == r[j-1] else mismatch_score)
            delete = H[i-1][j] + gap_open
            insert = H[i][j-1] + gap_open
            H[i][j] = max(0, match, delete, insert)
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)
    i, j = max_pos
    q_end, r_end = i-1, j-1
    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i-1][j-1] + (match_score if q[i-1] == r[j-1] else mismatch_score):
            i -= 1
            j -= 1
        elif H[i][j] == H[i-1][j] + gap_open:
            i -= 1
        else:
            j -= 1
    q_start, r_start = i, j
    return (q_start, q_end, r_start, r_end)

# gap区间局部比对，统计gap信息
def local_align(q: str, r: str, q_offset=0, r_offset=0, max_gap=1000, strand=1, gap_stats=None):
    results = []
    if len(q) == 0 or len(r) == 0:
        return results
    if gap_stats is not None:
        gap_stats['G'] += max(len(q), len(r))
        gap_stats['L'] = max(gap_stats['L'], len(q), len(r))
    if len(q) > max_gap or len(r) > max_gap:
        q1, r1 = q[:max_gap//2], r[:max_gap//2]
        q2, r2 = q[-max_gap//2:], r[-max_gap//2:]
        results.extend(local_align(q1, r1, q_offset, r_offset, max_gap, strand, gap_stats))
        results.extend(local_align(q2, r2, q_offset+len(q)-len(q2), r_offset+len(r)-len(r2), max_gap, strand, gap_stats))
        return results
    q_start, q_end, r_start, r_end = smith_waterman(q, r)
    if q_end >= q_start and r_end >= r_start:
        results.append((q_offset+q_start, q_offset+q_end, r_offset+r_start, r_offset+r_end))
    return results

# 合并连续区间
def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    cur = list(intervals[0])
    for t in intervals[1:]:
        if t[0] == cur[1] + 1 and t[2] == cur[3] + 1:
            cur[1] = t[1]
            cur[3] = t[3]
        else:
            merged.append(tuple(cur))
            cur = list(t)
    merged.append(tuple(cur))
    return merged

# 主流程
def match(query: str, ref: str, k: int = 8, w: int = 8) -> List[Tuple[int,int,int,int]]:
    """
    :param query: query序列
    :param ref: reference序列
    :param k: kmer长度
    :param w: minimizer窗口长度
    :return: 合并后的匹配区间列表
    """
    m, n = len(query), len(ref)
    # 1. 提取minimizer锚点并合并连续锚点
    anchors = build_anchors(query, ref, k, w)
    print(f"query长度m={m}, reference长度n={n}, 锚点总数A={len(anchors)}")

    # 2. O(A log A) DP求最长无重叠主链
    chain = dp(anchors)
    # 3. gap区间局部比对
    result = []
    last_q = last_r = None
    last_strand = None
    gap_stats = {'G': 0, 'L': 0}  # gap区间总长度G，单个gap最大长度L
    for anchor in chain:
        if last_q is not None and last_r is not None:
            # gap区间
            if anchor.strand == 1:
                q_gap = query[last_q+1:anchor.q_start]
            else:
                q_gap = reverse_complement(query[anchor.q_start+anchor.length:last_q]) if last_q > anchor.q_start+anchor.length else ''
            r_gap = ref[last_r+1:anchor.r_start]
            if q_gap and r_gap:
                gap_results = local_align(q_gap, r_gap, last_q+1 if anchor.strand==1 else anchor.q_start+anchor.length, last_r+1, 1000, anchor.strand, gap_stats)
                result.extend(gap_results)
        if anchor.strand == 1:
            result.append((anchor.q_start, anchor.q_start+anchor.length-1, anchor.r_start, anchor.r_start+anchor.length-1))
        else:
            q_len = len(query)
            rc_q_start = anchor.q_start
            rc_q_end = anchor.q_start+anchor.length-1
            result.append((rc_q_start, rc_q_end, anchor.r_start, anchor.r_start+anchor.length-1))
        last_q = anchor.q_start + anchor.length - 1 if anchor.strand == 1 else anchor.q_start
        last_r = anchor.r_start + anchor.length - 1
        last_strand = anchor.strand
    # 末尾gap不处理
    merged_result = merge_intervals(result)
    K = len(merged_result)
    print(f"gap区间总长度G={gap_stats['G']}, 单个gap最大长度L={gap_stats['L']}, 匹配区间数K={K}")
    return merged_result 