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

# 先长后短kmer锚点提取
def build_anchors(query: str, ref: str) -> List[Anchor]:
    """
    先长后短kmer锚点提取，优先长kmer，短kmer只补未覆盖区间。
    """
    m = len(query)
    covered = [False] * m
    anchors = []
    # k_list可调整
    k_list = []
    k_val = int(len(ref) * 2 / 3)
    while k_val >= 10:
        k_list.append(k_val)
        k_val = int(k_val * 2 / 3)
    for k in k_list:
        # 构建ref的kmer哈希表
        ref_kmers = collections.defaultdict(list)
        for i in range(len(ref) - k + 1):
            kmer = ref[i:i+k]
            ref_kmers[kmer].append(i)
        # 正向kmer
        for j in range(len(query) - k + 1):
            if any(covered[j:j+k]):
                continue
            kmer = query[j:j+k]
            for r_pos in ref_kmers.get(kmer, []):
                anchors.append(Anchor(j, r_pos, k, 1))
                for t in range(j, j+k):
                    covered[t] = True
        # 反向kmer
        rc_query = reverse_complement(query)
        for j in range(len(rc_query) - k + 1):
            if any(covered[len(query)-(j+k):len(query)-j]):
                continue
            kmer = rc_query[j:j+k]
            for r_pos in ref_kmers.get(kmer, []):
                q_start = len(query) - (j + k)
                anchors.append(Anchor(q_start, r_pos, k, -1))
                for t in range(q_start, q_start+k):
                    covered[t] = True
    # 合并连续/模糊锚点
    return merge_anchors(anchors)

# 合并连续锚点
def merge_anchors(anchors: List[Anchor], max_gap: int = 8) -> List[Anchor]:
    """
    模糊合并：只要锚点方向相同，且query和ref的起点距离在max_gap以内，就合并为一个更长的锚点。
    """
    if not anchors:
        return []
    anchors.sort(key=lambda a: (a.q_start, a.r_start, a.strand))
    merged = [anchors[0]]
    for a in anchors[1:]:
        last = merged[-1]
        # 方向相同，且query和ref的起点距离都在max_gap以内
        if (a.strand == last.strand and
            0 < a.q_start - (last.q_start + last.length) <= max_gap and
            0 < a.r_start - (last.r_start + last.length) <= max_gap):
            # 合并时把gap也算进长度
            gap_q = a.q_start - (last.q_start + last.length)
            last.length += gap_q + a.length  # query和ref的gap视为匹配
        else:
            merged.append(a)
    return merged

# 用二分查找优化的DP，O(A log A)
def dp(anchors: List[Anchor]) -> List[Anchor]:
    """
    用二分查找优化的DP，保证主链无重叠，支持倒位。
    """
    # 按ref坐标排序，便于二分查找
    anchors.sort(key=lambda a: (a.q_start, a.r_start, a.strand))
    n = len(anchors)
    dp = [a.length for a in anchors]  # 以第i个锚点结尾的最长链长度，初始化为锚点自身长度
    prev = [-1]*n # 第i个锚点前驱下标
    q_ends = []  # (query_end, idx)，每个锚点qry区间终点和下标
    for i, a in enumerate(anchors):
        max_dp = 0
        max_j = -1
        for j in range(len(q_ends)):
            q_end, idx = q_ends[j]
            if q_end < a.q_start: # qry区间无重叠
                if dp[idx] > max_dp:
                    max_dp = dp[idx]
                    max_j = idx
        if max_j != -1:
            dp[i] = max_dp + a.length
            prev[i] = max_j
        # 记录当前锚点的query_end
        q_ends.append((a.q_start + a.length - 1, i))
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
def merge_intervals(intervals, max_gap=5):
    if not intervals:
        return []
    merged = []
    cur = list(intervals[0])
    for t in intervals[1:]:
        # 模糊合并：只要query和ref的起点距离上一区间终点都在max_gap以内
        if 0 < t[0] - cur[1] <= max_gap and 0 < t[2] - cur[3] <= max_gap:
            cur[1] = t[1]
            cur[3] = t[3]
        else:
            merged.append(tuple(cur))
            cur = list(t)
    merged.append(tuple(cur))
    return merged

# 主流程
def match(query: str, ref: str) -> List[Tuple[int,int,int,int]]:
    """
    :param query: query序列
    :param ref: reference序列
    :return: 合并后的匹配区间列表
    """
    m, n = len(query), len(ref)
    # 1. 提取minimizer锚点并合并连续锚点
    anchors = build_anchors(query, ref)
    print(f"query长度m={m}, reference长度n={n}, 锚点总数A={len(anchors)}")

    # 2. O(A log A) DP求最长无重叠主链
    chain = dp(anchors)
    print(f"dp找到chain={chain}")


    # 3. gap区间局部比对
    result = []
    last_q = last_r = None
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
    # 末尾gap不处理
    merged_result = merge_intervals(result)
    K = len(merged_result)
    print(f"gap区间总长度G={gap_stats['G']}, 单个gap最大长度L={gap_stats['L']}, 匹配区间数K={K}")
    return merged_result 