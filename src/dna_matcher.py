import collections
from typing import List, Tuple
import math
import random
import edlib

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

# 简单minhash实现
HASH_SEEDS = [17, 31, 47, 73, 97, 131, 193, 257, 389, 521]
def kmer_minhash(s, num_hash=5):
    # s: kmer字符串
    # num_hash: 取前几个hash函数
    hashes = []
    for i in range(num_hash):
        h = 0
        for c in s:
            h = (h * HASH_SEEDS[i] + ord(c)) % 1000000007
        hashes.append(h)
    return tuple(hashes)

# 模糊合并锚点，允许小gap，方向一致且query/ref区间近连续
def merge_anchors_fuzzy(anchors, max_gap=5):
    if not anchors:
        return []
    anchors.sort(key=lambda a: (a.strand, a.q_start, a.r_start))
    merged = [anchors[0]]
    for a in anchors[1:]:
        last = merged[-1]
        if (a.strand == last.strand and
            0 <= a.q_start - (last.q_start + last.length) <= max_gap and
            0 <= a.r_start - (last.r_start + last.length) <= max_gap):
            gap_q = a.q_start - (last.q_start + last.length)
            last.length += gap_q + a.length
        else:
            merged.append(a)
    return merged

# 修改build_anchors，合并锚点后返回
def build_anchors(query: str, ref: str, max_gap=20) -> List[Anchor]:
    m = len(query)
    anchors = []
    covered = [False] * m
    k_list = []
    k_val = int(len(ref) * 2 / 3)
    while k_val >= 10:
        k_list.append(k_val)
        k_val = int(k_val * 2 / 3)
    for k in k_list:
        # print(f"start collect kmer for k={k}")
        ref_minhash = collections.defaultdict(list)
        for i in range(len(ref) - k + 1):
            kmer = ref[i:i+k]
            sig = kmer_minhash(kmer)
            ref_minhash[sig].append(i)
        for j in range(len(query) - k + 1):
            if any(covered[j:j+k]):
                continue
            q_kmer = query[j:j+k]
            sig = kmer_minhash(q_kmer)
            for r_pos in ref_minhash.get(sig, []):
                anchors.append(Anchor(j, r_pos, k, 1))
                for t in range(j, j+k):
                    covered[t] = True
        rc_query = reverse_complement(query)
        for j in range(len(rc_query) - k + 1):
            q_kmer = rc_query[j:j+k]
            sig = kmer_minhash(q_kmer)
            q_start = len(query) - (j + k)
            if any(covered[q_start:q_start+k]):
                continue
            for r_pos in ref_minhash.get(sig, []):
                anchors.append(Anchor(q_start, r_pos, k, -1))
                for t in range(q_start, q_start+k):
                    covered[t] = True
        # print(f"anchor num: {len(anchors)}")
    # 合并锚点
    anchors = merge_anchors_fuzzy(anchors, max_gap=max_gap)
    return anchors

# 新的gap penalty函数
def penalty(gap_q, gap_r, anchor_len, alpha=1, gamma=0.1, bonus=10, max_gap=20):
    # 顺滑延申奖励
    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap:
        return -bonus
    # 重复变异容忍
    if gap_q > 0 and gap_r < 0:
        return 0
    # 其它情况
    slope_penalty = gamma * abs(math.log((max(gap_q+1,1e-6))/(max(gap_r+1,1e-6))))
    # gap惩罚被长anchor"稀释"
    gap_penalty = alpha * max(0, gap_q) / (anchor_len + 1)
    return gap_penalty + slope_penalty

# 修改dp，传递anchor长度
def dp(anchors: List[Anchor], alpha=1, gamma=0.1, bonus=10, max_gap=20) -> List[Anchor]:
    anchors.sort(key=lambda a: (a.q_start, a.r_start, a.strand))
    n = len(anchors)
    dp = [a.length for a in anchors]
    prev = [-1]*n
    for j in range(n):
        for i in range(j):
            if anchors[i].q_start + anchors[i].length <= anchors[j].q_start:
                gap_q = anchors[j].q_start - (anchors[i].q_start + anchors[i].length)
                gap_r = anchors[j].r_start - (anchors[i].r_start + anchors[i].length)
                p = penalty(gap_q, gap_r, anchors[j].length, alpha=alpha, gamma=gamma, bonus=bonus, max_gap=max_gap)
                score = dp[i] + anchors[j].length - p
                if score > dp[j]:
                    dp[j] = score
                    prev[j] = i
    idx = max(range(n), key=lambda i: dp[i]) if n else -1
    chain = []
    while idx != -1:
        chain.append(anchors[idx])
        idx = prev[idx]
    return chain[::-1]

def calc_edit_dist_ratio(query, ref, q0, q1, r0, r1):
    # 计算edit distance比例，支持反向
    qseq = query[q0:q1+1] if q0 <= q1 else reverse_complement(query[q1:q0+1])
    rseq = ref[r0:r1+1] if r0 <= r1 else reverse_complement(ref[r1:r0+1])
    if not qseq or not rseq:
        return 1.0
    ed = edlib.align(qseq, rseq)['editDistance']
    return ed / max(1, len(qseq))

# 合并锚点，要求方向一致、近连续、edit distance低
def merge_anchors_edlib(anchors, query, ref, edit_dist_thresh=0.1, max_gap=20):
    if not anchors:
        return []
    anchors.sort(key=lambda a: (a.strand, a.q_start, a.r_start))
    merged = [anchors[0]]
    for a in anchors[1:]:
        last = merged[-1]
        # 判断是否可合并
        if (a.strand == last.strand and
            0 <= a.q_start - (last.q_start + last.length) <= max_gap and
            0 <= a.r_start - (last.r_start + last.length) <= max_gap):
            # 合并后区间edit distance比例
            q0 = last.q_start
            q1 = a.q_start + a.length - 1
            r0 = last.r_start
            r1 = a.r_start + a.length - 1
            if calc_edit_dist_ratio(query, ref, q0, q1, r0, r1) <= edit_dist_thresh:
                last.length = (a.q_start + a.length) - last.q_start
                continue
        merged.append(a)
    return merged

# 找主链gap区间
def find_gaps(chain, m, n):
    gaps = []
    last_q = -1
    last_r = -1
    for a in chain:
        q_start = a.q_start
        r_start = a.r_start
        if last_q != -1 and last_r != -1:
            if q_start > last_q + 1 and r_start > last_r + 1:
                gaps.append(((last_q+1, q_start-1), (last_r+1, r_start-1)))
        last_q = a.q_start + a.length - 1
        last_r = a.r_start + a.length - 1
    return gaps

# gap区间递归edlib局部比对，返回高质量区间
def edlib_local_align(qseq, rseq, q_offset, r_offset, edit_dist_thresh=0.1, min_len=30):
    results = []
    if len(qseq) < min_len or len(rseq) < min_len:
        return results
    # edlib全局比对
    res = edlib.align(qseq, rseq, mode='HW')
    if res['editDistance'] / max(1, len(qseq)) <= edit_dist_thresh and len(qseq) >= min_len:
        results.append((q_offset, q_offset+len(qseq)-1, r_offset, r_offset+len(rseq)-1))
        return results
    # 否则递归分割
    mid = len(qseq) // 2
    results += edlib_local_align(qseq[:mid], rseq[:mid], q_offset, r_offset, edit_dist_thresh, min_len)
    results += edlib_local_align(qseq[mid:], rseq[mid:], q_offset+mid, r_offset+mid, edit_dist_thresh, min_len)
    return results

# 主流程
# 保证接口不变
def match(query: str, ref: str, max_gap: int = 50, merge_chain: bool = True,
          alpha=43, gamma=6, bonus=16):
    m, n = len(query), len(ref)
    # 1. 先长后短提取anchors
    anchors_obj = build_anchors(query, ref, max_gap=max_gap)
    # print(f"query长度m={m}, reference长度n={n}, 锚点总数A={len(anchors_obj)}")
    # 2. DP主链
    chain = dp(anchors_obj, alpha=alpha, gamma=gamma, bonus=bonus, max_gap=max_gap)
    # 3. 合并主链区间（方向一致、近连续、edit distance低）
    chain = merge_anchors_edlib(chain, query, ref, edit_dist_thresh=0.1, max_gap=max_gap)
    # 4. gap区间补充
    gaps = find_gaps(chain, m, n)
    anchors_gap = []
    final_ans = []
    for a in chain:
        if a.strand == 1:
            final_ans.append((a.q_start, a.q_start+a.length-1, a.r_start, a.r_start+a.length-1))
        else:
            rc_q_start = a.q_start
            rc_q_end = a.q_start+a.length-1
            final_ans.append((rc_q_start, rc_q_end, a.r_start, a.r_start+a.length-1))
    for (qgap, rgap) in gaps:
        q0, q1 = qgap
        r0, r1 = rgap
        # 用更短kmer补充anchors
        anchors_short = build_anchors(query[q0:q1+1], ref[r0:r1+1], max_gap=5)
        # 坐标映射回全局
        for a in anchors_short:
            a.q_start += q0
            a.r_start += r0
        anchors_gap.extend(anchors_short)
        # 递归edlib局部比对补充高质量区间
        high_quality = edlib_local_align(query[q0:q1+1], ref[r0:r1+1], q0, r0, edit_dist_thresh=0.1, min_len=30)
        final_ans.extend(high_quality)
    # 5. anchors合并去重
    all_anchors = anchors_obj + anchors_gap
    all_anchors = merge_anchors_edlib(all_anchors, query, ref, edit_dist_thresh=0.1, max_gap=max_gap)
    anchors_out = []
    for anchor in all_anchors:
        if anchor.strand == 1:
            anchors_out.append((anchor.q_start, anchor.q_start+anchor.length-1, anchor.r_start, anchor.r_start+anchor.length-1))
        else:
            rc_q_start = anchor.q_start
            rc_q_end = anchor.q_start+anchor.length-1
            anchors_out.append((rc_q_start, rc_q_end, anchor.r_start, anchor.r_start+anchor.length-1))
    # 6. final_ans去重、排序
    final_ans = list(set(final_ans))
    final_ans.sort(key=lambda x: min(x[0], x[1]))
    # print(f"匹配区间数K={len(final_ans)}")
    return {'final_ans': final_ans, 'anchors': anchors_out} 