import collections
from typing import List
import math
import numpy as np
from numba import njit

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

def smart_k_list(query, ref, min_k=3, k_ratio=2/3, min_anchors=10):
    k_list = []
    k_val = int(len(ref) * 0.9)
    prev_k = None
    while k_val >= min_k:
        # 统计该k下的锚点数（正向+反向）
        ref_minhash = collections.defaultdict(list)
        for i in range(len(ref) - k_val + 1):
            kmer = ref[i:i+k_val]
            sig = kmer_minhash(kmer)
            ref_minhash[sig].append(i)
        anchor_count = 0
        for j in range(len(query) - k_val + 1):
            q_kmer = query[j:j+k_val]
            sig = kmer_minhash(q_kmer)
            anchor_count += len(ref_minhash.get(sig, []))
        rc_query = reverse_complement(query)
        for j in range(len(rc_query) - k_val + 1):
            q_kmer = rc_query[j:j+k_val]
            sig = kmer_minhash(q_kmer)
            anchor_count += len(ref_minhash.get(sig, []))
        if prev_k is not None and anchor_count < min_anchors and k_val > min_k:
            k_mid = int((k_val + prev_k) / 2)
            if k_mid not in k_list and k_mid >= min_k:
                k_list.append(k_mid)
        k_list.append(k_val)
        prev_k = k_val
        k_val = int(k_val * k_ratio)
    k_list = sorted(set(k_list), reverse=True)
    return k_list

# anchor延伸（numba加速）
@njit
def extend_anchor_numba(query_arr, ref_arr, q_start, r_start, length, q_len, r_len):
    # 向左延伸
    l = 0
    while q_start - l - 1 >= 0 and r_start - l - 1 >= 0:
        if query_arr[q_start - l - 1] != ref_arr[r_start - l - 1]:
            break
        l += 1
    # 向右延伸
    r = 0
    while q_start + length + r < q_len and r_start + length + r < r_len:
        if query_arr[q_start + length + r] != ref_arr[r_start + length + r]:
            break
        r += 1
    return q_start - l, r_start - l, length + l + r

def build_anchors(query: str, ref: str, max_gap=20, min_k=3, k_ratio=2/3, min_anchors=10) -> List[Anchor]:
    m = len(query)
    anchors = []
    covered = [False] * m
    k_list = smart_k_list(query, ref, min_k=min_k, k_ratio=k_ratio, min_anchors=min_anchors)
    long_anchors = []  # 存储长anchor区间 (q0, q1, r0, r1, strand)
    for idx, k in enumerate(k_list):
        ref_minhash = collections.defaultdict(list)
        for i in range(len(ref) - k + 1):
            kmer = ref[i:i+k]
            sig = tuple(kmer_minhash(kmer))
            ref_minhash[sig].append(i)
        for j in range(len(query) - k + 1):
            if any(covered[j:j+k]):
                continue
            q_kmer = query[j:j+k]
            sig = tuple(kmer_minhash(q_kmer))
            for r_pos in ref_minhash.get(sig, []):
                q0, q1 = j, j+k-1
                r0, r1 = r_pos, r_pos+k-1
                if idx > 0:  # 不是最长k时，判断是否被长anchor共线覆盖
                    found = False
                    for la in long_anchors:
                        if (la[4] == 1 and la[0] <= q0 and la[1] >= q1 and la[2] <= r0 and la[3] >= r1):
                            found = True
                            break
                    if found:
                        continue
                anchors.append(Anchor(j, r_pos, k, 1))
                for t in range(j, j+k):
                    covered[t] = True
                if idx == 0:
                    long_anchors.append((q0, q1, r0, r1, 1))
        rc_query = reverse_complement(query)
        for j in range(len(rc_query) - k + 1):
            if any(covered[len(query) - (j + k):len(query) - j]):
                continue
            q_kmer = rc_query[j:j+k]
            sig = tuple(kmer_minhash(q_kmer))
            q_start = len(query) - (j + k)
            for r_pos in ref_minhash.get(sig, []):
                q0, q1 = q_start, q_start+k-1
                r0, r1 = r_pos, r_pos+k-1
                if idx > 0:
                    found = False
                    for la in long_anchors:
                        if (la[4] == -1 and la[0] <= q0 and la[1] >= q1 and la[2] <= r0 and la[3] >= r1):
                            found = True
                            break
                    if found:
                        continue
                anchors.append(Anchor(q_start, r_pos, k, -1))
                for t in range(q_start, q_start+k):
                    covered[t] = True
                if idx == 0:
                    long_anchors.append((q0, q1, r0, r1, -1))
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

# 修改dp，放宽gap宽容度
def dp(anchors: List[Anchor], alpha=1, gamma=0.1, bonus=10, max_gap=40) -> List[Anchor]:
    anchors.sort(key=lambda a: (a.q_start, a.r_start, a.strand))
    n = len(anchors)
    dp = [a.length for a in anchors]
    prev = [-1]*n
    for j in range(n):
        for i in range(j):
            if anchors[i].q_start + anchors[i].length <= anchors[j].q_start:
                if anchors[i].strand == anchors[j].strand:
                    if anchors[j].strand == 1:
                        gap_q = anchors[j].q_start - (anchors[i].q_start + anchors[i].length)
                        gap_r = anchors[j].r_start - (anchors[i].r_start + anchors[i].length)
                    else:
                        gap_q = anchors[j].q_start - (anchors[i].q_start + anchors[i].length)
                        gap_r = anchors[i].r_start - (anchors[j].r_start + anchors[j].length)
                else:
                    gap_q = anchors[j].q_start - (anchors[i].q_start + anchors[i].length)
                    gap_r = abs(anchors[j].r_start - (anchors[i].r_start + anchors[i].length))
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

def edit_distance_simple(s1, s2, strand=1, r0=None, r1=None):
    # s1: query区间
    # s2: ref区间
    # strand: 1正向，-1反向
    # r0, r1: ref区间起止（反向时用）
    if strand == 1:
        return sum(a != b for a, b in zip(s1, s2))
    else:
        # 反向时s2为ref[r0:r1+1]，s1[i]应比对ref[r1-(i-q0)]
        return sum(s1[i] != s2[len(s2)-1-i] for i in range(len(s1)))

def match(query: str, ref: str, max_gap: int = 50, alpha=43, gamma=6, bonus=16, slope_eps=0.3, min_k=3, k_ratio=2/3, ex_max_edit_ratio=0.1):
    anchors_obj = build_anchors(query, ref, max_gap=max_gap, min_k=min_k, k_ratio=k_ratio)
    chain = dp(anchors_obj, alpha=alpha, gamma=gamma, bonus=bonus, max_gap=max_gap)
    merged_intervals = []
    for a in chain:
        q0 = a.q_start
        q1 = a.q_start + a.length - 1
        r0 = a.r_start
        r1 = a.r_start + a.length - 1
        strand = a.strand
        if not merged_intervals:
            merged_intervals.append([q0, q1, r0, r1, strand])
        else:
            last = merged_intervals[-1]
            if strand == last[4]:
                # 正向
                if strand == 1:
                    gap_q = q0 - (last[1] + 1)
                    gap_r = r0 - (last[3] + 1)
                    delta_q = q1 - last[0]
                    delta_r = r1 - last[2]
                    slope = delta_q / (delta_r + 1e-6) if delta_r != 0 else 0
                    slope_ok = abs(slope - 1) <= slope_eps
                    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap and slope_ok:
                        merged_intervals[-1] = [last[0], q1, last[2], r1, strand]
                        continue
                # 反向
                else:
                    gap_q = q0 - (last[1] + 1)
                    gap_r = last[2] - (r1 + 1)
                    # 合并后区间应为 (last[0], q1, r0, last[3])
                    delta_q = q1 - last[0]
                    delta_r = last[3] - r0
                    slope = delta_q / (delta_r + 1e-6) if delta_r != 0 else 0
                    slope_ok = abs(slope - 1) <= slope_eps
                    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap and slope_ok:
                        merged_intervals[-1] = [last[0], q1, r0, last[3], strand]
                        continue
            merged_intervals.append([q0, q1, r0, r1, strand])

    # 主链区间两端延伸，允许少量错配，不能导致query区间重叠
    extended_intervals = []
    prev_q1 = -1
    for idx, (q0, q1, r0, r1, strand) in enumerate(merged_intervals):
        orig_q0, orig_q1, orig_r0, orig_r1 = q0, q1, r0, r1
        # 向左延伸
        while True:
            if strand == 1:
                nq0 = q0 - 1
                nr0 = r0 - 1
                nr1 = r1
            else:
                nq0 = q0 - 1
                nr0 = r0
                nr1 = r1 + 1
            if nq0 <= prev_q1 or nq0 < 0:
                break
            if strand == 1:
                if nr0 < 0:
                    break
                new_query = query[nq0:q1+1]
                new_ref = ref[nr0:nr1+1]
            else:
                if nr1 >= len(ref):
                    break
                new_query = query[nq0:q1+1]
                new_ref = ref[nr0:nr1+1]
            if len(new_query) != len(new_ref):
                break
            ed = edit_distance_simple(new_query, new_ref, strand, nr0, nr1)
            if ed / len(new_query) > ex_max_edit_ratio:
                break
            q0 = nq0
            if strand == 1:
                r0 = nr0
            else:
                r1 = nr1
        # 向右延伸
        while True:
            if strand == 1:
                nq1 = q1 + 1
                nr0 = r0
                nr1 = r1 + 1
            else:
                nq1 = q1 + 1
                nr0 = r0 - 1
                nr1 = r1
            if nq1 >= len(query) or nq1 <= q1:
                break
            # 不能与下一区间重叠
            if idx + 1 < len(merged_intervals) and nq1 >= merged_intervals[idx + 1][0]:
                break
            if strand == 1:
                if nr1 >= len(ref):
                    break
                new_query = query[q0:nq1+1]
                new_ref = ref[nr0:nr1+1]
            else:
                if nr0 < 0:
                    break
                new_query = query[q0:nq1+1]
                new_ref = ref[nr0:nr1+1]
            if len(new_query) != len(new_ref):
                break
            ed = edit_distance_simple(new_query, new_ref, strand, nr0, nr1)
            if ed / len(new_query) > ex_max_edit_ratio:
                break
            q1 = nq1
            if strand == 1:
                r1 = nr1
            else:
                r0 = nr0
        extended_intervals.append([q0, q1, r0, r1, strand])
        prev_q1 = q1

    # 延伸后再执行一次merge操作，逻辑与前面一致
    merged2 = []
    for a in extended_intervals:
        q0, q1, r0, r1, strand = a
        if not merged2:
            merged2.append([q0, q1, r0, r1, strand])
        else:
            last = merged2[-1]
            if strand == last[4]:
                if strand == 1:
                    gap_q = q0 - (last[1] + 1)
                    gap_r = r0 - (last[3] + 1)
                    delta_q = q1 - last[0]
                    delta_r = r1 - last[2]
                    slope = delta_q / (delta_r + 1e-6) if delta_r != 0 else 0
                    slope_ok = abs(slope - 1) <= slope_eps
                    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap and slope_ok:
                        merged2[-1] = [last[0], q1, last[2], r1, strand]
                        continue
                else:
                    gap_q = q0 - (last[1] + 1)
                    gap_r = last[2] - (r1 + 1)
                    delta_q = q1 - last[0]
                    delta_r = last[3] - r0
                    slope = delta_q / (delta_r + 1e-6) if delta_r != 0 else 0
                    slope_ok = abs(slope - 1) <= slope_eps
                    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap and slope_ok:
                        merged2[-1] = [last[0], q1, r0, last[3], strand]
                        continue
            merged2.append([q0, q1, r0, r1, strand])

    final_ans = list(set(tuple(x) for x in merged2))
    final_ans.sort(key=lambda x: min(x[0], x[1]))
    anchors_out = []
    for a in anchors_obj:
        q0 = a.q_start
        q1 = a.q_start + a.length - 1
        r0 = a.r_start
        r1 = a.r_start + a.length - 1
        anchors_out.append((q0, q1, r0, r1, a.strand))
    return {'final_ans': final_ans, 'anchors': anchors_out}