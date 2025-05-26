import collections
from typing import List
import math

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

# 修改build_anchors，使用更丰富的k值列表，从len(ref)*0.9开始，逐步减少到5
def build_anchors(query: str, ref: str, max_gap=20, min_k=3, k_ratio=2/3) -> List[Anchor]:
    m = len(query)
    anchors = []
    covered = [False] * m
    k_list = []
    k_val = int(len(ref) * 0.9)
    while k_val >= min_k:
        k_list.append(k_val)
        k_val = int(k_val * k_ratio)
    long_anchors = []  # 存储长anchor区间 (q0, q1, r0, r1, strand)
    for idx, k in enumerate(k_list):
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
            sig = kmer_minhash(q_kmer)
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

def match(query: str, ref: str, max_gap: int = 50, alpha=43, gamma=6, bonus=16, slope_eps=0.3, min_k=3, k_ratio=2/3):
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
                    slope = abs(delta_q / (delta_r + 1e-6)) if delta_r != 0 else 0
                    slope_ok = abs(slope - 1) <= slope_eps
                    if 0 <= gap_q <= max_gap and 0 <= gap_r <= max_gap and slope_ok:
                        merged_intervals[-1] = [last[0], q1, r0, last[3], strand]
                        continue
            merged_intervals.append([q0, q1, r0, r1, strand])
    final_ans = [(q0, q1, r0, r1, strand) for q0, q1, r0, r1, strand in merged_intervals]
    final_ans = list(set(final_ans))
    final_ans.sort(key=lambda x: min(x[0], x[1]))
    anchors_out = []
    for a in anchors_obj:
        q0 = a.q_start
        q1 = a.q_start + a.length - 1
        r0 = a.r_start
        r1 = a.r_start + a.length - 1
        anchors_out.append((q0, q1, r0, r1, a.strand))
    return {'final_ans': final_ans, 'anchors': anchors_out}