'''
不指定输入输出文件：  python -m src
指定输入文件：       python -m src --input input.txt
指定输出文件：       python -m src --output output.txt
同时指定输入输出文件：python -m src --input input.txt --output output.txt
'''

import sys
import argparse
from .dna_matcher import match

# 读取input（txt）
def read_input_txt(filename):
    with open(filename, encoding='utf-8') as f:
        lines = f.readlines()
    ref = ''
    qry = ''
    for line in lines:
        if line.startswith('ref'):
            ref = line.split('=', 1)[1].strip()
        elif line.startswith('qry'):
            qry = line.split('=', 1)[1].strip()
    return qry, ref

def main():
    parser = argparse.ArgumentParser(description='复杂DNA比对')
    parser.add_argument('--input', type=str, default=None, help='输入文件')
    parser.add_argument('--output', type=str, default=None, help='输出文件')
    parser.add_argument('--max_gap', type=int, default=49, help='gap最大合并距离')
    parser.add_argument('--alpha', type=float, default=28.59, help='query gap惩罚系数')
    parser.add_argument('--gamma', type=float, default=5.28, help='斜率惩罚系数')
    parser.add_argument('--bonus', type=float, default=21.98, help='顺滑延伸奖励')
    parser.add_argument('--slope_eps', type=float, default=0.32, help='合并斜率偏差阈值')
    parser.add_argument('--min_k', type=int, default=9, help='最小kmer长度')
    parser.add_argument('--k_ratio', type=float, default=0.68, help='kmer递减公比')
    parser.add_argument('--ex_max_edit_ratio', type=float, default=0.13, help='主链区间延伸最大允许错配比例')
    args = parser.parse_args()

    print("algorithm start")
    if args.input:
        query, ref = read_input_txt(args.input)
    else:
        print('请输入query序列:')
        query = input().strip()
        print('请输入reference序列:')
        ref = input().strip()
    result = match(query, ref, max_gap=args.max_gap, alpha=args.alpha, gamma=args.gamma, bonus=args.bonus, slope_eps=args.slope_eps, min_k=args.min_k, k_ratio=args.k_ratio, ex_max_edit_ratio=args.ex_max_edit_ratio)
    final_ans = result['final_ans']
    anchors = result['anchors']
    # 输出的final_ans为答案，也就是主链。其它三个输出为锚点、主链正反向列表、锚点正反向列表，用于绘图和检查。
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(f"final_ans = {[t[:4] for t in final_ans]}\n")
            f.write(f"anchors = {[t[:4] for t in anchors]}\n")
            f.write(f"final_ans_strand = {[t[4] for t in final_ans]}\n")
            f.write(f"anchors_strand = {[t[4] for t in anchors]}\n")
    else:
        print('比对结果:')
        print(f"final_ans = {[t[:4] for t in final_ans]}")
        print(f"anchors = {[t[:4] for t in anchors]}")
        print(f"final_ans_strand = {[t[4] for t in final_ans]}")
        print(f"anchors_strand = {[t[4] for t in anchors]}")
    print("algorithm complete")

if __name__ == '__main__':
    main() 

