import sys
import argparse
from .dna_matcher import match

# 读取input.txt格式
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
    parser = argparse.ArgumentParser(description='复杂DNA序列比对工具')
    parser.add_argument('--input', type=str, default=None, help='输入文件（格式同input.txt）')
    parser.add_argument('--output', type=str, default=None, help='输出文件')
    parser.add_argument('--max_gap', type=int, default=34, help='gap区间最大合并距离')
    parser.add_argument('--bonus', type=float, default=40, help='顺滑延申奖励')
    parser.add_argument('--alpha', type=float, default=100, help='query gap惩罚系数')
    parser.add_argument('--gamma', type=float, default=3, help='斜率惩罚系数')
    args = parser.parse_args()

    print("algorithm start")
    # 1. 获取query和ref
    if args.input:
        query, ref = read_input_txt(args.input)
    else:
        print('请输入query序列:')
        query = input().strip()
        print('请输入reference序列:')
        ref = input().strip()
    # 2. 匹配
    result = match(query, ref, max_gap=args.max_gap, bonus=args.bonus, alpha=args.alpha, gamma=args.gamma)
    # 3. 输出
    final_ans = result['final_ans']
    anchors = result['anchors']
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(f"final_ans = {final_ans}\n")
            f.write(f"anchors = {anchors}\n")
    else:
        print('比对结果:')
        print(f"final_ans = {final_ans}")
        print(f"anchors = {anchors}")
    print("algorithm complete")

if __name__ == '__main__':
    main() 

'''
不指定输入输出文件：  python -m src
指定输入文件：       python -m src --input input.txt
指定输出文件：       python -m src --output output.txt
同时指定输入输出文件：python -m src --input input.txt --output output.txt
'''