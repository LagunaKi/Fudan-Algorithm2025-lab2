import edlib
import argparse

def rc(s):
    map_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    l = []
    for c in s:
        l.append(map_dict[c])
    l = l[::-1]
    return ''.join(l)

def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    A = ref[ref_st: ref_en]
    a = query[query_st: query_en]
    _a = rc(query[query_st: query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def get_first(x):
    return x[0]

def calculate_value(ans_list, ref, query):
    # ans_list: [(q0, q1, r0, r1), ...]
    if not ans_list or len(ans_list[0]) != 4:
        return 0
    editdistance = 0
    aligned = 0
    preend = 0
    points = sorted(ans_list, key=get_first)
    for onetuple in points:
        query_st, query_en, ref_st, ref_en = onetuple
        if preend > query_st:
            return 0
        if query_en - query_st < 30:
            continue
        preend = query_en
        ed = calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
        if ed / max(1, query_en - query_st) > 0.1:
            continue
        editdistance += ed
        aligned += query_en - query_st
    return max(aligned - editdistance, 0)

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
    parser = argparse.ArgumentParser(description='DNA区间答案评分工具')
    parser.add_argument('--input', type=str, default=None, help='输入文件（格式同input.txt）')
    args = parser.parse_args()
    if args.input:
        query, ref = read_input_txt(args.input)
    else:
        print("请输入ref序列:")
        ref = input().strip()
        print("请输入query序列:")
        query = input().strip()
    print("请输入区间答案（如[(0, 299, 0, 299), ...]）:")
    ans_str = input().strip()
    try:
        ans_list = eval(ans_str)
    except Exception as e:
        print("区间格式解析失败", e)
        return
    score = calculate_value(ans_list, ref, query)
    print(f"Score: {score}")

if __name__ == "__main__":
    main() 