'''
python optuna_param_search.py
'''

import subprocess
import sys
import optuna
import numpy as np
import edlib

def get_rc(s):
    map_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    l = []
    for c in s:
        l.append(map_dict[c])
    l = l[::-1]
    return ''.join(l)
def rc(s):
    map_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    l = []
    for c in s:
        l.append(map_dict[c])
    l = l[::-1]
    return ''.join(l)

def get_points(tuples_str):
    data = []
    num = 0
    for c in tuples_str:
        if(ord('0') <= c <= ord('9')):
            num = num * 10 + c - ord('0')
        elif(ord(',') == c):
            data.append(num)
            num = 0
    if(num != 0):
        data.append(num)
    return data

def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    A = ref[ref_st: ref_en]
    a = query[query_st: query_en]
    _a = rc(query[query_st: query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def get_first(x):
    return x[0]

def calculate_value(tuples_str, ref, query):  
    slicepoints = np.array(get_points(tuples_str.encode()))
    if(len(slicepoints) > 0 and len(slicepoints) % 4 == 0):
        editdistance = 0
        aligned = 0
        preend = 0
        points = np.array(slicepoints).reshape((-1, 4)).tolist()
        points.sort(key=get_first)
        for onetuple in points:
            query_st, query_en, ref_st, ref_en = onetuple
            if(preend > query_st):
                return 0
            if(query_en - query_st < 30):
                continue
            preend = query_en
            if((calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)/len(query[query_st:query_en])) > 0.1):
                continue
            editdistance += calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
            aligned += len(query[query_st:query_en])
        return max(aligned - editdistance, 0)
    else:
        return 0

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

def score_output(output_file, ref, query):
    with open(output_file, encoding='utf-8') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('final_ans'):
            tuples_str = line.split('=', 1)[1].strip()
            break
    else:
        return 0
    return calculate_value(tuples_str, ref, query)

def objective(trial):
    max_gap = trial.suggest_int('max_gap', 10, 120)
    bonus = trial.suggest_float('bonus', 10, 60)
    alpha = trial.suggest_float('alpha', 10.0, 200.0)
    gamma = trial.suggest_float('gamma', 1, 15)
    cmd = f"python -m src --input input2.txt --output output2.txt --max_gap {max_gap} --bonus {bonus} --alpha {alpha} --gamma {gamma}"
    subprocess.run(cmd, shell=True, check=True)
    query, ref = read_input_txt('input2.txt')
    score = score_output('output2.txt', ref, query)
    print(f"Params: max_gap={max_gap}, bonus={bonus}, alpha={alpha}, gamma={gamma}, Score={score}")
    return -score  # optuna默认最小化，这里取负

if __name__ == '__main__':
    study = optuna.create_study()
    study.optimize(objective, n_trials=80)
    print('Best params:', study.best_params)
    print('Best score:', -study.best_value) 