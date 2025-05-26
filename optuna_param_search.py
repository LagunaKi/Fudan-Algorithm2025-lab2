'''
python optuna_param_search.py
'''

import subprocess
import sys
import optuna
import numpy as np
import edlib
from optuna.samplers import TPESampler
import threading

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

def objective(trial, param_space=None):
    # param_space: dict or None
    if param_space is None:
        max_gap = trial.suggest_int('max_gap', 10, 60)
        alpha = trial.suggest_float('alpha', 1, 100)
        gamma = trial.suggest_float('gamma', 0.1, 10)
        bonus = trial.suggest_float('bonus', 1, 30)
        slope_eps = trial.suggest_float('slope_eps', 0.05, 0.5)
    else:
        max_gap = trial.suggest_int('max_gap', param_space['max_gap'][0], param_space['max_gap'][1])
        alpha = trial.suggest_float('alpha', param_space['alpha'][0], param_space['alpha'][1])
        gamma = trial.suggest_float('gamma', param_space['gamma'][0], param_space['gamma'][1])
        bonus = trial.suggest_float('bonus', param_space['bonus'][0], param_space['bonus'][1])
        slope_eps = trial.suggest_float('slope_eps', param_space['slope_eps'][0], param_space['slope_eps'][1])
    cmd = [sys.executable, '-m', 'src', '--input', 'input2.txt', '--output', 'output2.txt',
           '--max_gap', str(max_gap), '--alpha', str(alpha), '--gamma', str(gamma), '--bonus', str(bonus), '--slope_eps', str(slope_eps)]
    subprocess.run(cmd, check=True)
    query, ref = read_input_txt('input2.txt')
    score = score_output('output2.txt', ref, query)
    print(f"Params: max_gap={max_gap}, alpha={alpha}, gamma={gamma}, bonus={bonus}, slope_eps={slope_eps}, Score={score}")
    return score

def main():
    # 第一阶段：粗调
    print("[Optuna] 第一阶段：粗调参数空间...")
    study1 = optuna.create_study(direction='maximize', sampler=TPESampler())
    study1.optimize(lambda trial: objective(trial, None), n_trials=40, n_jobs=4)
    best_trials = sorted(study1.trials, key=lambda t: -t.value)[:max(4, int(0.1*len(study1.trials)))]
    # 取前10%最优trial的参数范围
    def get_range(key):
        vals = [t.params[key] for t in best_trials]
        return (min(vals), max(vals))
    param_space = {
        'max_gap': get_range('max_gap'),
        'alpha': get_range('alpha'),
        'gamma': get_range('gamma'),
        'bonus': get_range('bonus'),
        'slope_eps': get_range('slope_eps'),
    }
    print(f"[Optuna] 第一阶段最优参数区间: {param_space}")
    print(f"[Optuna] 第一阶段最优分数: {study1.best_value}")
    # 第二阶段：细调
    print("[Optuna] 第二阶段：细调参数空间...")
    study2 = optuna.create_study(direction='maximize', sampler=TPESampler())
    study2.optimize(lambda trial: objective(trial, param_space), n_trials=40, n_jobs=4)
    print(f"[Optuna] 第二阶段最优参数: {study2.best_params}")
    print(f"[Optuna] 第二阶段最优分数: {study2.best_value}")

if __name__ == '__main__':
    main() 