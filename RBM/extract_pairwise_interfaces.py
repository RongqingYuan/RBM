import os
import sys
import json
from multiprocessing import Pool
from utils import get_models


# ============================================================================
# Configuration Constants
# ============================================================================

# Score filter ratio: only keep interface matches where at least one score
# is greater than (best_score * SCORE_FILTER_RATIO)
SCORE_FILTER_RATIO = 0.5



def get_Rchain2resids_and_Rchain2lines(input_dir, target, name):
    fp = open(input_dir + '/' + target + '/' + name + '.pdb', 'r')
    Rchain2lines = {}
    Rchain2resids = {}
    for line in fp:
        if len(line) > 60:
            if line[:4] == 'ATOM':
                resid = int(line[22:26])
                chain = line[21]
                try:
                    Rchain2resids[chain].add(resid)
                    Rchain2lines[chain].append(line)
                except KeyError:
                    Rchain2resids[chain] = set([resid])
                    Rchain2lines[chain] = [line]
    fp.close()
    return Rchain2resids, Rchain2lines


def get_model2qsbest(input_dir, target):
    fp = open(input_dir + '/' + target + '/' + 'QS_best' + '/' + target + '.result','r')
    model2qsbest = {}
    for line in fp:
        words = line.split()
        model = words[0]
        try:
            model2qsbest[model]
        except KeyError:
            model2qsbest[model] = {}
        cate = words[1]
        pair1 = words[2].replace(':','-')
        pair2 = words[3].replace(':','-')
        model2qsbest[model][(cate, pair1, pair2)] = float(words[4])
    fp.close()
    return model2qsbest


def get_pair2scores(model, target, model2qsbest, output_dir):
    IPS_file = output_dir + '/' + target + '/' + 'IPS' + '/' + model + '.result'
    ICS_file = output_dir + '/' + target + '/' + 'ICS' + '/' + model + '.result'
    pair2results = {}
    cate = ''
    pair = ''
    results = []
    with open(IPS_file, 'r') as f1, open(ICS_file, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        assert len(lines1) == len(lines2)
        for line1, line2 in zip(lines1, lines2):
            if line1[0] == '>':
                if cate and pair and results:
                    pair2results[(cate, pair)] = results
                cate = line1[1:].split()[0]
                pair = line1[1:].split()[1].replace(':','-')
                results = []
            else:
                words1 = line1.split()
                words2 = line2.split()
                match_pair = words1[0] + '-' + words1[1]
                qsbest = model2qsbest[model][(cate, pair, match_pair)]
                results.append([words1[0], words1[1], float(words1[4]), float(words2[4]), qsbest])
    if cate and pair and results:
        pair2results[(cate, pair)] = results
    return pair2results


def process_model(input_dir, model, target, name, output_dir):
    Rchain2resids, Rchain2lines = get_Rchain2resids_and_Rchain2lines(input_dir, target, name)
    model2qsbest = get_model2qsbest(output_dir, target)
    pair2results = get_pair2scores(model, target, model2qsbest, output_dir)
    cases = set([])
    for (cate, pair1) in pair2results.keys():
        results = pair2results[(cate, pair1)]
        best_ips = 0
        best_ics = 0
        best_qs = 0
        for item in results:
            if item[2] > best_ips:
                best_ips = item[2]
            if item[3] > best_ics:
                best_ics = item[3]
            if item[4] > best_qs:
                best_qs = item[4]

        for item in results:
            if item[2] > best_ips * SCORE_FILTER_RATIO or item[3] > best_ics * SCORE_FILTER_RATIO or item[4] > best_qs * SCORE_FILTER_RATIO:
                chain1A = pair1.split('-')[0]
                chain1B = pair1.split('-')[1]
                if cate == 'reference':
                    cases.add((pair1, item[0] + '-' + item[1]))
                elif cate == 'prediction':
                    if item[0] < item[1]:
                        cases.add((item[0] + '-' + item[1], chain1A + '-' + chain1B))
                    else:
                        cases.add((item[1] + '-' + item[0], chain1B + '-' + chain1A))

    fp = open(input_dir + '/' + target + '/' + 'model' + '/' + model, 'r')
    Mchain2lines = {}
    for line in fp:
        if len(line) > 60:
            if line[:4] == 'ATOM':
                resid = int(line[22:26])
                chain = line[21]
                try:
                    Mchain2lines[chain].append([resid, line])
                except KeyError:
                    Mchain2lines[chain] = [[resid, line]]
    fp.close()
    
    os.makedirs(output_dir + '/' + target  + '/' + 'pairwise_interfaces' + '/' + model, exist_ok=True)
    need_cases = []
    for case in cases:
        pair1 = case[0]
        pair2 = case[1]
        check1 = 0
        check2 = 0
        check3 = 0
        check4 = 0
        check5 = 0
        check6 = 0

        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_target.pdb'):
            check1 = 1
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_model.pdb'):
            check2 = 1
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.dockq'):
            check3 = 1
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.lddt'):
            check4 = 1
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.log'):
            check5 = 1
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.tm'):
            check6 = 1

        if check1 and check2 and check3 and check4 and check5 and check6:
            pass
        else:
            if check1:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_target.pdb')
            if check2:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_model.pdb')
            if check3:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.dockq')
            if check4:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.lddt')
            if check5:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.log')
            if check6:
                os.remove(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.tm')

            rp1 = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_target.pdb','w')
            rp2 = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '_model.pdb','w')
            Rchain1 = pair1.split('-')[0]
            Rchain2 = pair1.split('-')[1]
            Rresids1 = Rchain2resids[Rchain1]
            Rresids2 = Rchain2resids[Rchain2]
            Mchain1 = pair2.split('-')[0]
            Mchain2 = pair2.split('-')[1]
            for line in Rchain2lines[Rchain1]:
                rp1.write(line[:21] + 'A' + line[22:])
            for line in Rchain2lines[Rchain2]:
                rp1.write(line[:21] + 'B' + line[22:])
            for item in Mchain2lines[Mchain1]:
                if item[0] in Rresids1:
                    rp2.write(item[1][:21] + 'A' + item[1][22:])
            for item in Mchain2lines[Mchain2]:
                if item[0] in Rresids2:
                    rp2.write(item[1][:21] + 'B' + item[1][22:])
            rp1.close()
            rp2.close()
            need_cases.append(case)
    return model, need_cases



def process_model_same_chain(input_dir, model, target, name,output_dir):
    Rchain2resids, Rchain2lines = get_Rchain2resids_and_Rchain2lines(input_dir, target, name)
    model2qsbest = get_model2qsbest(output_dir, target)
    pair2results = get_pair2scores(model, target, model2qsbest, output_dir)
    cases = set([])
    for (cate, pair1) in pair2results.keys():
        results = pair2results[(cate, pair1)]
        best_ips = 0
        best_ics = 0
        best_qs = 0
        for item in results:
            if item[2] > best_ips:
                best_ips = item[2]
            if item[3] > best_ics:
                best_ics = item[3]
            if item[4] > best_qs:
                best_qs = item[4]

        for item in results:
            if item[2] > best_ips * SCORE_FILTER_RATIO or item[3] > best_ics * SCORE_FILTER_RATIO or item[4] > best_qs * SCORE_FILTER_RATIO:
                chain1A = pair1.split('-')[0]
                chain1B = pair1.split('-')[1]
                if cate == 'reference':
                    cases.add((pair1, item[0] + '-' + item[1]))
                elif cate == 'prediction':
                    if item[0] < item[1]:
                        cases.add((item[0] + '-' + item[1], chain1A + '-' + chain1B))
                    else:
                        cases.add((item[1] + '-' + item[0], chain1B + '-' + chain1A))

    fp = open(input_dir + '/' + target + '/' + 'model' + '/' + model, 'r')
    Mchain2lines = {}
    for line in fp:
        if len(line) > 60:
            if line[:4] == 'ATOM':
                resid = int(line[22:26])
                chain = line[21]
                try:
                    Mchain2lines[chain].append([resid, line])
                except KeyError:
                    Mchain2lines[chain] = [[resid, line]]
    fp.close()

    os.makedirs(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + model, exist_ok=True)
    need_cases = []
    for case in cases:
        pair1 = case[0]
        pair2 = case[1]
        rp1 = open(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + model + '/' + pair1 + '_' + pair2 + '_target.pdb','w')
        rp2 = open(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + model + '/' + pair1 + '_' + pair2 + '_model.pdb','w')
        Rchain1 = pair1.split('-')[0]
        Rchain2 = pair1.split('-')[1]
        Rresids1 = Rchain2resids[Rchain1]
        Rresids2 = Rchain2resids[Rchain2]
        Mchain1 = pair2.split('-')[0]
        Mchain2 = pair2.split('-')[1]

        for line in Rchain2lines[Rchain1]:
            rp1.write(line[:21] + 'A' + line[22:])
        for item in Mchain2lines[Mchain1]:
            if item[0] in Rresids1:
                rp2.write(item[1][:21] + 'A' + item[1][22:])
        shift = max(Rresids1) + 100

        for line in Rchain2lines[Rchain2]:
            resid = int(line[22:26]) + shift
            rp1.write(line[:21] + 'A' + str(resid).rjust(4, ' ') + line[26:])
        for item in Mchain2lines[Mchain2]:
            if item[0] in Rresids2:
                resid = int(item[1][22:26]) + shift
                rp2.write(item[1][:21] + 'A' + str(resid).rjust(4, ' ') + item[1][26:])
        rp1.close()
        rp2.close()
        need_cases.append(case)
    return [model, need_cases]

def save_pairwise_interfaces(input_dir, target, name, output_dir, score_filter_ratio):
    # Override the global constant with the passed parameter
    global SCORE_FILTER_RATIO
    SCORE_FILTER_RATIO = score_filter_ratio
    
    models = get_models(input_dir, target, name)
    if not os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces'):
        os.makedirs(output_dir + '/' + target + '/' + 'pairwise_interfaces')
    with open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + target + '.list','w') as f:
        for model in models:
            model, need_cases = process_model(input_dir, model, target, name, output_dir)
            for case in need_cases:
                f.write(model + '\t' + case[0] + '\t' + case[1] + '\n')
    if not os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt'):
        os.makedirs(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt')
    with open(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + target + '.list','w') as f:
        for model in models:
            model, need_cases = process_model_same_chain(input_dir, model, target, name, output_dir)
            for case in need_cases:
                f.write(model + '\t' + case[0] + '\t' + case[1] + '\n')


if __name__ == "__main__":
    input_dir = sys.argv[1]
    target = sys.argv[2]
    name = sys.argv[3]
    output_dir = sys.argv[4]
    save_pairwise_interfaces(input_dir, target, name, output_dir)



# models = get_models(input_dir, target, name)
# Rchain2resids, Rchain2lines = get_Rchain2resids_and_Rchain2lines(input_dir, target, name)
# model2qsbest = get_model2qsbest(output_dir, target)
# if not os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces'):
#     os.makedirs(output_dir + '/' + target + '/' + 'pairwise_interfaces')
# with open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + target + '.list','w') as f:
#     for model in models:
#         pair2results = get_pair2scores(model, target, model2qsbest, output_dir)
#         model, need_cases = process_model(input_dir, model, target, pair2results, Rchain2resids, Rchain2lines, output_dir)
#         for case in need_cases:
#             f.write(model + '\t' + case[0] + '\t' + case[1] + '\n')




# pool = Pool(processes = 32)
# processes = []
# for model in models:
# 	process = pool.apply_async(process_model, [model])
# 	processes.append(process)
# listp = open('step9/' + target + '.list','w')
# for process in processes:
#     [model, cases] = process.get()
#     for case in cases:
#         pair1 = case[0]
#         pair2 = case[1]
#         listp.write(model + '\t' + pair1 + '\t' + pair2 + '\n')
# listp.close()