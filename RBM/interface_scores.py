"""
Aggregate per-interface scores and identify best matches.

This module consolidates scores from multiple sources (IPS, ICS, QS_best, DockQ, lDDT, TM-score)
for each interface pair and identifies the best matching interface pairs based on maximum scores.
"""

import os
import sys
import numpy as np


def get_pair2scores(model_name, model2qsbest, output_dir):
    """
    Load IPS and ICS scores for a model and combine with QS_best scores.
    
    Reads both IPS and ICS result files for a model, merges them with QS_best scores,
    and calculates interface weights based on interface size. Returns a dictionary
    mapping interface pairs to their combined score results and weights.
    
    Args:
        model_name: Model name
        model2qsbest: Dictionary of QS_best scores
        output_dir: Path to output directory
    """
    model_output_dir = os.path.join(output_dir, model_name)
    IPS_file = os.path.join(model_output_dir, model_name + '.ips')
    ICS_file = os.path.join(model_output_dir, model_name + '.ics')
    pair2results = {}
    
    # Read IPS scores - new unified format: model_name category chainpair matchpair size1 size2 score
    ips_data = {}
    with open(IPS_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # parts: [model_name, category, chainpair, matchpair, size1, size2, score]
            category = parts[1]
            chainpair = parts[2]
            matchpair = parts[3]
            size1 = int(parts[4])
            size2 = int(parts[5])
            score = float(parts[6])
            key = (category, chainpair, matchpair)
            ips_data[key] = {'score': score, 'size1': size1, 'size2': size2}
    
    # Read ICS scores - new unified format: model_name category chainpair matchpair size1 size2 score
    ics_data = {}
    with open(ICS_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            category = parts[1]
            chainpair = parts[2]
            matchpair = parts[3]
            score = float(parts[6])
            key = (category, chainpair, matchpair)
            ics_data[key] = score
    
    # Combine IPS and ICS data
    for key in ips_data:
        category, chainpair, matchpair = key
        ips_score = ips_data[key]['score']
        ics_score = ics_data.get(key, 0.0)
        size1 = ips_data[key]['size1']
        size2 = ips_data[key]['size2']
        weight = (size1 + size2) / 2
        
        # Get QS_best score
        match_chain1, match_chain2 = matchpair.split(':')
        qsbest = model2qsbest[model_name][(category, chainpair, matchpair)]
        
        # Store in pair2results
        if (category, chainpair) not in pair2results:
            pair2results[(category, chainpair)] = [weight, []]
        
        pair2results[(category, chainpair)][1].append([match_chain1, match_chain2, ips_score, ics_score, qsbest])
    
    return pair2results


def save_interface_scores(model_name, output_dir):
    """
    Aggregate and save per-interface scores for a single model.
    
    Reads scores from IPS, ICS, QS_best, DockQ, lDDT, and TM-score result files.
    For each interface pair, identifies the best matching pair based on maximum
    scores across all metrics. Writes consolidated per-interface scores with
    best matches and interface weights. Handles both reference (REF) and prediction (MOD) interfaces.
    
    Args:
        model_name: Model name to process
        output_dir: Path to output directory (will use output_dir/model_name/)
    """
    models = [model_name]
    model_output_dir = os.path.join(output_dir, model_name)
    
    # Read QS_best scores
    fp = open(os.path.join(model_output_dir, model_name + '.qs'), 'r')
    model2qsbest = {}
    for line in fp:
        words = line.split()
        model = words[0]
        try:
            model2qsbest[model]
        except KeyError:
            model2qsbest[model] = {}
        cate = words[1]
        pair1 = words[2]
        pair2 = words[3]
        model2qsbest[model][(cate, pair1, pair2)] = float(words[4])
    fp.close()

    # Read DockQ scores
    fp = open(os.path.join(model_output_dir, model_name + '.dockq'), 'r')
    Rpair2dockq = {}
    Mpair2dockq = {}
    for line in fp:
        words = line.split()
        model = words[0]
        # New format: model category Rpair Mpair score
        category = words[1]
        Rpair = words[2]
        Mpair = words[3]
        score = float(words[4])
        
        try:
            Rpair2dockq[model]
        except KeyError:
            Rpair2dockq[model] = {}
        try:
            Mpair2dockq[model]
        except KeyError:
            Mpair2dockq[model] = {}

        try:
            Rpair2dockq[model][Rpair].append([Mpair, score])
        except KeyError:
            Rpair2dockq[model][Rpair] = [[Mpair, score]]
        try:
            Mpair2dockq[model][Mpair].append([Rpair, score])
        except KeyError:
            Mpair2dockq[model][Mpair] = [[Rpair, score]]
    fp.close()

    # Read lDDT scores
    fp = open(os.path.join(model_output_dir, model_name + '.lddt'), 'r')
    Rpair2lddt = {}
    Mpair2lddt = {}
    for line in fp:
        words = line.split()
        model = words[0]
        # New format: model category Rpair Mpair score
        category = words[1]
        Rpair = words[2]
        Mpair = words[3]
        score = float(words[4])
        
        try:
            Rpair2lddt[model]
        except KeyError:
            Rpair2lddt[model] = {}
        try:
            Mpair2lddt[model]
        except KeyError:
            Mpair2lddt[model] = {}
        
        try:
            Rpair2lddt[model][Rpair].append([Mpair, score])
        except KeyError:
            Rpair2lddt[model][Rpair] = [[Mpair, score]]
        try:
            Mpair2lddt[model][Mpair].append([Rpair, score])
        except KeyError:
            Mpair2lddt[model][Mpair] = [[Rpair, score]]
    fp.close()

    # Read TMscore scores
    fp = open(os.path.join(model_output_dir, model_name + '.tm'), 'r')
    Rpair2tm = {}
    Mpair2tm = {}
    for line in fp:
        words = line.split()
        model = words[0]
        # New format: model category Rpair Mpair score
        category = words[1]
        Rpair = words[2]
        Mpair = words[3]
        score = float(words[4])
        
        try:
            Rpair2tm[model]
        except KeyError:
            Rpair2tm[model] = {}
        try:
            Mpair2tm[model]
        except KeyError:
            Mpair2tm[model] = {}

        try:
            Rpair2tm[model][Rpair].append([Mpair, score])
        except KeyError:
            Rpair2tm[model][Rpair] = [[Mpair, score]]
        try:
            Mpair2tm[model][Mpair].append([Rpair, score])
        except KeyError:
            Mpair2tm[model][Mpair] = [[Rpair, score]]
    fp.close()


    good_count = 0
    bad_count = 0
    rp = open(os.path.join(model_output_dir, model_name + '.interface_scores'), 'w')
    rp.write('model\tcategory\tchainpair\tmatchpair\tweight\tips\tics\tqsbest\tdockq\tlddt\ttm\n')
    for model in models:        
        pair2results = get_pair2scores(model_name, model2qsbest, output_dir)
        for (cate, pair) in pair2results.keys():
            [weight, results] = pair2results[(cate, pair)]
            if not results:
                if cate == 'reference':
                    rp.write(model + '\tREF\t' + pair + '\tna\t' + str(weight) + '\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n')
                elif cate == 'prediction':
                    rp.write(model + '\tMOD\t' + pair + '\tna\t' + str(weight) + '\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n')
            
            else:
                best_ips = 0.0
                best_ics = 0.0
                best_qs = 0.0
                best_matches = set([])
                for item in results:
                    if item[2] > best_ips:
                        best_ips = item[2]
                        best_matches.add(item[0] + ':' + item[1])
                    if item[3] > best_ics:
                        best_ics = item[3]
                        best_matches.add(item[0] + ':' + item[1])
                    if item[4] > best_qs:
                        best_qs = item[4]
                        best_matches.add(item[0] + ':' + item[1])

                if cate == 'reference':
                    check = 0
                    try:
                        Rpair2dockq[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_dockq = Rpair2dockq[model][pair][0][1]
                        best_matches.add(Rpair2dockq[model][pair][0][0])
                        check += 1
                    except KeyError:
                        best_dockq = 0.0

                    try:
                        Rpair2lddt[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_lddt = Rpair2lddt[model][pair][0][1]
                        best_matches.add(Rpair2lddt[model][pair][0][0])
                        check += 1 
                    except KeyError:
                        best_lddt = 0.0

                    try:
                        Rpair2tm[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_tm = Rpair2tm[model][pair][0][1]
                        best_matches.add(Rpair2tm[model][pair][0][0])
                        check += 1
                    except KeyError:
                        best_tm = 0.0

                    if check == 3 or check == 0:
                        if check:
                            good_count += 1
                        else:
                            bad_count += 1
                        if best_matches:
                            rp.write(model + '\tREF\t' + pair + '\t' + ','.join(best_matches) + '\t' + str(weight) + '\t' + str(best_ips) + '\t' + str(best_ics) + '\t' + str(best_qs) + '\t' + str(best_dockq) + '\t' + str(best_lddt) + '\t' + str(best_tm) + '\n')
                        else:
                            rp.write(model + '\tREF\t' + pair + '\tna\t' + str(weight) + '\t' + str(best_ips) + '\t' + str(best_ics) + '\t' + str(best_qs) + '\t' + str(best_dockq) + '\t' + str(best_lddt) + '\t' + str(best_tm) + '\n')
                    else:
                        print('error\t' + model + '\t' + cate + '\t' + pair)

                elif cate == 'prediction':
                    check = 0
                    try:
                        Mpair2dockq[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_dockq = Mpair2dockq[model][pair][0][1]
                        best_matches.add(Mpair2dockq[model][pair][0][0])
                        check += 1
                    except KeyError:
                        best_dockq = 0.0

                    try:
                        Mpair2lddt[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_lddt = Mpair2lddt[model][pair][0][1]
                        best_matches.add(Mpair2lddt[model][pair][0][0])
                        check += 1
                    except KeyError:
                        best_lddt = 0.0

                    try:
                        Mpair2tm[model][pair].sort(key = lambda x:x[1], reverse = True)
                        best_tm = Mpair2tm[model][pair][0][1]
                        best_matches.add(Mpair2tm[model][pair][0][0])
                        check += 1
                    except KeyError:
                        best_tm = 0.0

                    if check == 3 or check == 0:
                        if check:
                            good_count += 1
                        else:
                            bad_count += 1
                        if best_matches:
                            rp.write(model + '\tMOD\t' + pair + '\t' + ','.join(best_matches) + '\t' + str(weight) + '\t' + str(best_ips) + '\t' + str(best_ics) + '\t' + str(best_qs) + '\t' + str(best_dockq) + '\t' + str(best_lddt) + '\t' + str(best_tm) + '\n')
                        else:
                            rp.write(model + '\tMOD\t' + pair + '\tna\t' + str(weight) + '\t' + str(best_ips) + '\t' + str(best_ics) + '\t' + str(best_qs) + '\t' + str(best_dockq) + '\t' + str(best_lddt) + '\t' + str(best_tm) + '\n')
                    else:
                        print('error\t' + model + '\t' + cate + '\t' + pair)
    rp.close()
    print(model_name + '\t' + str(good_count) + '\t' + str(bad_count))

# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     save_interface_scores(input_dir, target, name, output_dir)