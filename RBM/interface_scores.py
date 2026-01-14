"""
Aggregate per-interface scores and identify best matches.

This module consolidates scores from multiple sources (IPS, ICS, QS_best, DockQ, lDDT, TM-score)
for each interface pair and identifies the best matching interface pairs based on maximum scores.
"""

import os
import sys
import numpy as np
from .utils import get_models


def get_pair2scores(model, target, model2qsbest, output_dir):
    """
    Load IPS and ICS scores for a model and combine with QS_best scores.
    
    Reads both IPS and ICS result files for a model, merges them with QS_best scores,
    and calculates interface weights based on interface size. Returns a dictionary
    mapping interface pairs to their combined score results and weights.
    """
    IPS_file = output_dir + '/' + target + '/' + 'IPS' + '/' + model + '.result'
    ICS_file = output_dir + '/' + target + '/' + 'ICS' + '/' + model + '.result'
    pair2results = {}
    cate = ''
    pair = ''
    weight = 0
    results = []
    with open(IPS_file, 'r') as f1, open(ICS_file, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        assert len(lines1) == len(lines2)
        for line1, line2 in zip(lines1, lines2):
            if line1[0] == '>':
                if cate and pair and weight:
                    pair2results[(cate, pair)] = [weight, results]
                cate = line1[1:].split()[0]
                pair = line1[1:].split()[1]
                weight = (int(line1[1:].split()[2].split(':')[0]) + int(line1[1:].split()[2].split(':')[1])) / 2
                results = []
            else:
                words1 = line1.split()
                words2 = line2.split()
                match_pair = words1[0] + ':' + words1[1]
                qsbest = model2qsbest[model][(cate, pair, match_pair)]
                results.append([words1[0], words1[1], float(words1[4]), float(words2[4]), qsbest])
    if cate and pair and weight:
        pair2results[(cate, pair)] = [weight, results]
    return pair2results


def save_interface_scores(input_dir, target, name, output_dir):
    """
    Aggregate and save per-interface scores for all models.
    
    Reads scores from IPS, ICS, QS_best, DockQ, lDDT, and TM-score result files.
    For each interface pair, identifies the best matching pair based on maximum
    scores across all metrics. Writes consolidated per-interface scores with
    best matches and interface weights. Handles both reference (REF) and prediction (MOD) interfaces.
    """
    models = get_models(input_dir, target, name)
    fp = open(output_dir + '/' + target + '/QS_best/' + target + '.result','r')
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

    fp = open(output_dir + '/' + target + '/DockQ/' + target + '.result','r')
    Rpair2dockq = {}
    Mpair2dockq = {}
    for line in fp:
        words = line.split()
        model = words[0]
        try:
            Rpair2dockq[model]
        except KeyError:
            Rpair2dockq[model] = {}
        try:
            Mpair2dockq[model]
        except KeyError:
            Mpair2dockq[model] = {}

        Rpair = words[1]
        Mpair = words[2]
        try:
            Rpair2dockq[model][Rpair].append([Mpair, float(words[3])])
        except KeyError:
            Rpair2dockq[model][Rpair] = [[Mpair, float(words[3])]]
        try:
            Mpair2dockq[model][Mpair].append([Rpair, float(words[3])])
        except KeyError:
            Mpair2dockq[model][Mpair] = [[Rpair, float(words[3])]]
    fp.close()

    fp = open(output_dir + '/' + target + '/lDDT/' + target + '.result','r')
    Rpair2lddt = {}
    Mpair2lddt = {}
    for line in fp:
        words = line.split()
        model = words[0]
        try:
            Rpair2lddt[model]
        except KeyError:
            Rpair2lddt[model] = {}
        try:
            Mpair2lddt[model]
        except KeyError:
            Mpair2lddt[model] = {}
        
        Rpair = words[1]
        Mpair = words[2]
        try:
            Rpair2lddt[model][Rpair].append([Mpair, float(words[3])])
        except KeyError:
            Rpair2lddt[model][Rpair] = [[Mpair, float(words[3])]]
        try:
            Mpair2lddt[model][Mpair].append([Rpair, float(words[3])])
        except KeyError:
            Mpair2lddt[model][Mpair] = [[Rpair, float(words[3])]]
    fp.close()

    fp = open(output_dir + '/' + target + '/TMscore/' + target + '.result','r')
    Rpair2tm = {}
    Mpair2tm = {}
    for line in fp:
        words = line.split()
        model = words[0]
        try:
            Rpair2tm[model]
        except KeyError:
            Rpair2tm[model] = {}
        try:
            Mpair2tm[model]
        except KeyError:
            Mpair2tm[model] = {}

        Rpair = words[1]
        Mpair = words[2]
        try:
            Rpair2tm[model][Rpair].append([Mpair, float(words[3])])
        except KeyError:
            Rpair2tm[model][Rpair] = [[Mpair, float(words[3])]]
        try:
            Mpair2tm[model][Mpair].append([Rpair, float(words[3])])
        except KeyError:
            Mpair2tm[model][Mpair] = [[Rpair, float(words[3])]]
    fp.close()


    good_count = 0
    bad_count = 0
    os.makedirs(output_dir + '/' + target + '/per_interface_scores', exist_ok=True)
    rp = open(output_dir + '/' + target + '/per_interface_scores/' + target + '.result','w')
    rp.write('model\tcategory\tchainpair\tmatchpair\tweight\tips\tics\tqsbest\tdockq\tlddt\ttm\n')
    for model in models:        
        pair2results = get_pair2scores(model, target, model2qsbest, output_dir)
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
                        print ('error\t' + target + '\t' + model + '\t' + cate + '\t' + pair)

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
                        print ('error\t' + target + '\t' + model + '\t' + cate + '\t' + pair)
    rp.close()
    print (target + '\t' + str(good_count) + '\t' + str(bad_count))

# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     save_interface_scores(input_dir, target, name, output_dir)