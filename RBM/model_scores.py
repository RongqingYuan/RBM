"""
Calculate per-model scores from per-interface scores.

This module aggregates interface-level scores (IPS, ICS, QS_best, DockQ, lDDT, TM-score)
into model-level scores using weighted averaging. Multiple versions (v1, v2, v3) implement
different strategies for combining reference and prediction interface scores.
"""

import sys
import math
import os
from utils import get_models



def get_weighted_average(scores):
    """
    Calculate weighted average for multiple score types.
    
    Takes a list of score entries, each containing a weight and six score values
    (IPS, ICS, QS_best, DockQ, lDDT, TM-score). Returns the weighted average
    for each score type.
    """
    weight = []
    ips = []
    ics = []
    qsbest = []
    dockq = []
    lddt = []
    tm = []
    for item in scores:
        weight.append(item[0])
        ips.append(float(item[1]))
        ics.append(float(item[2]))
        qsbest.append(float(item[3]))
        dockq.append(float(item[4]))
        lddt.append(float(item[5]))
        tm.append(float(item[6]))
    ave_ips = round(sum(v * w for v, w in zip(ips, weight)) / sum(weight), 4)
    ave_ics = round(sum(v * w for v, w in zip(ics, weight)) / sum(weight), 4)
    ave_qsbest = round(sum(v * w for v, w in zip(qsbest, weight)) / sum(weight), 4)
    ave_dockq = round(sum(v * w for v, w in zip(dockq, weight)) / sum(weight), 4)
    ave_lddt = round(sum(v * w for v, w in zip(lddt, weight)) / sum(weight), 4)
    ave_tm = round(sum(v * w for v, w in zip(tm, weight)) / sum(weight), 4)
    return [ave_ips, ave_ics, ave_qsbest, ave_dockq, ave_lddt, ave_tm]

def save_model_scores_v1(input_dir, target, name, output_dir, interface_weight):
    """
    Calculate and save per-model scores using v1 aggregation method (average of both directions).
    
    Reads per-interface scores and calculates weighted averages separately for
    reference (REF) and prediction (MOD) interfaces. The final model score is the
    average of REF and MOD scores. Missing MOD scores are set to 0.0.
    """
    models = get_models(input_dir, target, name)
    fp = open(output_dir + '/' + target + '/per_interface_scores/' + target + '.result', 'r')
    model2Rscores = {}
    model2Mscores = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            model = words[0]
            if interface_weight == 'log2':
                weight = math.log2(float(words[4]))
            elif interface_weight == 'log10':
                weight = math.log10(float(words[4]))
            elif interface_weight == 'linear':
                weight = float(words[4])
            else:
                raise ValueError("Invalid interface weight")
            if weight:
                if words[1] == 'REF':
                    try:
                        model2Rscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Rscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
                elif words[1] == 'MOD':
                    try:
                        model2Mscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Mscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
    fp.close()

    os.makedirs(output_dir + '/' + target + '/per_model_scores', exist_ok=True)
    rp = open(output_dir + '/' + target + '/per_model_scores/' + target + '.result', 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\tlddt\ttm\n')
    for model in models:
        checkR = 0
        checkM = 0
        try:
            model2Rscores[model]
            checkR = 1
        except KeyError:
            print (target + ' miss REF scores for ' + model)
        try:
            model2Mscores[model]
            checkM = 1
        except KeyError:
            print (target + ' miss MOD scores for ' + model)

        if checkR:
            [Rips, Rics, Rqsbest, Rdockq, Rlddt, Rtm] = get_weighted_average(model2Rscores[model])
            if checkM:
                [Mips, Mics, Mqsbest, Mdockq, Mlddt, Mtm] = get_weighted_average(model2Mscores[model])
            else:
                Mips, Mics, Mqsbest, Mdockq, Mlddt, Mtm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            ips = round((Rips + Mips) / 2, 4)
            ics = round((Rics + Mics) / 2, 4)
            qsbest = round((Rqsbest + Mqsbest) / 2, 4)
            dockq = round((Rdockq + Mdockq) / 2, 4)
            lddt = round((Rlddt + Mlddt) / 2, 4)
            tm = round((Rtm + Mtm) / 2, 4)
            rp.write(model + '\t' + str(ips) + '\t' + str(ics) + '\t' + str(qsbest) + '\t' + str(dockq) + '\t' + str(lddt) + '\t' + str(tm) + '\n')
    rp.close()
def save_model_scores_v2(input_dir, target, name, output_dir, interface_weight):
    """
    Calculate and save per-model scores using v2 aggregation method (all interfaces together).
    
    Reads per-interface scores and combines all REF and MOD interfaces together,
    then calculates a single weighted average across all interfaces. This treats
    reference and prediction interfaces equally in the aggregation.
    """
    models = get_models(input_dir, target, name)
    fp = open(output_dir + '/' + target + '/per_interface_scores/' + target + '.result', 'r')
    model2Rscores = {}
    model2Mscores = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            model = words[0]
            if interface_weight == 'log2':
                weight = math.log2(float(words[4]))
            elif interface_weight == 'log10':
                weight = math.log10(float(words[4]))
            elif interface_weight == 'linear':
                weight = float(words[4])
            else:
                raise ValueError("Invalid interface weight")
            if weight:
                if words[1] == 'REF':
                    try:
                        model2Rscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Rscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
                elif words[1] == 'MOD':
                    try:
                        model2Mscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Mscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
    fp.close()

    os.makedirs(output_dir + '/' + target + '/per_model_scores_v2', exist_ok=True)
    rp = open(output_dir + '/' + target + '/per_model_scores_v2/' + target + '.result', 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\tlddt\ttm\n')
    for model in models:
        checkR = 0
        checkM = 0
        try:
            model2Rscores[model]
            checkR = 1
        except KeyError:
            print (target + ' miss REF scores for ' + model)
        try:
            model2Mscores[model]
            checkM = 1
        except KeyError:
            print (target + ' miss MOD scores for ' + model)

        if checkR:
            all_R_interface_scores = model2Rscores[model]
            if checkM:
                all_M_interface_scores = model2Mscores[model]
            else:
                all_M_interface_scores = []

            all_interface_scores = all_R_interface_scores + all_M_interface_scores
            [ips, ics, qsbest, dockq, lddt, tm] = get_weighted_average(all_interface_scores)
            rp.write(model + '\t' + str(ips) + '\t' + str(ics) + '\t' + str(qsbest) + '\t' + str(dockq) + '\t' + str(lddt) + '\t' + str(tm) + '\n')
    rp.close()
def save_model_scores_v3(input_dir, target, name, output_dir, interface_weight):
    """
    Calculate and save per-model scores using v3 aggregation method (minimum score).
    
    Reads per-interface scores and calculates weighted averages separately for
    reference (REF) and prediction (MOD) interfaces. The final model score is the
    minimum of REF and MOD scores, penalizing models that perform poorly on either
    reference or prediction interfaces.
    """
    models = get_models(input_dir, target, name)
    fp = open(output_dir + '/' + target + '/per_interface_scores/' + target + '.result', 'r')
    model2Rscores = {}
    model2Mscores = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            model = words[0]
            if interface_weight == 'log2':
                weight = math.log2(float(words[4]))
            elif interface_weight == 'log10':
                weight = math.log10(float(words[4]))
            elif interface_weight == 'linear':
                weight = float(words[4])
            else:
                raise ValueError("Invalid interface weight")
            if weight:
                if words[1] == 'REF':
                    try:
                        model2Rscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Rscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
                elif words[1] == 'MOD':
                    try:
                        model2Mscores[model].append([weight, words[5], words[6], words[7], words[8], words[9], words[10]])
                    except KeyError:
                        model2Mscores[model] = [[weight, words[5], words[6], words[7], words[8], words[9], words[10]]]
    fp.close()

    os.makedirs(output_dir + '/' + target + '/per_model_scores_v3', exist_ok=True)
    rp = open(output_dir + '/' + target + '/per_model_scores_v3/' + target + '.result', 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\tlddt\ttm\n')
    for model in models:
        checkR = 0
        checkM = 0
        try:
            model2Rscores[model]
            checkR = 1
        except KeyError:
            print (target + ' miss REF scores for ' + model)
        try:
            model2Mscores[model]
            checkM = 1
        except KeyError:
            print (target + ' miss MOD scores for ' + model)

        if checkR:
            [Rips, Rics, Rqsbest, Rdockq, Rlddt, Rtm] = get_weighted_average(model2Rscores[model])
            if checkM:
                [Mips, Mics, Mqsbest, Mdockq, Mlddt, Mtm] = get_weighted_average(model2Mscores[model])
            else:
                Mips, Mics, Mqsbest, Mdockq, Mlddt, Mtm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            ips = round(min(Rips, Mips), 4)
            ics = round(min(Rics, Mics), 4)
            qsbest = round(min(Rqsbest, Mqsbest), 4)
            dockq = round(min(Rdockq, Mdockq), 4)
            lddt = round(min(Rlddt, Mlddt), 4)
            tm = round(min(Rtm, Mtm), 4)
            rp.write(model + '\t' + str(ips) + '\t' + str(ics) + '\t' + str(qsbest) + '\t' + str(dockq) + '\t' + str(lddt) + '\t' + str(tm) + '\n')
    rp.close()

    
# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     save_model_scores_v1(input_dir, target, name, output_dir)
#     save_model_scores_v2(input_dir, target, name, output_dir)
#     save_model_scores_v3(input_dir, target, name, output_dir)