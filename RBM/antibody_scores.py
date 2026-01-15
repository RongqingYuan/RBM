"""
Calculate per-model scores for antibody structures.

This module aggregates interface-level scores specifically for antibody-antigen
interfaces. It filters interfaces based on antibody chain IDs (chainAs and chainBs)
and calculates model-level scores using weighted averaging. Only uses a subset
of scores (IPS, ICS, QS_best, DockQ) suitable for antibody evaluation.
"""

import os
import sys
import math
import numpy as np


def get_weighted_average(scores):
    """
    Calculate weighted average for antibody score types.
    
    Takes a list of score entries, each containing a weight and four score values
    (IPS, ICS, QS_best, DockQ). Returns the weighted average for each score type.
    """
    weight = []
    ips = []
    ics = []
    qsbest = []
    dockq = []
    for item in scores:
        weight.append(item[0])
        ips.append(float(item[1]))
        ics.append(float(item[2]))
        qsbest.append(float(item[3]))
        dockq.append(float(item[4]))
    ave_ips = sum(v * w for v, w in zip(ips, weight)) / sum(weight)
    ave_ics = sum(v * w for v, w in zip(ics, weight)) / sum(weight)
    ave_qsbest = sum(v * w for v, w in zip(qsbest, weight)) / sum(weight)
    ave_dockq = sum(v * w for v, w in zip(dockq, weight)) / sum(weight)
    return [ave_ips, ave_ics, ave_qsbest, ave_dockq]

def save_antibody_scores(model_name, output_dir, chainAs, chainBs):
    """
    Calculate and save per-model scores for antibody structures (average method).
    
    Filters per-interface scores to only include antibody-antigen interfaces
    based on chain IDs. Calculates weighted averages separately for reference (REF)
    and prediction (MOD) interfaces. The final model score is the average of REF
    and MOD scores. Uses log10 weighting for interface sizes.
    
    Args:
        model_name: Model name to process
        output_dir: Path to output directory (will use output_dir/model_name/)
        chainAs: List of antibody chain IDs
        chainBs: List of antigen chain IDs
    """
    models = [model_name]
    model_output_dir = os.path.join(output_dir, model_name)
    fp = open(os.path.join(model_output_dir, model_name + '.interface_scores'), 'r')
    model2Rscores = {}
    model2Mscores = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            model = words[0]
            weight = math.log10(float(words[4]))
            
            if weight:    
                getit = 0
                if words[1] == 'REF':
                    refpairs = words[2].split(',')
                    for refpair in refpairs:
                        if refpair.split(':')[0] in chainAs and refpair.split(':')[1] in chainBs:
                            getit = 1
                        if refpair.split(':')[0] in chainBs and refpair.split(':')[1] in chainAs:
                            getit = 1
                    if getit:
                        try:
                            model2Rscores[model].append([weight, words[5], words[6], words[7], words[8]])
                        except KeyError:
                            model2Rscores[model] = [[weight, words[5], words[6], words[7], words[8]]]

                elif words[1] == 'MOD':
                    refpairs = words[3].split(',')
                    for refpair in refpairs:
                        if refpair.split(':')[0] in chainAs and refpair.split(':')[1] in chainBs:
                            getit = 1
                        if refpair.split(':')[0] in chainBs and refpair.split(':')[1] in chainAs:
                            getit = 1
                    if getit:
                        try:
                            model2Mscores[model].append([weight, words[5], words[6], words[7], words[8]])
                        except KeyError:
                            model2Mscores[model] = [[weight, words[5], words[6], words[7], words[8]]]
    fp.close()

    rp = open(os.path.join(model_output_dir, model_name + '.antibody_scores'), 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\n')
    for model in models:
        checkR = 0
        checkM = 0
        try:
            model2Rscores[model]
            checkR = 1
        except KeyError:
            print(model_name + ' miss REF scores for ' + model)
        try:
            model2Mscores[model]
            checkM = 1
        except KeyError:
            print(model_name + ' miss MOD scores for ' + model)

        if checkR:
            [Rips, Rics, Rqsbest, Rdockq] = get_weighted_average(model2Rscores[model])
            if checkM:
                [Mips, Mics, Mqsbest, Mdockq] = get_weighted_average(model2Mscores[model])
            else:
                Mips, Mics, Mqsbest, Mdockq = 0.0, 0.0, 0.0, 0.0
            ips = round((Rips + Mips) / 2, 4)
            ics = round((Rics + Mics) / 2, 4)
            qsbest = round((Rqsbest + Mqsbest) / 2, 4)
            dockq = round((Rdockq + Mdockq) / 2, 4)
            rp.write(model + '\t' + str(ips) + '\t' + str(ics) + '\t' + str(qsbest) + '\t' + str(dockq) + '\n')
    rp.close()

def save_antibody_scores_v3(model_name, output_dir, chainAs, chainBs):
    """
    Calculate and save per-model scores for antibody structures (minimum method).
    
    Filters per-interface scores to only include antibody-antigen interfaces
    based on chain IDs. Calculates weighted averages separately for reference (REF)
    and prediction (MOD) interfaces. The final model score is the minimum of REF
    and MOD scores, penalizing models that perform poorly on either direction.
    
    Args:
        model_name: Model name to process
        output_dir: Path to output directory (will use output_dir/model_name/)
        chainAs: List of antibody chain IDs
        chainBs: List of antigen chain IDs
    """
    models = [model_name]
    model_output_dir = os.path.join(output_dir, model_name)
    fp = open(os.path.join(model_output_dir, model_name + '.interface_scores'), 'r')
    model2Rscores = {}
    model2Mscores = {}
    for countl, line in enumerate(fp):
        if countl:
            words = line.split()
            model = words[0]
            weight = math.log10(float(words[4]))
            
            if weight:    
                getit = 0
                if words[1] == 'REF':
                    refpairs = words[2].split(',')
                    for refpair in refpairs:
                        if refpair.split(':')[0] in chainAs and refpair.split(':')[1] in chainBs:
                            getit = 1
                        if refpair.split(':')[0] in chainBs and refpair.split(':')[1] in chainAs:
                            getit = 1
                    if getit:
                        try:
                            model2Rscores[model].append([weight, words[5], words[6], words[7], words[8]])
                        except KeyError:
                            model2Rscores[model] = [[weight, words[5], words[6], words[7], words[8]]]

                elif words[1] == 'MOD':
                    refpairs = words[3].split(',')
                    for refpair in refpairs:
                        if refpair.split(':')[0] in chainAs and refpair.split(':')[1] in chainBs:
                            getit = 1
                        if refpair.split(':')[0] in chainBs and refpair.split(':')[1] in chainAs:
                            getit = 1
                    if getit:
                        try:
                            model2Mscores[model].append([weight, words[5], words[6], words[7], words[8]])
                        except KeyError:
                            model2Mscores[model] = [[weight, words[5], words[6], words[7], words[8]]]
    fp.close()

    rp = open(os.path.join(model_output_dir, model_name + '.antibody_scores'), 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\n')
    for model in models:
        checkR = 0
        checkM = 0
        try:
            model2Rscores[model]
            checkR = 1
        except KeyError:
            print(model_name + ' miss REF scores for ' + model)
        try:
            model2Mscores[model]
            checkM = 1
        except KeyError:
            print(model_name + ' miss MOD scores for ' + model)

        if checkR:
            [Rips, Rics, Rqsbest, Rdockq] = get_weighted_average(model2Rscores[model])
            if checkM:
                [Mips, Mics, Mqsbest, Mdockq] = get_weighted_average(model2Mscores[model])
            else:
                Mips, Mics, Mqsbest, Mdockq = 0.0, 0.0, 0.0, 0.0
            ips = round(min(Rips, Mips), 4)
            ics = round(min(Rics, Mics), 4)
            qsbest = round(min(Rqsbest, Mqsbest), 4)
            dockq = round(min(Rdockq, Mdockq), 4)
            rp.write(model + '\t' + str(ips) + '\t' + str(ics) + '\t' + str(qsbest) + '\t' + str(dockq) + '\n')
    rp.close()

    
# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     chainAs = []
#     for chain in sys.argv[5]:
#         chainAs.append(chain)
#     chainBs = []
#     for chain in sys.argv[6]:
#         chainBs.append(chain)
#     save_antibody_scores(input_dir, target, name, output_dir, chainAs, chainBs)
