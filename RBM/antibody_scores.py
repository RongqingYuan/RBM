import re
import os
import sys
import math
import numpy as np



# if target[1] == '0':
#     fp = open('unknown_stoic/' + target + '/' + target + '.txt', 'r')
# else:
#     fp = open('known_stoic/' + target + '/' + target + '.txt', 'r')
def get_models(input_dir, target, name):
    fp = open(input_dir + '/' + target + '/' + name + '.txt', 'r')
    start = 0
    models = set([])
    for line in fp:
        words = line.split()
        if words:
            if words[0] == '#':
                start = 1
            elif start and len(words) > 1:
                models.add(words[1])
    fp.close()
    return models


def get_weighted_average(scores):
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

def save_antibody_scores(input_dir, target, name, output_dir, chainAs, chainBs):
    models = get_models(input_dir, target, name)
    fp = open(output_dir + '/' + target + '/' + 'per_interface_scores' + '/' + target + '.result', 'r')
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
    if not os.path.exists(output_dir + '/' + target + '/' + 'per_model_scores'):
        os.makedirs(output_dir + '/' + target + '/' + 'per_model_scores')

    rp = open(output_dir + '/' + target + '/' + 'per_model_scores' + '/' + target + '.result', 'w')
    rp.write('model\tips\tics\tqsbest\tdockq\n')
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

if __name__ == "__main__":
    input_dir = sys.argv[1]
    target = sys.argv[2]
    name = sys.argv[3]
    output_dir = sys.argv[4]
    chainAs = []
    for chain in sys.argv[5]:
        chainAs.append(chain)
    chainBs = []
    for chain in sys.argv[6]:
        chainBs.append(chain)

    save_antibody_scores(input_dir, target, name, output_dir, chainAs, chainBs)

# else:
#     fp = open('step16/' + target + '.result', 'r')
#     rp = open('step18/' + target + '.result', 'w')
#     rp.write('model\tips\tics\tqsbest\tdockq\n')
#     counts = set([])
#     for line in fp:
#         words = line.split()
#         model = words[0]

#         if model in models:
#             dockqs = []
#             qsbests = []
#             icss = []
#             ipss = []

#             good_inds = set([])
#             for counti, item in enumerate(re.split(r'[;,]', words[3])):
#                 chain1 = item[0]
#                 chain2 = item[1]
#                 if chain1 in chainAs and chain2 in chainBs:
#                     good_inds.add(counti)
#                 if chain1 in chainBs and chain2 in chainAs:
#                     good_inds.add(counti)

#             for counti, item in enumerate(re.split(r'[;,]', words[6])):
#                 if counti in good_inds:
#                     if item == 'None':
#                         dockqs.append(0.0)
#                     else:
#                         dockqs.append(float(item))
#             for counti, item in enumerate(re.split(r'[;,]', words[7])):
#                 if counti in good_inds:
#                     if item == 'None':
#                         qsbests.append(0.0)
#                     else:
#                         qsbests.append(float(item))
#             for counti, item in enumerate(re.split(r'[;,]', words[8])):
#                 if counti in good_inds:
#                     if item == 'None':
#                         icss.append(0.0)
#                     else:
#                         icss.append(float(item))
#             for counti, item in enumerate(re.split(r'[;,]', words[9])):
#                 if counti in good_inds:
#                     if item == 'None':
#                         ipss.append(0.0)
#                     else:
#                         ipss.append(float(item))

#             if not ipss or not icss or not qsbests or not dockqs:
#                 print (len(ipss), len(icss), len(qsbests), len(dockqs))
#             else:
#                 ave_ips = round(np.mean(ipss), 4)
#                 ave_ics = round(np.mean(icss), 4)
#                 ave_qsbest = round(np.mean(qsbests), 4)
#                 ave_dockq = round(np.mean(dockqs), 4)
#                 rp.write(model + '\t' + str(ave_ips) + '\t' + str(ave_ics) + '\t' + str(ave_qsbest) + '\t' + str(ave_dockq) + '\n')
#                 counts.add(len(dockqs))
#                 counts.add(len(qsbests))
#                 counts.add(len(icss))
#                 counts.add(len(ipss))
#     fp.close()
#     rp.close()
#     if len(counts) != 1:
#         print (target, counts)
