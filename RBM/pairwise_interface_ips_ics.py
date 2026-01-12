import sys, json
from itertools import permutations
from multiprocessing import Pool
import os
from utils import get_models

def get_Rchain2resids(input_dir, target, name):
    fp = open(input_dir + '/' + target + '/' + name + '.pdb', 'r')
    Rchain2resids = {}
    for line in fp:
        if len(line) > 60:
            if line[:4] == 'ATOM':
                resid = int(line[22:26])
                chain = line[21]
                try:
                    Rchain2resids[chain].add(resid)
                except KeyError:
                    Rchain2resids[chain] = set([resid])
    fp.close()
    return Rchain2resids

def get_ips_and_ics(input_dir, target, name, model, output_dir):
    Rchain2resids = get_Rchain2resids(input_dir, target, name)
    with open(input_dir + '/' + target + '/' + 'ost' + '/' + model + '.json', 'r') as file:
        data = json.load(file)
    ref_chem_groups = data["chem_groups"]
    mod_chem_groups = data["chem_mapping"]
    if target[0] == 'T':
        groups = [[]]
        for i in range(len(ref_chem_groups)):
            if ref_chem_groups[i]:
                for chain in ref_chem_groups[i]:
                    groups[0].append(['ref', chain])
            if mod_chem_groups[i]:
                for chain in mod_chem_groups[i]:
                    groups[0].append(['mod', chain])

    else:
        groups = []
        for i in range(len(ref_chem_groups)):
            groups.append([])
            if ref_chem_groups[i]:
                for chain in ref_chem_groups[i]:
                    groups[i].append(['ref', chain])
            if mod_chem_groups[i]:
                for chain in mod_chem_groups[i]:
                    groups[i].append(['mod', chain])

    Rchain_pairs = set([])
    Rpair2contacts = {}
    reference_contacts = data['reference_contacts']
    for refcon in reference_contacts:
        chain1 = refcon[0].split('.')[0]
        resid1 = int(refcon[0].split('.')[1])
        chain2 = refcon[1].split('.')[0]
        resid2 = int(refcon[1].split('.')[1])
        if chain1 < chain2:
            Rchain_pairs.add((chain1, chain2))
        elif chain2 < chain1:
            Rchain_pairs.add((chain2, chain1))
        try:
            Rpair2contacts[(chain1, chain2)].add((resid1, resid2))
        except KeyError:
            Rpair2contacts[(chain1, chain2)] = set([(resid1, resid2)])
        try:
            Rpair2contacts[(chain2, chain1)].add((resid2, resid1))
        except KeyError:
            Rpair2contacts[(chain2, chain1)] = set([(resid2, resid1)])

    Mchain_pairs = set([])
    Mpair2contacts = {}
    for modcon in data['model_contacts']:
        chain1 = modcon[0].split('.')[0]
        resid1 = int(modcon[0].split('.')[1])
        chain2 = modcon[1].split('.')[0]
        resid2 = int(modcon[1].split('.')[1])
        if chain1 < chain2:
            Mchain_pairs.add((chain1, chain2))
        elif chain2 < chain1:
            Mchain_pairs.add((chain2, chain1))
        try:
            Mpair2contacts[(chain1, chain2)].add((resid1, resid2))
        except KeyError:
            Mpair2contacts[(chain1, chain2)] = set([(resid1, resid2)])
        try:
            Mpair2contacts[(chain2, chain1)].add((resid2, resid1))
        except KeyError:
            Mpair2contacts[(chain2, chain1)] = set([(resid2, resid1)])

    get_Rchains = set([])
    get_Mchains = set([])
    Rchain2Mchain = {}
    Mchain2Rchain = {}
    for group in groups:
        refchains = []
        modchains = []
        for item in group:
            if item[0] == 'ref':
                refchains.append(item[1])
            elif item[0] == 'mod':
                modchains.append(item[1])
        for refchain in refchains:
            for modchain in modchains:
                try:
                    Rchain2Mchain[refchain].add(modchain)
                except KeyError:
                    get_Rchains.add(refchain)
                    Rchain2Mchain[refchain] = set([modchain])
                try:
                    Mchain2Rchain[modchain].add(refchain)
                except KeyError:
                    get_Mchains.add(modchain)
                    Mchain2Rchain[modchain] = set([refchain])

    all_results = []
    for Rchain_pair in Rchain_pairs:
        Rcontacts = Rpair2contacts[Rchain_pair]
        Rresid1s = set([])
        Rresid2s = set([])
        for Rcontact in Rcontacts:
            Rresid1s.add(Rcontact[0])
            Rresid2s.add(Rcontact[1])
        Rchain1 = Rchain_pair[0]
        Rchain2 = Rchain_pair[1]

        map_results = []
        if Rchain1 in get_Rchains and Rchain2 in get_Rchains:
            for Mchain1 in Rchain2Mchain[Rchain1]:
                for Mchain2 in Rchain2Mchain[Rchain2]:
                    try:
                        Mcontacts = Mpair2contacts[(Mchain1, Mchain2)]
                    except KeyError:
                        Mcontacts = set([])
                    good_Mcontacts = set([])
                    good_Mresid1s = set([])
                    good_Mresid2s = set([])
                    for Mcontact in Mcontacts:
                        Mresid1 = Mcontact[0]
                        Mresid2 = Mcontact[1]
                        if Mresid1 in Rchain2resids[Rchain1] and Mresid2 in Rchain2resids[Rchain2]:
                            good_Mcontacts.add(Mcontact)
                            good_Mresid1s.add(Mresid1)
                            good_Mresid2s.add(Mresid2)

                    if good_Mcontacts:
                        Rcount = len(Rresid1s) + len(Rresid2s)
                        Mcount = len(good_Mresid1s) + len(good_Mresid2s)
                        RMcount = len(Rresid1s.intersection(good_Mresid1s)) + len(Rresid2s.intersection(good_Mresid2s))
                        ips = round(RMcount / (Rcount + Mcount - RMcount), 4)
                        overlap = Rcontacts.intersection(good_Mcontacts)
                        if overlap:
                            precision = len(overlap) / len(Rcontacts)
                            recall = len(overlap) / len(good_Mcontacts)
                            ics = round(2 * precision * recall / (precision + recall), 4)
                        else:
                            ics = 0.0
                        map_results.append([Mchain1, Mchain2, len(Rresid1s), len(Rresid2s), ips, ics])
                    else:
                        map_results.append([Mchain1, Mchain2, len(Rresid1s), len(Rresid2s), 0.0, 0.0])
        all_results.append(['reference', Rchain1, Rchain2, len(Rresid1s), len(Rresid2s), map_results])

    for Mchain_pair in Mchain_pairs:
        Mcontacts = Mpair2contacts[Mchain_pair]
        Mresid1s = set([])
        Mresid2s = set([])
        for Mcontact in Mcontacts:
            Mresid1s.add(Mcontact[0])
            Mresid2s.add(Mcontact[1])
        Mchain1 = Mchain_pair[0]
        Mchain2 = Mchain_pair[1]

        map_results = []
        if Mchain1 in get_Mchains and Mchain2 in get_Mchains:
            for Rchain1 in Mchain2Rchain[Mchain1]:
                for Rchain2 in Mchain2Rchain[Mchain2]:
                    try:
                        Rcontacts = Rpair2contacts[(Rchain1, Rchain2)]
                    except KeyError:
                        Rcontacts = set([])
                    Rresid1s = set([])
                    Rresid2s = set([])
                    for Rcontact in Rcontacts:
                        Rresid1s.add(Rcontact[0])
                        Rresid2s.add(Rcontact[1])

                    good_Mcontacts = set([])
                    good_Mresid1s = set([])
                    good_Mresid2s = set([])
                    for Mcontact in Mcontacts:
                        Mresid1 = Mcontact[0]
                        Mresid2 = Mcontact[1]
                        if Mresid1 in Rchain2resids[Rchain1] and Mresid2 in Rchain2resids[Rchain2]:
                            good_Mcontacts.add(Mcontact)
                            good_Mresid1s.add(Mresid1)
                            good_Mresid2s.add(Mresid2)

                    if Rcontacts and good_Mcontacts:
                        Rcount = len(Rresid1s) + len(Rresid2s)
                        Mcount = len(good_Mresid1s) + len(good_Mresid2s)
                        RMcount = len(Rresid1s.intersection(good_Mresid1s)) + len(Rresid2s.intersection(good_Mresid2s))
                        ips = round(RMcount / (Rcount + Mcount - RMcount), 4)
                        overlap = Rcontacts.intersection(good_Mcontacts)
                        if overlap:
                            precision = len(overlap) / len(Rcontacts)
                            recall = len(overlap) / len(good_Mcontacts)
                            ics = round(2 * precision * recall / (precision + recall), 4)
                        else:
                            ics = 0.0
                        map_results.append([Rchain1, Rchain2, len(good_Mresid1s), len(good_Mresid2s), ips, ics])
                    else:
                        map_results.append([Rchain1, Rchain2, len(good_Mresid1s), len(good_Mresid2s), 0.0, 0.0])
        all_results.append(['prediction', Mchain1, Mchain2, len(Mresid1s), len(Mresid2s), map_results])
    
    if not os.path.exists(output_dir + '/' + target + '/' + 'IPS'):
        os.makedirs(output_dir + '/' + target + '/' + 'IPS')
    if not os.path.exists(output_dir + '/' + target + '/' + 'ICS'):
        os.makedirs(output_dir + '/' + target + '/' + 'ICS')
    rp = open(output_dir + '/' + target + '/' + 'IPS' + '/' + model + '.result','w')
    for result in all_results:
        category = result[0]
        chain1 = result[1]
        chain2 = result[2]
        size1 = str(result[3])
        size2 = str(result[4])
        rp.write('>' + category + ' ' + chain1 + ':' + chain2 + ' ' + size1 + ':' + size2 + '\n')
        for item in result[5]:
            rp.write(item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\t' + str(item[3]) + '\t' + str(item[4])  + '\n')
    rp.close()
    rp = open(output_dir + '/' + target + '/' + 'ICS' + '/' + model + '.result','w')
    for result in all_results:
        category = result[0]
        chain1 = result[1]
        chain2 = result[2]
        size1 = str(result[3])
        size2 = str(result[4])
        rp.write('>' + category + ' ' + chain1 + ':' + chain2 + ' ' + size1 + ':' + size2 + '\n')
        for item in result[5]:
            rp.write(item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\t' + str(item[3]) + '\t' + str(item[5])  + '\n')
    rp.close()
    return 0

def save_ips_and_ics(input_dir, target, name, output_dir):
    models = get_models(input_dir, target, name)
    for model in models:
        get_ips_and_ics(input_dir, target, name, model, output_dir)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    target = sys.argv[2]
    name = sys.argv[3]
    output_dir = sys.argv[4]
    save_ips_and_ics(input_dir, target, name, output_dir)


