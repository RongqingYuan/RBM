"""
Calculate QS_best scores for pairwise interfaces.

This module computes QS_best (Quality Score best) by comparing interface contacts
between reference and model structures. It uses a distance-based approach to identify
contacts and calculates the best matching score between reference and model interfaces.
"""

import sys
import json
import os
import numpy as np
from itertools import permutations
from multiprocessing import Pool
from .utils import get_model_name_from_path, parse_reference_pdb


def get_qs_best(reference_pdb_path, model_pdb_path, ost_json_path, model_name, cutoff):
    """
    Calculate QS_best score for a single model.
    
    This function processes both reference and model structures to compute QS_best
    (Quality Score best) for interface pairs. It first extracts contacts from both
    structures using distance-based criteria, then identifies chain mappings from
    OpenStructure JSON files.
    
    The function handles two types of interface pairs:
    - Reference interfaces: compares model chain pairs that match reference chain pairs
    - Prediction interfaces: compares reference chain pairs that match model chain pairs
    
    QS_best is calculated as:
    overlap / max(len(reference_contacts), len(model_contacts))
    where overlap is the intersection of reference and model contacts. Only contacts
    within the reference chain residue ranges are considered for model contacts.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        model_pdb_path: Full path to model PDB file
        ost_json_path: Full path to OST JSON file with chain mappings and contacts
        model_name: Model name to use in output (e.g., 'H0208TS014_1')
        target: Target name for categorization (e.g., 'H0208')
        cutoff: Distance cutoff in Angstroms for identifying contacts
    
    Returns:
        tuple: (model_name, all_results) where all_results contains interface pairs
               and their QS_best scores organized by category
    """
    Rchain2resids, all_Rcontacts, _ = parse_reference_pdb(reference_pdb_path, cutoff)
    with open(ost_json_path, 'r') as file:
        data = json.load(file)

    ref_chem_groups = data["chem_groups"]
    mod_chem_groups = data["chem_mapping"]
    
    # Determine grouping strategy based on model name (T-series uses single group, H-series uses multiple)
    use_single_group = model_name[0] == 'T' if model_name else False
    
    if use_single_group:
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
    for refcon in all_Rcontacts:
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

    fp = open(model_pdb_path, 'r')

    Mchains = []
    Mresid2CAcoor = {}
    Mresid2coors = {}
    all_Mcontacts = []
    Mchain2resids = {}
    for line in fp:
        if len(line) > 60:
            if line[:4] == 'ATOM':
                resid = int(line[22:26])
                chain = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom = line[12:16].strip()
                element = line[76:78].strip()
                # residue_name = line[17:20].strip()
                try:
                    Mchain2resids[chain].add(resid)
                    Mresid2coors[chain]
                    Mresid2CAcoor[chain]
                except KeyError:
                    Mchains.append(chain)
                    Mchain2resids[chain] = set([resid])
                    Mresid2coors[chain] = {}
                    Mresid2CAcoor[chain] = {}
                # if residue_name == 'GLY':
                #     if atom == 'CA':
                #         Mresid2CAcoor[chain][resid] = [x, y, z]
                # else:
                #     if atom == 'CB':
                #         Mresid2CAcoor[chain][resid] = [x, y, z]
                if atom == 'CA':
                    Mresid2CAcoor[chain][resid] = [x, y, z]
                if element != 'H':
                    try:
                        Mresid2coors[chain][resid].append([x, y, z])
                    except KeyError:
                        Mresid2coors[chain][resid] = [[x, y, z]]
    fp.close()

    Mchains.sort()
    for c1, chain1 in enumerate(Mchains):
        for c2, chain2 in enumerate(Mchains):
            if c1 < c2:
                for resid1 in Mchain2resids[chain1]:
                    for resid2 in Mchain2resids[chain2]:
                        CAcoor1 = Mresid2CAcoor[chain1][resid1]
                        CAcoor2 = Mresid2CAcoor[chain2][resid2]
                        CAdist = ((CAcoor1[0] - CAcoor2[0]) ** 2 + (CAcoor1[1] - CAcoor2[1]) ** 2 + (CAcoor1[2] - CAcoor2[2]) ** 2) ** 0.5
                        if CAdist < CA_DISTANCE_PREFILTER:
                            coor1s = Mresid2coors[chain1][resid1]
                            coor2s = Mresid2coors[chain2][resid2]
                            dists = []
                            for coor1 in coor1s:
                                for coor2 in coor2s:
                                    dist = ((coor1[0] - coor2[0]) ** 2 + (coor1[1] - coor2[1]) ** 2 + (coor1[2] - coor2[2]) ** 2) ** 0.5
                                    dists.append(dist)
                            mindist = min(dists)
                            if mindist <= cutoff:
                                all_Mcontacts.append([chain1 + '.' + str(resid1), chain2 + '.' + str(resid2)])

    Mchain_pairs = set([])
    Mpair2contacts = {}
    for modcon in all_Mcontacts:
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
                    for Mcontact in Mcontacts:
                        if Mcontact[0] in Rchain2resids[Rchain1] and Mcontact[1] in Rchain2resids[Rchain2]:
                            good_Mcontacts.add(Mcontact)

                    if good_Mcontacts:
                        overlap = Rcontacts.intersection(good_Mcontacts)
                        if overlap:
                            qs_best = round(len(overlap) / max(len(Rcontacts), len(good_Mcontacts)), 4)
                        else:
                            qs_best = 0.0
                        map_results.append([Mchain1, Mchain2, qs_best])
                    else:
                        map_results.append([Mchain1, Mchain2, 0.0])
        all_results.append(['reference', Rchain1, Rchain2, map_results])

    for Mchain_pair in Mchain_pairs:
        Mcontacts = Mpair2contacts[Mchain_pair]
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

                    good_Mcontacts = set([])
                    for Mcontact in Mcontacts:
                        if Mcontact[0] in Rchain2resids[Rchain1] and Mcontact[1] in Rchain2resids[Rchain2]:
                            good_Mcontacts.add(Mcontact)

                    if Rcontacts and good_Mcontacts:
                        overlap = Rcontacts.intersection(good_Mcontacts)
                        if overlap:
                            qs_best = round(len(overlap) / max(len(Rcontacts), len(good_Mcontacts)), 4)
                        else:
                            qs_best = 0.0
                        map_results.append([Rchain1, Rchain2, qs_best])
                    else:
                        map_results.append([Rchain1, Rchain2, 0.0])
        all_results.append(['prediction', Mchain1, Mchain2, map_results])
    return model_name, all_results

def save_qs_best(reference_pdb_path, model_pdb_path, ost_json_path, model_name, output_dir, cutoff=10, ca_distance_prefilter=20):
    """
    Calculate and save QS_best scores for a single model.
    
    Computes QS_best scores for the model and writes results to output file.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        model_pdb_path: Full path to model PDB file
        ost_json_path: Full path to OST JSON file
        model_name: Model name to use in output
        output_dir: Path to output directory (will create output_dir/model_name/)
        cutoff: Distance cutoff in Angstroms for identifying contacts (default: 10)
        ca_distance_prefilter: CA distance pre-filter threshold (default: 20)
    """
    # Override the global constant with the passed parameter
    global CA_DISTANCE_PREFILTER
    CA_DISTANCE_PREFILTER = ca_distance_prefilter
    
    # Create output directory for this model
    model_output_dir = os.path.join(output_dir, model_name)
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)
    
    model_name, all_results = get_qs_best(reference_pdb_path, model_pdb_path, ost_json_path, model_name, cutoff)
    
    rp = open(os.path.join(model_output_dir, model_name + '.qs'), 'w')
    for result in all_results:
        category = result[0]
        chain1 = result[1]
        chain2 = result[2]
        for item in result[3]:
            rp.write(model_name + '\t' + category + '\t' + chain1 + ':' + chain2 + '\t' + item[0] + ':' + item[1] + '\t' + str(item[2]) + '\n')
    rp.close()





if __name__ == "__main__":
    input_dir = sys.argv[1]
    target = sys.argv[2]
    name = sys.argv[3]
    output_dir = sys.argv[4]
    save_qs_best(input_dir, target, name, output_dir)



