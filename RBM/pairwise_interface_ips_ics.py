"""
Calculate Interface Patch Similarity (IPS) and Interface Contact Similarity (ICS) scores.

This module computes pairwise interface quality scores by comparing reference and model
protein structures. IPS measures the overlap of interface residues, while ICS measures
the precision and recall of interface contacts using the F1 score.
"""

import sys
import json
from itertools import permutations
from multiprocessing import Pool
import os
from .utils import get_model_name_from_path

def get_Rchain2resids(reference_pdb_path):
    """
    Extract residue IDs for each chain from the reference (target) PDB file.
    
    Args:
        reference_pdb_path: Full path to the reference PDB file
    """
    fp = open(reference_pdb_path, 'r')
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

def get_ips_and_ics(reference_pdb_path, ost_json_path, model_name, output_dir):
    """
    Calculate IPS and ICS scores for a single model and save results.
    
    This function processes a model structure by comparing it with the reference structure.
    It loads chain mappings and contact information from OpenStructure (OST) JSON files,
    identifies interface pairs between reference and model chains, and computes quality scores.
    
    The function handles two types of interface pairs:
    - Reference interfaces: compares model chain pairs that match reference chain pairs
    - Prediction interfaces: compares reference chain pairs that match model chain pairs
    
    IPS (Interface Patch Similarity) is calculated as:
    RMcount / (Rcount + Mcount - RMcount)
    where RMcount is the number of overlapping residues, Rcount and Mcount are total
    residues in reference and model interfaces respectively.
    
    ICS (Interface Contact Similarity) is calculated as the F1 score:
    2 * precision * recall / (precision + recall)
    where precision = overlap / reference_contacts and recall = overlap / model_contacts.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        ost_json_path: Full path to OST JSON file with chain mappings and contacts
        model_name: Model name to use in output files (e.g., 'H0208TS014_1')
        output_dir: Path to output directory where results will be saved
        target: Target name for organizing output files (e.g., 'H0208')
    
    Returns:
        int: Always returns 0 on success
    
    Output Files:
        Creates two result files:
        - {output_dir}/{target}/IPS/{model_name}.result
        - {output_dir}/{target}/ICS/{model_name}.result
        
        Each file contains interface pairs with their corresponding scores.
        Format: '>category chain1:chain2 size1:size2' followed by match results.
    """
    Rchain2resids = get_Rchain2resids(reference_pdb_path)
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
    
    # Create output directory for this model
    model_output_dir = os.path.join(output_dir, model_name)
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)
    
    # Write IPS results in unified format (matching QS format with interface sizes)
    rp = open(os.path.join(model_output_dir, model_name + '.ips'), 'w')
    for result in all_results:
        category = result[0]
        chain1 = result[1]
        chain2 = result[2]
        size1 = str(result[3])
        size2 = str(result[4])
        for item in result[5]:
            # Format: model_name category chainpair matchpair size1 size2 score
            rp.write(model_name + '\t' + category + '\t' + chain1 + ':' + chain2 + '\t' + 
                    item[0] + ':' + item[1] + '\t' + size1 + '\t' + size2 + '\t' + str(item[4]) + '\n')
    rp.close()
    
    # Write ICS results in unified format (matching QS format with interface sizes)
    rp = open(os.path.join(model_output_dir, model_name + '.ics'), 'w')
    for result in all_results:
        category = result[0]
        chain1 = result[1]
        chain2 = result[2]
        size1 = str(result[3])
        size2 = str(result[4])
        for item in result[5]:
            # Format: model_name category chainpair matchpair size1 size2 score
            rp.write(model_name + '\t' + category + '\t' + chain1 + ':' + chain2 + '\t' + 
                    item[0] + ':' + item[1] + '\t' + size1 + '\t' + size2 + '\t' + str(item[5]) + '\n')
    rp.close()
    return 0

def save_ips_and_ics(reference_pdb_path, ost_json_path, model_name, output_dir):
    """
    Calculate and save IPS and ICS scores for a single model.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        ost_json_path: Full path to OST JSON file
        model_name: Model name to use in output files
        output_dir: Path to output directory (will create output_dir/model_name/)
    """
    get_ips_and_ics(reference_pdb_path, ost_json_path, model_name, output_dir)

# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     save_ips_and_ics(input_dir, target, name, output_dir)


