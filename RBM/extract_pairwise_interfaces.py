"""
Extract pairwise interface structures from reference and model PDB files.

This module processes interface scores (IPS, ICS, QS_best) to identify high-quality
interface pairs, then extracts the corresponding interface regions from reference
and model structures into separate PDB files for further analysis.
"""

import os
import sys
import json
from multiprocessing import Pool
from .utils import get_model_name_from_path, parse_reference_pdb


# ============================================================================
# Configuration Constants
# ============================================================================

# Score filter ratio: only keep interface matches where at least one score
# is greater than (best_score * SCORE_FILTER_RATIO)
SCORE_FILTER_RATIO = 0.5


def get_model2qsbest(model_name, output_dir):
    """
    Load QS_best scores from result file.
    
    Reads the QS_best result file and organizes scores by category and interface pair.
    
    Args:
        model_name: Model name
        output_dir: Path to output directory
    """
    qs_file = os.path.join(output_dir, model_name, model_name + '.qs')
    model2qsbest = {model_name: {}}
    
    with open(qs_file, 'r') as fp:
        for line in fp:
            words = line.split()
            model = words[0]
            cate = words[1]
            pair1 = words[2].replace(':','-')
            pair2 = words[3].replace(':','-')
            model2qsbest[model][(cate, pair1, pair2)] = float(words[4])
    
    return model2qsbest


def get_pair2scores(model_name, model2qsbest, output_dir):
    """
    Load IPS and ICS scores for a model and combine with QS_best scores.
    
    Reads both IPS and ICS result files for a model and merges them with
    QS_best scores. Returns a dictionary mapping interface pairs to their
    combined score results.
    """
    IPS_file = os.path.join(output_dir, model_name, model_name + '.ips')
    ICS_file = os.path.join(output_dir, model_name, model_name + '.ics')
    pair2results = {}
    
    # Read IPS scores - new unified format: model_name category chainpair matchpair size1 size2 score
    ips_data = {}
    with open(IPS_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # parts: [model_name, category, chainpair, matchpair, size1, size2, score]
            category = parts[1]
            chainpair = parts[2].replace(':', '-')
            matchpair = parts[3].replace(':', '-')
            score = float(parts[6])
            key = (category, chainpair, matchpair)
            ips_data[key] = score
    
    # Read ICS scores - new unified format: model_name category chainpair matchpair size1 size2 score
    ics_data = {}
    with open(ICS_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            category = parts[1]
            chainpair = parts[2].replace(':', '-')
            matchpair = parts[3].replace(':', '-')
            score = float(parts[6])
            key = (category, chainpair, matchpair)
            ics_data[key] = score
    
    # Combine IPS and ICS data
    for key in ips_data:
        category, chainpair, matchpair = key
        ips_score = ips_data[key]
        ics_score = ics_data.get(key, 0.0)
        
        # Get QS_best score
        qsbest = model2qsbest[model_name][(category, chainpair, matchpair)]
        
        # Store in pair2results
        if (category, chainpair) not in pair2results:
            pair2results[(category, chainpair)] = []
        
        match_chain1, match_chain2 = matchpair.split('-')
        pair2results[(category, chainpair)].append([match_chain1, match_chain2, ips_score, ics_score, qsbest])
    
    return pair2results


def process_model(reference_pdb_path, model_pdb_path, model_name, output_dir):
    """
    Extract pairwise interface structures for a model.
    
    Identifies high-quality interface pairs based on IPS, ICS, and QS_best scores
    using a score filter ratio. Extracts the interface regions from reference and
    model structures and writes them to separate PDB files. Only extracts interfaces
    that don't already have complete result files (DockQ, lDDT, TMscore).
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        model_pdb_path: Full path to model PDB file
        model_name: Model name to use in output
        output_dir: Path to output directory (will use output_dir/model_name/)
    """
    Rchain2resids, _, Rchain2lines = parse_reference_pdb(reference_pdb_path)
    model2qsbest = get_model2qsbest(model_name, output_dir)
    pair2results = get_pair2scores(model_name, model2qsbest, output_dir)
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
                    cases.add((cate, pair1, item[0] + '-' + item[1]))
                elif cate == 'prediction':
                    if item[0] < item[1]:
                        cases.add((cate, item[0] + '-' + item[1], chain1A + '-' + chain1B))
                    else:
                        cases.add((cate, item[1] + '-' + item[0], chain1B + '-' + chain1A))

    fp = open(model_pdb_path, 'r')
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
    
    # Create pairwise_interfaces directory in interface_tmp
    interface_tmp_dir = os.path.join(output_dir, model_name, 'interface_tmp')
    pairwise_dir = os.path.join(interface_tmp_dir, 'pairwise_interfaces')
    os.makedirs(pairwise_dir, exist_ok=True)
    
    need_cases = []
    for case in cases:
        category = case[0]
        pair1 = case[1]
        pair2 = case[2]
        
        check1 = 0
        check2 = 0
        check3 = 0
        check4 = 0
        check5 = 0
        check6 = 0

        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_target.pdb')):
            check1 = 1
        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_model.pdb')):
            check2 = 1
        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.dockq')):
            check3 = 1
        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.lddt')):
            check4 = 1
        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.log')):
            check5 = 1
        if os.path.exists(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.tm')):
            check6 = 1

        if check1 and check2 and check3 and check4 and check5 and check6:
            pass
        else:
            if check1:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_target.pdb'))
            if check2:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_model.pdb'))
            if check3:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.dockq'))
            if check4:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.lddt'))
            if check5:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.log'))
            if check6:
                os.remove(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '.tm'))

            rp1 = open(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_target.pdb'), 'w')
            rp2 = open(os.path.join(pairwise_dir, pair1 + '_' + pair2 + '_model.pdb'), 'w')
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
    return model_name, need_cases



def process_model_same_chain(reference_pdb_path, model_pdb_path, model_name, output_dir):
    """
    Extract pairwise interface structures with same chain ID for lDDT analysis.
    
    Similar to process_model but writes interfaces with both chains assigned to
    chain 'A' and applies residue number shifting to the second chain. This format
    is required for lDDT calculation which expects single-chain structures.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        model_pdb_path: Full path to model PDB file
        model_name: Model name to use in output
        output_dir: Path to output directory (will use output_dir/model_name/)
    """
    Rchain2resids, _, Rchain2lines = parse_reference_pdb(reference_pdb_path)
    model2qsbest = get_model2qsbest(model_name, output_dir)
    pair2results = get_pair2scores(model_name, model2qsbest, output_dir)
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
                    cases.add((cate, pair1, item[0] + '-' + item[1]))
                elif cate == 'prediction':
                    if item[0] < item[1]:
                        cases.add((cate, item[0] + '-' + item[1], chain1A + '-' + chain1B))
                    else:
                        cases.add((cate, item[1] + '-' + item[0], chain1B + '-' + chain1A))

    fp = open(model_pdb_path, 'r')
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

    # Create pairwise_interfaces_for_lddt directory in interface_tmp
    interface_tmp_dir = os.path.join(output_dir, model_name, 'interface_tmp')
    pairwise_lddt_dir = os.path.join(interface_tmp_dir, 'pairwise_interfaces_for_lddt')
    os.makedirs(pairwise_lddt_dir, exist_ok=True)
    
    need_cases = []
    for case in cases:
        category = case[0]
        pair1 = case[1]
        pair2 = case[2]
        rp1 = open(os.path.join(pairwise_lddt_dir, pair1 + '_' + pair2 + '_target.pdb'), 'w')
        rp2 = open(os.path.join(pairwise_lddt_dir, pair1 + '_' + pair2 + '_model.pdb'), 'w')
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
    return model_name, need_cases

def save_pairwise_interfaces(reference_pdb_path, model_pdb_path, model_name, output_dir, score_filter_ratio):
    """
    Extract and save pairwise interface structures for a single model.
    
    Extracts interface pairs that meet the score filter criteria. Creates two sets
    of interface PDB files: one for general analysis (with separate chain IDs) and
    one for lDDT analysis (with same chain ID). Generates list files containing all
    extracted interface pairs.
    
    Args:
        reference_pdb_path: Full path to reference PDB file
        model_pdb_path: Full path to model PDB file
        model_name: Model name to use in output
        output_dir: Path to output directory (will use output_dir/model_name/)
        score_filter_ratio: Score filter ratio threshold
    """
    # Override the global constant with the passed parameter
    global SCORE_FILTER_RATIO
    SCORE_FILTER_RATIO = score_filter_ratio
    
    # Create model output directory and interface_tmp directory
    model_output_dir = os.path.join(output_dir, model_name)
    interface_tmp_dir = os.path.join(model_output_dir, 'interface_tmp')
    os.makedirs(interface_tmp_dir, exist_ok=True)
    
    # Process and save pairwise interfaces
    with open(os.path.join(interface_tmp_dir, 'pairwise_interfaces.list'), 'w') as f:
        model_name, need_cases = process_model(reference_pdb_path, model_pdb_path, model_name, output_dir)
        for case in need_cases:
            # case is (category, pair1, pair2)
            f.write(model_name + '\t' + case[0] + '\t' + case[1] + '\t' + case[2] + '\n')
    
    # Process and save pairwise interfaces for lDDT
    with open(os.path.join(interface_tmp_dir, 'pairwise_interfaces_for_lddt.list'), 'w') as f:
        model_name, need_cases = process_model_same_chain(reference_pdb_path, model_pdb_path, model_name, output_dir)
        for case in need_cases:
            # case is (category, pair1, pair2)
            f.write(model_name + '\t' + case[0] + '\t' + case[1] + '\t' + case[2] + '\n')


# if __name__ == "__main__":
#     input_dir = sys.argv[1]
#     target = sys.argv[2]
#     name = sys.argv[3]
#     output_dir = sys.argv[4]
#     save_pairwise_interfaces(input_dir, target, name, output_dir)


