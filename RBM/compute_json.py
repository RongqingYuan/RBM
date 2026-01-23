"""
Generate OST-compatible JSON from PDB files.

This module creates chain mappings and contacts between target and model structures
in a format compatible with the RBM pipeline (similar to OpenStructure output).
"""

import json
import math

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def get_chain_sequences(pdb_lines):
    sequences = {}
    for line in pdb_lines:
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            chain = line[21]
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()
            if resname in THREE_TO_ONE:
                if chain not in sequences:
                    sequences[chain] = {}
                if resnum not in sequences[chain]:
                    sequences[chain][resnum] = THREE_TO_ONE[resname]
    return sequences


def sequence_similarity(seq1, seq2):
    """Calculate similarity score between two residue dictionaries."""
    if not seq1 or not seq2:
        return 0.0
    
    common = set(seq1.keys()) & set(seq2.keys())
    if not common:
        return 0.0
    
    # Identity in overlap
    identity = sum(1 for r in common if seq1[r] == seq2[r]) / len(common)
    # coverage = len(common) / max(len(seq1), len(seq2))
    # length_sim = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))
    
    return identity


def identical_in_overlap(seq1, seq2, threshold=0.99):
    common = set(seq1.keys()) & set(seq2.keys())
    if not common:
        return False
    identity = sum(1 for r in common if seq1[r] == seq2[r]) / len(common)
    return identity >= threshold


def cluster_chains(sequences):
    chains = sorted(sequences.keys())
    clusters = []
    assigned = set()
    
    for chain in chains:
        if chain in assigned:
            continue
        cluster = [chain]
        assigned.add(chain)
        for other in chains:
            if other not in assigned and identical_in_overlap(sequences[chain], sequences[other]):
                cluster.append(other)
                assigned.add(other)
        clusters.append(cluster)
    
    return clusters


def map_chains(target_seqs, model_seqs):
    """Map model chains to target chain clusters based on sequence similarity."""
    # Cluster target chains
    target_clusters = cluster_chains(target_seqs)
    
    # Find best target for each model chain
    model_to_target = {}
    for m_chain, m_seq in model_seqs.items():
        best_target, best_score = None, 0
        for t_chain, t_seq in target_seqs.items():
            score = sequence_similarity(m_seq, t_seq)
            if score > best_score:
                best_score, best_target = score, t_chain
        
        # Check if best score meets 99% threshold
        if best_score < 0.99:
            raise ValueError(
                f"Model chain '{m_chain}' cannot be mapped to any target chain with >= 99% "
                f"sequence identity (best match: {best_score:.2%}). "
                f"Please use run_blast.py to verify chain mapping."
            )
        
        model_to_target[m_chain] = best_target
    
    # Map clusters to their model chains
    chain_to_cluster = {c: tuple(cluster) for cluster in target_clusters for c in cluster}
    cluster_to_models = {}
    for m_chain, t_chain in model_to_target.items():
        cluster = chain_to_cluster[t_chain]
        cluster_to_models.setdefault(cluster, []).append(m_chain)
    
    # Return sorted groups
    target_groups = [list(c) for c in sorted(cluster_to_models.keys())]
    model_groups = [cluster_to_models[tuple(tg)] for tg in target_groups]
    
    return target_groups, model_groups


def get_contacts(pdb_lines, cutoff=5.0):
    """Find inter-chain contacts within distance cutoff (heavy atoms only)."""
    CA_DISTANCE_PREFILTER = 25.0
    
    # Parse atoms - separate CA and heavy atoms
    ca_coords = {}  # (chain, resnum) -> [x, y, z]
    heavy_atoms = {}  # (chain, resnum) -> [[x, y, z], ...]
    chains = []
    chain_residues = {}  # chain -> set of resnum
    
    for line in pdb_lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            chain = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            # Track chains
            if chain not in chain_residues:
                chains.append(chain)
                chain_residues[chain] = set()
            chain_residues[chain].add(resnum)
            
            # Store CA coordinates
            if atom_name == 'CA':
                ca_coords[(chain, resnum)] = [x, y, z]
            
            # Store heavy atom coordinates
            if not atom_name.startswith('H'):
                key = (chain, resnum)
                if key not in heavy_atoms:
                    heavy_atoms[key] = []
                heavy_atoms[key].append([x, y, z])
    
    # Find inter-chain contacts with CA prefilter
    contacts = set()
    chains.sort()
    
    for i, chain1 in enumerate(chains):
        for chain2 in chains[i+1:]:
            # Only check different chains
            for resnum1 in chain_residues[chain1]:
                for resnum2 in chain_residues[chain2]:
                    # Prefilter using CA distance
                    ca1 = ca_coords.get((chain1, resnum1))
                    ca2 = ca_coords.get((chain2, resnum2))
                    
                    if ca1 and ca2:
                        ca_dist = math.sqrt(
                            (ca1[0] - ca2[0])**2 + 
                            (ca1[1] - ca2[1])**2 + 
                            (ca1[2] - ca2[2])**2
                        )
                        
                        if ca_dist < CA_DISTANCE_PREFILTER:
                            # Check actual heavy atom distances
                            atoms1 = heavy_atoms.get((chain1, resnum1), [])
                            atoms2 = heavy_atoms.get((chain2, resnum2), [])
                            
                            min_dist = float('inf')
                            for coord1 in atoms1:
                                for coord2 in atoms2:
                                    dist = math.sqrt(
                                        (coord1[0] - coord2[0])**2 + 
                                        (coord1[1] - coord2[1])**2 + 
                                        (coord1[2] - coord2[2])**2
                                    )
                                    min_dist = min(min_dist, dist)
                            
                            if min_dist <= cutoff:
                                contact = [f"{chain1}.{resnum1}.", f"{chain2}.{resnum2}."]
                                contacts.add(tuple(sorted(contact)))
    
    return [list(c) for c in sorted(contacts)]


def generate_json(target_pdb, model_pdb, output_json=None, cutoff=5.0):
    # Read ATOM lines
    with open(target_pdb) as f:
        target_lines = [line for line in f if line.startswith('ATOM')]
    with open(model_pdb) as f:
        model_lines = [line for line in f if line.startswith('ATOM')]
    
    # Get chain mapping and contacts
    target_seqs = get_chain_sequences(target_lines)
    model_seqs = get_chain_sequences(model_lines)
    target_groups, model_groups = map_chains(target_seqs, model_seqs)
    target_contacts = get_contacts(target_lines, cutoff)
    model_contacts = get_contacts(model_lines, cutoff)
    
    result = {
        'chem_groups': target_groups,
        'chem_mapping': model_groups,
        'reference_contacts': target_contacts,
        'model_contacts': model_contacts
    }
    
    if output_json:
        with open(output_json, 'w') as f:
            json.dump(result, f, indent=2)
    
    return result


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate OST-compatible JSON from target and model PDB files'
    )
    parser.add_argument('--target_pdb', type=str, required=True,
                       help='Path to target/reference PDB file')
    parser.add_argument('--model_pdb', type=str, required=True,
                       help='Path to model PDB file')
    parser.add_argument('--output_json', type=str, required=True,
                       help='Path to output JSON file')
    parser.add_argument('--contact_cutoff', type=float, default=5.0,
                       help='Distance cutoff for contact detection in Angstroms (default: 5.0)')
    
    args = parser.parse_args()
    
    # Verify input files exist
    import os
    if not os.path.exists(args.target_pdb):
        raise FileNotFoundError(f"Target PDB file not found: {args.target_pdb}")
    if not os.path.exists(args.model_pdb):
        raise FileNotFoundError(f"Model PDB file not found: {args.model_pdb}")
    
    print(f"Generating OST-compatible JSON...")
    print(f"  Target PDB: {args.target_pdb}")
    print(f"  Model PDB: {args.model_pdb}")
    print(f"  Output JSON: {args.output_json}")
    print(f"  Contact cutoff: {args.contact_cutoff} Å")
    print()
    
    # Generate JSON with custom cutoff if provided
    result = generate_json(args.target_pdb, args.model_pdb, args.output_json, args.contact_cutoff)
    
    print("="*60)
    print("JSON generation completed successfully!")
    print(f"  Target groups: {result['chem_groups']}")
    print(f"  Model groups: {result['chem_mapping']}")
    print(f"  Target contacts: {len(result['reference_contacts'])}")
    print(f"  Model contacts: {len(result['model_contacts'])}")
    print(f"  Saved to: {args.output_json}")
    print("="*60)

