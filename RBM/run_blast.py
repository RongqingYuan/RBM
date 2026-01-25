import os
import argparse
from Bio import Align
from Bio.Seq import Seq

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


def run_blast_alignment(query_seq, target_seq):
    """Run pairwise alignment using BioPython."""
    # Create aligner with default parameters
    aligner = Align.PairwiseAligner()
    
    # Perform alignment
    alignments = aligner.align(query_seq, target_seq)
    
    if not alignments:
        return None
    
    # Get best alignment
    best_alignment = alignments[0]
    
    # Extract aligned sequences
    aligned_query = str(best_alignment[0])
    aligned_target = str(best_alignment[1])
    
    # Find start positions (first non-gap position)
    query_start = 0
    for i, char in enumerate(aligned_query):
        if char != '-':
            query_start = i
            break
    
    target_start = 0
    for i, char in enumerate(aligned_target):
        if char != '-':
            target_start = i
            break
    
    # Calculate identity excluding gaps
    matches = sum(1 for q, t in zip(aligned_query, aligned_target) 
                 if q != '-' and t != '-' and q == t)
    aligned_positions = sum(1 for q, t in zip(aligned_query, aligned_target) 
                           if q != '-' and t != '-')
    
    identity = matches / aligned_positions if aligned_positions > 0 else 0.0
    
    return {
        'identity': identity,
        'query_aligned': aligned_query,
        'target_aligned': aligned_target,
        'query_start': query_start,
        'target_start': target_start
    }


def blast_chain_mapping(target_pdb, model_pdb, model_chain, output_pdb):
    # Read PDB files (all lines for target, ATOM lines for model)
    with open(target_pdb) as f:
        target_lines = f.readlines()
    with open(model_pdb) as f:
        model_lines = [line for line in f if line.startswith('ATOM')]
    
    # Get sequences from ATOM lines
    target_atom_lines = [line for line in target_lines if line.startswith('ATOM')]
    target_seqs = get_chain_sequences(target_atom_lines)
    model_seqs = get_chain_sequences(model_lines)
    
    if model_chain not in model_seqs:
        raise ValueError(f"Model chain '{model_chain}' not found in {model_pdb}")
    
    # Convert model chain sequence dict to string
    model_residues = sorted(model_seqs[model_chain].items())
    model_seq_str = ''.join(res for _, res in model_residues)
    
    print(f"Model chain {model_chain}: {len(model_seq_str)} residues")
    print(f"Target chains: {list(target_seqs.keys())}")
    print()
    
    # BLAST against all target chains
    best_match = None
    best_identity = 0.0
    best_target_chain = None
    
    for target_chain, target_seq_dict in target_seqs.items():
        target_residues = sorted(target_seq_dict.items())
        target_seq_str = ''.join(res for _, res in target_residues)
        
        print(f"BLASTing against target chain {target_chain} ({len(target_seq_str)} residues)...")
        
        alignment = run_blast_alignment(model_seq_str, target_seq_str)
        
        if alignment:
            identity = alignment['identity']
            print(f"  Identity (excluding gaps): {identity:.2%}")
            
            if identity > best_identity:
                best_identity = identity
                best_match = alignment
                best_target_chain = target_chain
                best_match['target_chain'] = target_chain
                best_match['target_residues'] = target_residues
                best_match['model_residues'] = model_residues
    
    print()
    print("="*60)
    
    if best_match is None:
        print("ERROR: No BLAST alignment found!")
        print("Unable to align sequences. Please check the input PDB files.")
        raise ValueError(f"No BLAST alignment found for model chain '{model_chain}'")
    
    print(f"Best match: Target chain {best_target_chain}")
    print(f"Identity (excluding gaps): {best_identity:.2%}")
    print()
    
    if best_identity < 0.99:
        print(f"ERROR: Best identity ({best_identity:.2%}) is below 99% threshold!")
        print(f"Model chain '{model_chain}' does not match any target chain well enough.")
        print("sequences are too different for reliable chain mapping")
        print("Please check the input PDB files. Make sure the residues and numberings are correct.")
        raise ValueError(
            f"Model chain '{model_chain}' has only {best_identity:.2%} identity with "
            f"best target match (chain {best_target_chain}). Required: >= 99%"
        )
    
    # Create residue mapping
    print("Creating residue number mapping...")
    
    # Build mapping from alignment
    residue_mapping = {}
    query_idx = best_match['query_start']
    target_idx = best_match['target_start']
    
    for q_aa, t_aa in zip(best_match['query_aligned'], best_match['target_aligned']):
        if q_aa != '-' and t_aa != '-':
            # Both have residues - create mapping from target to model
            target_resnum = best_match['target_residues'][target_idx][0]
            model_resnum = best_match['model_residues'][query_idx][0]
            residue_mapping[target_resnum] = model_resnum
        
        if q_aa != '-':
            query_idx += 1
        if t_aa != '-':
            target_idx += 1
    
    # Renumber target PDB residues to match model
    print(f"Renumbering residues in target PDB to match model...")
    
    renumbered_lines = []
    for line in target_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21]
            if chain == best_target_chain:
                try:
                    old_resnum = int(line[22:26].strip())
                    if old_resnum in residue_mapping:
                        new_resnum = residue_mapping[old_resnum]
                        # Replace residue number (columns 23-26, right-justified)
                        new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                        renumbered_lines.append(new_line)
                    else:
                        # Residue not in mapping (gap in alignment), keep original
                        renumbered_lines.append(line)
                except ValueError:
                    # Can't parse residue number, keep original
                    renumbered_lines.append(line)
            else:
                # Different chain, keep original
                renumbered_lines.append(line)
        else:
            # Not ATOM/HETATM, keep original
            renumbered_lines.append(line)
    
    # Write renumbered PDB
    with open(output_pdb, 'w') as f:
        f.writelines(renumbered_lines)
    
    print(f"Renumbered target PDB saved to: {output_pdb}")
    print(f"  Target chain {best_target_chain}: {len(residue_mapping)} residues renumbered")
    print(f"  Renumbered to match model chain: {model_chain}")
    print("="*60)



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Run BLAST to find best chain mapping between target and model'
    )
    parser.add_argument('--target_pdb', type=str, required=True,
                       help='Path to target/reference PDB file')
    parser.add_argument('--model_pdb', type=str, required=True,
                       help='Path to model PDB file')
    parser.add_argument('--model_chain', type=str, required=True,
                       help='Model chain ID that needs to be mapped to target')
    parser.add_argument('--output_pdb', type=str, required=True,
                       help='Path to output PDB file with renumbered residues')
    
    args = parser.parse_args()
    
    # Verify input files exist
    if not os.path.exists(args.target_pdb):
        raise FileNotFoundError(f"Target PDB file not found: {args.target_pdb}")
    if not os.path.exists(args.model_pdb):
        raise FileNotFoundError(f"Model PDB file not found: {args.model_pdb}")
    
    print(f"Running BLAST chain mapping...")
    print(f"  Target PDB: {args.target_pdb}")
    print(f"  Model PDB: {args.model_pdb}")
    print(f"  Model chain: {args.model_chain}")
    print(f"  Output PDB: {args.output_pdb}")
    print()
    
    blast_chain_mapping(args.target_pdb, args.model_pdb, args.model_chain, args.output_pdb)