"""Utility functions shared across the RBM module."""

import os
from typing import List, Tuple, Dict, Set


# CA distance pre-filter for computational optimization (Angstroms)
CA_DISTANCE_PREFILTER = 20.0


def get_model_name_from_path(model_path: str) -> str:
    """
    Extract model name from the full file path.
    
    Args:
        model_path: Full path to the model file
        
    Returns:
        The base name of the model file without directory path
    """
    return os.path.basename(model_path)


def parse_reference_pdb(reference_pdb_path: str, cutoff: float = 10.0) -> Tuple[Dict[str, Set[int]], list, Dict[str, list]]:
    """
    Parse reference PDB file and extract all chain-level information.
    
    This comprehensive function parses the reference PDB file once and returns:
    - Residue IDs for each chain
    - Inter-chain contacts based on distance cutoff
    - Original PDB lines for each chain
    
    Different modules can use what they need from the returned values.
    
    Args:
        reference_pdb_path: Full path to the reference PDB file
        cutoff: Distance cutoff in Angstroms for identifying contacts (default: 10.0)
    
    Returns:
        Tuple of (Rchain2resids, all_Rcontacts, Rchain2lines):
        - Rchain2resids: Dict mapping chain ID to set of residue IDs
        - all_Rcontacts: List of inter-chain contacts as [chain1.resid1, chain2.resid2]
        - Rchain2lines: Dict mapping chain ID to list of PDB ATOM lines
    """
    fp = open(reference_pdb_path, 'r')
    Rchains = []
    Rresid2CAcoor = {}
    Rresid2coors = {}
    all_Rcontacts = []
    Rchain2resids = {}
    Rchain2lines = {}
    
    # First pass: collect all chain information
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
                
                # Initialize chain data structures
                try:
                    Rchain2resids[chain].add(resid)
                    Rchain2lines[chain].append(line)
                    Rresid2coors[chain]
                    Rresid2CAcoor[chain]
                except KeyError:
                    Rchains.append(chain)
                    Rchain2resids[chain] = set([resid])
                    Rchain2lines[chain] = [line]
                    Rresid2coors[chain] = {}
                    Rresid2CAcoor[chain] = {}
                
                # Store coordinates for contact calculation
                if atom == 'CA':
                    Rresid2CAcoor[chain][resid] = [x, y, z]
                if element != 'H':
                    try:
                        Rresid2coors[chain][resid].append([x, y, z])
                    except KeyError:
                        Rresid2coors[chain][resid] = [[x, y, z]]
    fp.close()
    
    # Second pass: calculate inter-chain contacts
    Rchains.sort()
    for c1, chain1 in enumerate(Rchains):
        for c2, chain2 in enumerate(Rchains):
            if c1 < c2:
                for resid1 in Rchain2resids[chain1]:
                    for resid2 in Rchain2resids[chain2]:
                        CAcoor1 = Rresid2CAcoor[chain1][resid1]
                        CAcoor2 = Rresid2CAcoor[chain2][resid2]
                        CAdist = ((CAcoor1[0] - CAcoor2[0]) ** 2 + (CAcoor1[1] - CAcoor2[1]) ** 2 + (CAcoor1[2] - CAcoor2[2]) ** 2) ** 0.5
                        if CAdist < CA_DISTANCE_PREFILTER:
                            coor1s = Rresid2coors[chain1][resid1]
                            coor2s = Rresid2coors[chain2][resid2]
                            dists = []
                            for coor1 in coor1s:
                                for coor2 in coor2s:
                                    dist = ((coor1[0] - coor2[0]) ** 2 + (coor1[1] - coor2[1]) ** 2 + (coor1[2] - coor2[2]) ** 2) ** 0.5
                                    dists.append(dist)
                            mindist = min(dists)
                            if mindist <= cutoff:
                                all_Rcontacts.append([chain1 + '.' + str(resid1), chain2 + '.' + str(resid2)])
    
    return Rchain2resids, all_Rcontacts, Rchain2lines
