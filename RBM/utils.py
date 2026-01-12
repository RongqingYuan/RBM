"""
Utility functions shared across the RBM module.

This module contains common helper functions used by multiple scripts
in the oligomer evaluation pipeline.
"""

import os
from typing import List, Tuple, Dict, Set, Optional


def get_models(input_dir: str, target: str, name: str) -> List[str]:
    """
    Read model names from the model list file.
    
    The model list file should contain model names, one per line.
    Lines starting with '#' are treated as comments and skipped.
    
    Args:
        input_dir: Path to input directory containing target folders
        target: Target name (folder name)
        name: Name for the model list file (usually same as target)
    
    Returns:
        List of model names to evaluate
    
    Example:
        >>> models = get_models('./input', 'H0208', 'H0208')
        >>> print(len(models))
        315
    """
    model_list_path = os.path.join(input_dir, target, f'{name}.txt')
    
    models = []
    start = False
    
    with open(model_list_path, 'r') as fp:
        for line in fp:
            words = line.split()
            if not words:
                continue
            
            # Skip comment lines starting with #
            if words[0] == '#':
                start = True
                continue
            
            # After seeing a comment, start reading model names
            if start and len(words) > 1:
                models.append(words[1])
    
    return models


def get_models_set(input_dir: str, target: str, name: str) -> Set[str]:
    """
    Read model names from the model list file and return as a set.
    
    Same as get_models() but returns a set for faster membership testing.
    
    Args:
        input_dir: Path to input directory
        target: Target name
        name: Name for the model list file
    
    Returns:
        Set of model names
    """
    return set(get_models(input_dir, target, name))


def get_best_score_and_match(
    pair2scores: Dict, 
    model: str, 
    pair: str
) -> Tuple[float, Optional[str]]:
    """
    Extract best score and corresponding match from score dictionary.
    
    Args:
        pair2scores: Dictionary mapping model->pair to list of [match, score]
        model: Model name
        pair: Interface pair (e.g., 'A:B')
    
    Returns:
        Tuple of (best_score, best_match) or (0.0, None) if not found
    
    Example:
        >>> scores = {'model1': {'A:B': [['C:D', 0.85], ['E:F', 0.72]]}}
        >>> best_score, best_match = get_best_score_and_match(scores, 'model1', 'A:B')
        >>> print(best_score, best_match)
        0.85 C:D
    """
    try:
        # Sort by score (second element) in descending order
        pair2scores[model][pair].sort(key=lambda x: x[1], reverse=True)
        best_score = pair2scores[model][pair][0][1]
        best_match = pair2scores[model][pair][0][0]
        return best_score, best_match
    except KeyError:
        return 0.0, None


def process_best_scores_for_metrics(
    results: List,
    Rpair2dockq: Dict,
    Rpair2lddt: Dict,
    Rpair2tm: Dict,
    Mpair2dockq: Dict,
    Mpair2lddt: Dict,
    Mpair2tm: Dict,
    model: str,
    pair: str,
    cate: str
) -> Tuple[float, float, float, float, float, float, Set[str], int]:
    """
    Process and extract best scores for all metrics for a single interface.
    
    This function consolidates the logic for finding the best IPS, ICS, QS,
    DockQ, lDDT, and TM-score across all possible matches for an interface.
    
    Args:
        results: List of [cate, pair, ips, ics, qs] for all matches
        Rpair2dockq: Reference interface to DockQ scores mapping
        Rpair2lddt: Reference interface to lDDT scores mapping
        Rpair2tm: Reference interface to TM-scores mapping
        Mpair2dockq: Model interface to DockQ scores mapping
        Mpair2lddt: Model interface to lDDT scores mapping
        Mpair2tm: Model interface to TM-scores mapping
        model: Model name
        pair: Interface pair identifier
        cate: Category ('reference' or 'prediction')
    
    Returns:
        Tuple containing:
        - best_ips: Best Interface Patch Similarity score
        - best_ics: Best Interface Contact Similarity score
        - best_qs: Best QS-score
        - best_dockq: Best DockQ score
        - best_lddt: Best lDDT score
        - best_tm: Best TM-score
        - best_matches: Set of all matches that contributed to any best score
        - check: Count of successfully found structural scores (0 or 3)
    """
    # Find best IPS, ICS, QS scores from results
    best_ips = 0.0
    best_ics = 0.0
    best_qs = 0.0
    best_matches = set()
    
    for item in results:
        # item format: [cate, pair, ips, ics, qs]
        if item[2] > best_ips:
            best_ips = item[2]
            best_matches.add(f"{item[0]}:{item[1]}")
        if item[3] > best_ics:
            best_ics = item[3]
            best_matches.add(f"{item[0]}:{item[1]}")
        if item[4] > best_qs:
            best_qs = item[4]
            best_matches.add(f"{item[0]}:{item[1]}")
    
    # Select appropriate dictionaries based on category
    if cate == 'reference':
        pair2dockq = Rpair2dockq
        pair2lddt = Rpair2lddt
        pair2tm = Rpair2tm
    else:  # prediction
        pair2dockq = Mpair2dockq
        pair2lddt = Mpair2lddt
        pair2tm = Mpair2tm
    
    # Get best structural scores
    check = 0
    
    # DockQ
    best_dockq, match_dockq = get_best_score_and_match(pair2dockq, model, pair)
    if match_dockq:
        best_matches.add(match_dockq)
        check += 1
    
    # lDDT
    best_lddt, match_lddt = get_best_score_and_match(pair2lddt, model, pair)
    if match_lddt:
        best_matches.add(match_lddt)
        check += 1
    
    # TM-score
    best_tm, match_tm = get_best_score_and_match(pair2tm, model, pair)
    if match_tm:
        best_matches.add(match_tm)
        check += 1
    
    return (best_ips, best_ics, best_qs, best_dockq, best_lddt, best_tm,
            best_matches, check)


def format_interface_score_line(
    model: str,
    category: str,
    pair: str,
    best_matches: Set[str],
    weight: float,
    best_ips: float,
    best_ics: float,
    best_qs: float,
    best_dockq: float,
    best_lddt: float,
    best_tm: float
) -> str:
    """
    Format a line for the interface scores output file.
    
    Args:
        model: Model name
        category: 'REF' or 'MOD'
        pair: Interface pair (e.g., 'A:B')
        best_matches: Set of matching interfaces
        weight: Interface weight
        best_ips: Best IPS score
        best_ics: Best ICS score
        best_qs: Best QS score
        best_dockq: Best DockQ score
        best_lddt: Best lDDT score
        best_tm: Best TM-score
    
    Returns:
        Formatted output line with tab-separated values
    """
    match_str = ','.join(best_matches) if best_matches else 'na'
    
    return (f"{model}\t{category}\t{pair}\t{match_str}\t{weight}\t"
            f"{best_ips}\t{best_ics}\t{best_qs}\t{best_dockq}\t"
            f"{best_lddt}\t{best_tm}\n")


def safe_file_open(filepath: str, mode: str = 'r'):
    """
    Safely open a file with better error handling.
    
    Args:
        filepath: Path to file
        mode: File open mode ('r', 'w', etc.)
    
    Returns:
        File handle
    
    Raises:
        FileNotFoundError: If file doesn't exist (for read mode)
        IOError: If file cannot be opened
    """
    try:
        return open(filepath, mode)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")
    except IOError as e:
        raise IOError(f"Cannot open file {filepath}: {e}")


def ensure_dir(directory: str):
    """
    Create directory if it doesn't exist.
    
    Args:
        directory: Path to directory
    """
    os.makedirs(directory, exist_ok=True)

