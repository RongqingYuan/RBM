"""
Utility functions shared across the RBM module.

This module contains common helper functions used by multiple scripts
in the oligomer evaluation pipeline.
"""

import os
from typing import List


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
