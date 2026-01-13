"""Utility functions shared across the RBM module."""

import os
from typing import List


def get_models(input_dir: str, target: str, name: str) -> List[str]:
    """
    Read model names from the model list file.
    
    Parses the model list file and extracts model names. Lines starting with '#'
    are treated as comments. Model names are read after the first comment line.
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
