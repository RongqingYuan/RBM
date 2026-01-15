"""Utility functions shared across the RBM module."""

import os
from typing import List


def get_model_name_from_path(model_path: str) -> str:
    """
    Extract model name from the full file path.
    
    Args:
        model_path: Full path to the model file
        
    Returns:
        The base name of the model file without directory path
    """
    return os.path.basename(model_path)
