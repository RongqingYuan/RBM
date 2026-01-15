"""
Generate and execute commands for external structure analysis tools.

This module creates and executes commands for running DockQ, lDDT, and TMscore
on pairwise interface structures.
"""

import os
import sys


def save_inputs_and_run(model_name, output_dir, dockq_path, lddt_path, tmscore_path):
    """
    Generate and execute commands for running external analysis tools.
    
    Reads the list of pairwise interfaces and directly executes DockQ, lDDT, and
    TMscore commands. Only runs for interfaces that don't already have result files.
    
    Args:
        model_name: Model name
        output_dir: Path to output directory (will use output_dir/model_name/)
        dockq_path: Path to DockQ executable
        lddt_path: Path to lDDT executable
        tmscore_path: Path to TMscore executable
    """
    print(f'Running DockQ, lDDT, and TMscore for {model_name}...')
    
    model_output_dir = os.path.join(output_dir, model_name)
    pairwise_dir = os.path.join(model_output_dir, 'pairwise_interfaces')
    pairwise_lddt_dir = os.path.join(model_output_dir, 'pairwise_interfaces_for_lddt')
    
    fp = open(os.path.join(model_output_dir, 'pairwise_interfaces.list'), 'r')
    for line in fp:
        words = line.split()
        # model_dir = words[0]  # Not needed anymore
        pair1 = words[1]
        pair2 = words[2]
        
        # Run DockQ if results don't exist
        dockq_out = os.path.join(pairwise_dir, f'{pair1}_{pair2}.dockq')
        if not os.path.exists(dockq_out):
            print(f'  Running DockQ for {pair1}_{pair2}...')
            cmd = f'{dockq_path} --no_align --n_cpu 1 {os.path.join(pairwise_dir, pair1 + "_" + pair2 + "_model.pdb")} {os.path.join(pairwise_dir, pair1 + "_" + pair2 + "_target.pdb")} > {dockq_out} 2> {os.path.join(pairwise_dir, pair1 + "_" + pair2 + ".log")}'
            os.system(cmd)
        
        # Run lDDT if results don't exist
        lddt_out = os.path.join(pairwise_lddt_dir, f'{pair1}_{pair2}.lddt')
        if not os.path.exists(lddt_out):
            print(f'  Running lDDT for {pair1}_{pair2}...')
            cmd = f'{lddt_path} -x {os.path.join(pairwise_lddt_dir, pair1 + "_" + pair2 + "_model.pdb")} {os.path.join(pairwise_lddt_dir, pair1 + "_" + pair2 + "_target.pdb")} > {lddt_out}'
            os.system(cmd)
        
        # Run TMscore if results don't exist
        tm_out = os.path.join(pairwise_dir, f'{pair1}_{pair2}.tm')
        if not os.path.exists(tm_out):
            print(f'  Running TMscore for {pair1}_{pair2}...')
            cmd = f'{tmscore_path} -c -ter 0 {os.path.join(pairwise_dir, pair1 + "_" + pair2 + "_model.pdb")} {os.path.join(pairwise_dir, pair1 + "_" + pair2 + "_target.pdb")} > {tm_out}'
            os.system(cmd)
    fp.close()
    
    print(f'Completed all external tools for {model_name}')



# if __name__ == "__main__":
#     output_dir = sys.argv[1]
#     target = sys.argv[2]
#     dockq_path = sys.argv[3]
#     lddt_path = sys.argv[4]
#     tmscore_path = sys.argv[5]
#     n_cpu = sys.argv[6]
#     save_inputs(output_dir, target, dockq_path, lddt_path, tmscore_path, n_cpu)