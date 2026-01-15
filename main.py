import argparse
import os
from RBM.pairwise_interface_ips_ics import save_ips_and_ics
from RBM.pairwise_interface_qs import save_qs_best
from RBM.extract_pairwise_interfaces import save_pairwise_interfaces
from RBM.commands import save_inputs_and_run
from RBM.pairwise_dockq import get_dockq_scores
from RBM.pairwise_lddt import get_lddt_scores
from RBM.pairwise_tmscore import get_tm_scores
from RBM.interface_scores import save_interface_scores
from RBM.model_scores import save_model_scores_v1, save_model_scores_v2, save_model_scores_v3
from RBM.antibody_scores import save_antibody_scores_v3


# ============================================================================
# Internal Constants (Do NOT modify unless you know what you're doing)
# ============================================================================

# Score filter ratio: only keep interface matches where at least one score
# is greater than (best_score * SCORE_FILTER_RATIO)
# This is an internal threshold that should not be changed in normal usage
SCORE_FILTER_RATIO = 0.5

# The cutoff for pre-filtering contacts based on CA distance
CA_DISTANCE_PREFILTER = 20.0

# The cutoff for QS-score
QS_CUTOFF = 10

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='RBM (Reciprocal Best Matching) - Evaluate protein complex structure predictions'
    )
    
    # File-based arguments (new approach)
    parser.add_argument('--reference_pdb', type=str, required=True,
                       help='Path to reference/target PDB file')
    parser.add_argument('--model_pdb', type=str, required=True,
                       help='Path to model PDB file')
    parser.add_argument('--ost_json', type=str, required=True,
                       help='Path to OpenStructure (OST) JSON file with chain mappings and contacts')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Path to output directory')
    parser.add_argument('--target_name', type=str, required=True,
                       help='Target name for organizing output files (e.g., H0208)')
    parser.add_argument('--model_name', type=str, required=False,
                       help='Model name to use in output files (default: basename of model_pdb)')
    
    # Tool paths
    parser.add_argument('--dockq_path', type=str, default="DockQ",
                       help='Path to DockQ executable (default: DockQ)')
    parser.add_argument('--lddt_path', type=str, default="lddt",
                       help='Path to lDDT executable (default: lddt)')
    parser.add_argument('--tmscore_path', type=str, default="TMscore",
                       help='Path to TMscore executable (default: TMscore)')
    
    # Optional arguments
    parser.add_argument('--scores', nargs='+', choices=['ICS', 'IPS', 'QS_best', 'DockQ', 'lDDT', 'TMscore'], 
                       default=['ICS', 'IPS', 'QS_best', 'DockQ', 'lDDT', 'TMscore'],
                       help='Choose 1-6 scores from: ICS, IPS, QS_best, DockQ, lDDT, TMscore')
    parser.add_argument('--antibody', action='store_true',
                       help='Enable antibody-specific analysis')
    parser.add_argument('--chainAs', type=str, default="",
                       help='Chain IDs for antibody component (for antibody mode)')
    parser.add_argument('--chainBs', type=str, default="",
                       help='Chain IDs for antigen component (for antibody mode)')
    parser.add_argument('--rbm_version', type=str, choices=['min', 'all', 'average'], default='min',
                       help='RBM version: min, all, or average (default: min)')
    parser.add_argument('--interface_weight', type=str, choices=['log2', 'log10', 'linear'], default='log10',
                       help='Weighting method: log2, log10, or linear (default: log10)')
    parser.add_argument('--keep_tmp', action='store_true',
                       help='Keep temporary interface files in interface_tmp directory (default: False, will remove after completion)')

    args = parser.parse_args()
    
    # Extract arguments
    reference_pdb_path = args.reference_pdb
    model_pdb_path = args.model_pdb
    ost_json_path = args.ost_json
    output_dir = args.output_dir
    target_name = args.target_name
    model_name = args.model_name if args.model_name else os.path.basename(model_pdb_path)
    dockq_path = args.dockq_path
    lddt_path = args.lddt_path
    tmscore_path = args.tmscore_path
    antibody = args.antibody
    chainAs = args.chainAs
    chainBs = args.chainBs
    rbm_version = args.rbm_version
    interface_weight = args.interface_weight
    scores = args.scores
    keep_tmp = args.keep_tmp
    
    # Verify input files exist
    if not os.path.exists(reference_pdb_path):
        raise FileNotFoundError(f"Reference PDB file not found: {reference_pdb_path}")
    if not os.path.exists(model_pdb_path):
        raise FileNotFoundError(f"Model PDB file not found: {model_pdb_path}")
    if not os.path.exists(ost_json_path):
        raise FileNotFoundError(f"OST JSON file not found: {ost_json_path}")
    
    print(f"Starting RBM analysis...")
    print(f"  Reference: {reference_pdb_path}")
    print(f"  Model: {model_pdb_path}")
    print(f"  OST JSON: {ost_json_path}")
    print(f"  Target name: {target_name}")
    print(f"  Model name: {model_name}")
    print(f"  Output directory: {output_dir}")
    print()
    
    # Run the pipeline
    print("Step 1: Calculating IPS and ICS scores...")
    save_ips_and_ics(reference_pdb_path, ost_json_path, model_name, output_dir)
    
    print("Step 2: Calculating QS_best scores...")
    save_qs_best(reference_pdb_path, model_pdb_path, ost_json_path, model_name, output_dir, QS_CUTOFF, CA_DISTANCE_PREFILTER)
    
    print("Step 3: Extracting pairwise interfaces...")
    save_pairwise_interfaces(reference_pdb_path, model_pdb_path, model_name, output_dir, SCORE_FILTER_RATIO)
    
    print("Step 4: Running external tools (DockQ, lDDT, TMscore)...")
    save_inputs_and_run(model_name, output_dir, dockq_path, lddt_path, tmscore_path)
    
    print("Step 5: Collecting DockQ scores...")
    get_dockq_scores(model_name, output_dir)
    
    print("Step 6: Collecting lDDT scores...")
    get_lddt_scores(model_name, output_dir)
    
    print("Step 7: Collecting TMscore scores...")
    get_tm_scores(model_name, output_dir)
    
    print("Step 8: Aggregating interface scores...")
    save_interface_scores(model_name, output_dir)
    
    if not antibody:
        print(f"Step 9: Calculating model scores (RBM version: {rbm_version}, interface weight: {interface_weight})...")
        if rbm_version == 'average':
            save_model_scores_v1(model_name, output_dir, interface_weight)
        elif rbm_version == 'all':
            save_model_scores_v2(model_name, output_dir, interface_weight)
        elif rbm_version == 'min':
            save_model_scores_v3(model_name, output_dir, interface_weight)
        else:
            raise ValueError("Invalid RBM version")
        print()
        print("="*60)
        print("RBM scoring completed successfully!")
        print(f"  RBM version: {rbm_version}")
        print(f"  Interface weight: {interface_weight}")
        print(f"  Results saved to: {output_dir}/{model_name}/")
        print("="*60)
    else:
        if chainAs and chainBs:
            print("Step 9: Calculating antibody-specific scores...")
            chainA_list = list(chainAs)
            chainB_list = list(chainBs)
            save_antibody_scores_v3(model_name, output_dir, chainA_list, chainB_list)
            print()
            print("="*60)
            print("Antibody scoring completed successfully!")
            print(f"  Antibody chains: {chainAs}")
            print(f"  Antigen chains: {chainBs}")
            print(f"  Results saved to: {output_dir}/{model_name}/")
            print("="*60)
        else:
            raise ValueError("Please provide valid chainAs and chainBs for antibody mode")
    
    # Cleanup temporary interface files if requested
    if not keep_tmp:
        import shutil
        interface_tmp_dir = os.path.join(output_dir, model_name, 'interface_tmp')
        if os.path.exists(interface_tmp_dir):
            print()
            print("Cleaning up temporary interface files...")
            shutil.rmtree(interface_tmp_dir)
            print(f"  Removed: {interface_tmp_dir}")