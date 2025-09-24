import argparse
from RBM.pairwise_interface_ips_ics import save_ips_and_ics
from RBM.pairwise_interface_qs import save_qs_best
from RBM.extract_pairwise_interfaces import save_pairwise_interfaces
from RBM.commands import save_inputs, run_dockq, run_lddt, run_tmscore
from RBM.pairwise_dockq import get_dockq_scores
from RBM.pairwise_lddt import get_lddt_scores
from RBM.pairwise_tmscore import get_tm_scores
from RBM.interface_scores import save_interface_scores
from RBM.model_scores import save_model_scores_v1, save_model_scores_v2, save_model_scores_v3
from RBM.antibody_scores import save_antibody_scores_v3


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--target', type=str, required=True)
    parser.add_argument('--name', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--qs_cutoff', type=int, default=10)
    parser.add_argument('--dockq_path', type=str, default="DockQ")
    parser.add_argument('--lddt_path', type=str, default="lddt")
    parser.add_argument('--tmscore_path', type=str, default="TMscore")
    parser.add_argument('--n_cpu', type=int, default=48)
    parser.add_argument('--antibody', action='store_true')
    parser.add_argument('--chainAs', type=str, default="")
    parser.add_argument('--chainBs', type=str, default="")
    parser.add_argument('--rbm_version', type=str, choices=['min', 'all', 'average'], default='min', help='RBM version: min, all, or average')
    parser.add_argument('--interface_weight', type=str, choices=['log2', 'log10', 'linear'], default='log10', help='Weighting method: log2, log10, or linear')

    args = parser.parse_args()
    input_dir = args.input_dir
    target = args.target
    name = args.name
    output_dir = args.output_dir
    qs_cutoff = args.qs_cutoff
    dockq_path = args.dockq_path
    lddt_path = args.lddt_path
    tmscore_path = args.tmscore_path
    n_cpu = args.n_cpu
    antibody = args.antibody
    chainAs = args.chainAs
    chainBs = args.chainBs
    rbm_version = args.rbm_version
    interface_weight = args.interface_weight

    save_ips_and_ics(input_dir, target, name, output_dir)
    save_qs_best(input_dir, target, name, output_dir, qs_cutoff, n_cpu)
    save_pairwise_interfaces(input_dir, target, name, output_dir)
    save_inputs(output_dir, target, dockq_path, lddt_path, tmscore_path, n_cpu)
    run_dockq(output_dir, target)
    run_lddt(output_dir, target)
    run_tmscore(output_dir, target)
    get_dockq_scores(target, output_dir)
    get_lddt_scores(target, output_dir)
    get_tm_scores(target, output_dir)
    save_interface_scores(input_dir, target, name, output_dir)

    if not antibody:
        if rbm_version == 'average':
            save_model_scores_v1(input_dir, target, name, output_dir, interface_weight)
        elif rbm_version == 'all':
            save_model_scores_v2(input_dir, target, name, output_dir, interface_weight)
        elif rbm_version == 'min':
            save_model_scores_v3(input_dir, target, name, output_dir, interface_weight)
        else:
            raise ValueError("Invalid RBM version")
        print("Completed RBM scoring")
        print("RBM version: " + rbm_version + " interface weight: " + interface_weight)
    else:
        if chainAs and chainBs:
            chainA_list = []
            chainB_list = []
            for chain in chainAs:
                chainA_list.append(chain)
            for chain in chainBs:
                chainB_list.append(chain)
            save_antibody_scores_v3(input_dir, target, name, output_dir, chainA_list, chainB_list)
        else:
            raise ValueError("Please provide valid chainAs and chainBs")