import os
import sys





def save_inputs(output_dir, target, dockq_path, lddt_path, tmscore_path, n_cpu):
    print (f'Saving inputs for {target}...')
    if not os.path.exists(output_dir + '/' + target + '/' + 'cmds'):
        os.makedirs(output_dir + '/' + target + '/' + 'cmds')
    rp1 = open(output_dir + '/' + target + '/' + 'cmds' + '/DockQ_inputs','w')
    rp2 = open(output_dir + '/' + target + '/' + 'cmds' + '/lDDT_inputs', 'w')
    rp3 = open(output_dir + '/' + target + '/' + 'cmds' + '/TMscore_inputs', 'w')
    fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + target + '.list','r')
    for line in fp:
        words = line.split()
        # print (f'{output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.dockq')
        if os.path.exists(f'{output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.dockq'):
            pass
        else:
            rp1.write(f'{dockq_path} --no_align --n_cpu 1 {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}_model.pdb {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}_target.pdb > {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.dockq 2> {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.log\n')
        # print (f'{output_dir}/{target}/pairwise_interfaces_for_lddt/{words[0]}/{words[1]}_{words[2]}.lddt')
        if os.path.exists(f'{output_dir}/{target}/pairwise_interfaces_for_lddt/{words[0]}/{words[1]}_{words[2]}.lddt'):
            pass
        else:
            rp2.write(f'{lddt_path} -x {output_dir}/{target}/pairwise_interfaces_for_lddt/{words[0]}/{words[1]}_{words[2]}_model.pdb {output_dir}/{target}/pairwise_interfaces_for_lddt/{words[0]}/{words[1]}_{words[2]}_target.pdb > {output_dir}/{target}/pairwise_interfaces_for_lddt/{words[0]}/{words[1]}_{words[2]}.lddt\n')
        # print (f'{output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.tm')
        if os.path.exists(f'{output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.tm'):
            pass
        else:
            rp3.write(f'{tmscore_path} -c -ter 0 {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}_model.pdb {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}_target.pdb > {output_dir}/{target}/pairwise_interfaces/{words[0]}/{words[1]}_{words[2]}.tm\n')
    fp.close()
    rp1.close()
    rp2.close()
    rp3.close()
    if not os.path.exists(output_dir + '/' + target + '/' + 'cmds'):
        os.makedirs(output_dir + '/' + target + '/' + 'cmds')
    with open(output_dir + '/' + target + '/' + 'cmds' + '/DockQ_cmds.sh','w') as rp:
        rp.write(f'python batch_run.py {output_dir}/{target}/cmds/DockQ_inputs {n_cpu}\n')
    with open(output_dir + '/' + target + '/' + 'cmds' + '/lDDT_cmds.sh','w') as rp:
        rp.write(f'python batch_run.py {output_dir}/{target}/cmds/lDDT_inputs {n_cpu}\n')
    with open(output_dir + '/' + target + '/' + 'cmds' + '/TMscore_cmds.sh','w') as rp:
        rp.write(f'python batch_run.py {output_dir}/{target}/cmds/TMscore_inputs {n_cpu}\n')


def run_dockq(output_dir, target):
    print (f'Running DockQ for {target}...')
    os.system(f"bash {output_dir}/{target}/cmds/DockQ_cmds.sh")

def run_lddt(output_dir, target):
    print (f'Running lDDT for {target}...')
    os.system(f"bash {output_dir}/{target}/cmds/lDDT_cmds.sh")

def run_tmscore(output_dir, target):
    print (f'Running TMscore for {target}...')
    os.system(f"bash {output_dir}/{target}/cmds/TMscore_cmds.sh")

if __name__ == "__main__":
    output_dir = sys.argv[1]
    target = sys.argv[2]
    dockq_path = sys.argv[3]
    lddt_path = sys.argv[4]
    tmscore_path = sys.argv[5]
    n_cpu = sys.argv[6]
    save_inputs(output_dir, target, dockq_path, lddt_path, tmscore_path, n_cpu)