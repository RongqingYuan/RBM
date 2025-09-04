import os
import sys

target_list = sys.argv[1]
output_dir = sys.argv[2]

fp = open(target_list,'r')
targets = []
for line in fp:
    words = line.split()
    targets.append(words[0])
fp.close()
targets = ['H0208']

with open('DockQ_cmds','w') as rp:
    rp.write('python batch_run.py DockQ_inputs 64\n')
with open('lDDT_cmds','w') as rp:
    rp.write('python batch_run.py lDDT_inputs 64\n')
with open('TMscore_cmds','w') as rp:
    rp.write('python batch_run.py TMscore_inputs 64\n')

rp1 = open('DockQ_inputs','w')
rp2 = open('lDDT_inputs', 'w')
rp3 = open('TMscore_inputs', 'w')
for target in targets:
    fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + target + '.list','r')
    for line in fp:
        words = line.split()

        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.dockq'):
            pass
        else:
            rp1.write('DockQ --no_align --n_cpu 4 ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_model.pdb ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_target.pdb > ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.dockq 2> ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.log\n')
        
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces_1' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.lddt'):
            pass
        else:
            rp2.write('lddt -x ' + output_dir + '/' + target + '/' + 'pairwise_interfaces_1' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_model.pdb ' + output_dir + '/' + target + '/' + 'pairwise_interfaces_1' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_target.pdb > ' + output_dir + '/' + target + '/' + 'pairwise_interfaces_1' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.lddt\n')
        
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.tm'):
            pass
        else:
            rp3.write('/home2/s439906/software/USalign/TMscore -c -ter 0 ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_model.pdb ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '_target.pdb > ' + output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + words[0] + '/' + words[1] + '_' + words[2] + '.tm\n')
    fp.close()
rp1.close()
rp2.close()
rp3.close()


