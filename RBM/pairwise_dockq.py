import os
import sys

def get_dockq_scores(target, output_dir):
    fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + target + '.list', 'r')
    cases = []
    for line in fp:
        words = line.split()
        model = words[0]
        pair1 = words[1]
        pair2 = words[2]
        cases.append([model, pair1, pair2])
    fp.close()
    os.makedirs(output_dir + '/' + target + '/DockQ', exist_ok=True)
    rp = open(output_dir + '/' + target + '/DockQ/' + target + '.result', 'w')
    for case in cases:
        model = case[0]
        pair1 = case[1]
        prot1A = pair1.split('-')[0]
        prot1B = pair1.split('-')[1]
        newpair1 = prot1B + ':' + prot1A
        pair2 = case[2]
        prot2A = pair2.split('-')[0]
        prot2B = pair2.split('-')[1]
        newpair2 = prot2B + ':' + prot2A

        dockq = ''
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.dockq'):
            fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces' + '/' + model + '/' + pair1 + '_' + pair2 + '.dockq', 'r')
            for line in fp:
                if line[0] == '*':
                    pass
                else:
                    words = line.split()
                    if len(words) > 7:
                        if words[0] == 'Total' and words[1] == 'DockQ' and words[2] == 'over':
                            dockq = words[6]
                            break
            fp.close()
        if dockq:
            rp.write(model + '\t' + pair1.replace('-',':') + '\t' + pair2.replace('-',':') + '\t' + dockq + '\n')
            rp.write(model + '\t' + newpair1 + '\t' + newpair2 + '\t' + dockq + '\n')
        else:
            rp.write(model + '\t' + pair1.replace('-',':') + '\t' + pair2.replace('-',':') + '\t0.0\n')
            rp.write(model + '\t' + newpair1 + '\t' + newpair2 + '\t0.0\n')
            print (target, model, pair1, pair2)
    rp.close()

# if __name__ == "__main__":
#     target = sys.argv[1]
#     output_dir = sys.argv[2]
#     get_dockq_scores(target, output_dir)