import os
import sys

def get_lddt_scores(target, output_dir):
    fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + target + '.list', 'r')
    cases = []
    for line in fp:
        words = line.split()
        model = words[0]
        pair1 = words[1]
        pair2 = words[2]
        cases.append([model, pair1, pair2])
    fp.close()

    os.makedirs(output_dir + '/' + target + '/lDDT', exist_ok=True)
    rp = open(output_dir + '/' + target + '/lDDT/' + target + '.result', 'w')
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

        lddt = ''
        if os.path.exists(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + model + '/' + pair1 + '_' + pair2 + '.lddt'):
            fp = open(output_dir + '/' + target + '/' + 'pairwise_interfaces_for_lddt' + '/' + model + '/' + pair1 + '_' + pair2 + '.lddt', 'r')
            for line in fp:
                words = line.split()
                if len(words) >= 4:
                    if words[0] == 'Global' and words[1] == 'LDDT' and words[2] == 'score:':
                        lddt = words[3]
                        break

            fp.close()
        if lddt:
            rp.write(model + '\t' + pair1.replace('-',':') + '\t' + pair2.replace('-',':') + '\t' + lddt + '\n')
            rp.write(model + '\t' + newpair1 + '\t' + newpair2 + '\t' + lddt + '\n')
        else:
            print (target, model, pair1, pair2)
    rp.close()

if __name__ == "__main__":
    target = sys.argv[1]
    output_dir = sys.argv[2]
    get_lddt_scores(target, output_dir)