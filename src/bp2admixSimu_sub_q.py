if __name__ == '__main__':

    import argparse
    import pandas as pd 
    import numpy as np 
    import admixString_generators as strGen
    import parsers as prs
    import re
    import random

    parser = argparse.ArgumentParser(description="takes in a bp file and creates admixsimu format bp files with ancestry states for all possible ancestry strings and t1,t2values consistent with the input bpgen file. Can be run with or without subadmixture.Also creates a .txt table file with bp-headers, admixture strings, p1, p2, t1,t2 columns.")
    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file - must be a bp file")
    parser.add_argument("--t1", "-t1", required=True, type=int, help="number of observed meioses on the larger t branch")
    parser.add_argument("--t2", "-t2", required=True, type=int, help="number of observed meioses on the smaller t branch")
    #parser.add_argument('-o','--options', nargs='+', help='options for subadmixture - e.g -o n s - for normal and subadmixed on the respective branches', required=True)

    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 
    
#{{{ 

    s_list=strGen.generateAdmixStrings_branch_subadmix_50(args.t1-1,args.t2-1) #t-1 for subbranch because strings are duplicated below
    admixstr_list=[x+y for x,y in zip(s_list,s_list)]    
    print('Generated strings!')
    with open(args.inputfile, 'r') as f:
        lines_all = f.readlines()
        
    nfams= (len(lines_all)//6) #Batch size of individuals for each admixtype. 2 haplotypes for the individual, and 4 for his/her parents

    tablelines = [['header\tadmixstring\tpA\tpB\tt']]
    writelines = []

    for fam in range(nfams):
        admixstr = admixstr_list[fam%len(admixstr_list)] 
        start = fam*6
        end = start + 6
        famlines = lines_all[start:end]

        famlines_anc = [] 
        for j in range(len(famlines)):

            if j%2 != 0: 
                continue
            famlines_anc.append(prs.ancestryParse([famlines[j], famlines[j+1]],list(admixstr)))

        famlines_anc_flat = [x for y in famlines_anc for x in y]
        famlines_haplodicts = [prs.haplo_parse(x) for x in famlines_anc_flat]

        print('famlines haplo parsed')

        ps_0,ps_1 = prs.computeAncprops(famlines_haplodicts[:4])
        assert 0.99 <= ps_0[0] + ps_1[0] <= 1.01, "t={}, admixstr = {}, header={}".format(str(args.t1), admixstr, famlines_anc_flat[0].split(' ')[0])

        pA = sum(ps_0[:2])/2
        pB = sum(ps_0[2:])/2

        tablelines.append([x.split(' ')[0] +'\t'+ admixstr + '\t' +str(pA) + '\t' + str(pB) + '\t' + str(args.t1)  for x in famlines_anc_flat])
        writelines.append(famlines_anc_flat[4:])


    tablefile = str(args.t1) + '_' + 'intf_sxspf.inds'

    with open(tablefile, 'w+') as f:
        f.write("\n".join([x for y in tablelines for x in y]))

    writelines_flat = [x for y in writelines for x in y]
    outfile = str(args.t1) + '_' + 'intf_sxspf.asu.bp'
    with open(outfile, 'w+') as f:
        f.write("".join(writelines_flat))

        #Write file for current typ
#        infile = re.search('\/([^\/]*$)', args.inputfile).group(1)
#}}}
