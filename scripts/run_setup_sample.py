#!/bin/python
#file to run pathways with MC

import sys
from optparse import OptionParser
from pathway_restraints import *
from sample_graph import Paths
from tables import open_file
from reaction_calculations import set_up_data
from reaction_calculations import ligandData
import pickle
import pathway_analysis as pa



transferase1="[C:1][O:2][H]>>[C:1][O:2][P]([O])([O])=O" #,
transferase2="[C:1][O:2][H]>>[C:1][O:2][S]([O])([O])=O"


def glyco_set_up(ligandfile, datafile):
    dockfiles = {}
    dockfiles['1'] = 'scores/bt0183_re2_scores.txt'
    dockfiles['2'] = 'scores/bt0184_re2_scores.txt'
    dockfiles['3'] = 'scores/bt0185_re2_scores.txt'
    lig = ligandData(datafile=datafile, smilesfile=ligandfile, delimiter=';', keepmultiples=True)
    lig.save_smiles()

#save docking scores to data file, index of column starts from 1, check delimiter,will carry forward multiples
#Normalize docking scores
#Saved docking scores corresponding to each enzyme

    lig.save_docking_scores(dockfiles, id_column=1, score_column=2, delimiter=';')
    #lig.save_reaction_scores(reactions)

    #set_up_data(ligandfile, ',', dockfiles,None,None, datafile, useMorganFP=True)
    #results = pd.read_hdf(datafile, 'rsim')
    print datafile

#run glycolysis run 
def run_standard_glycolysis(number, steps, prots, ligandfile, datafile, fileprefix,
                            seed=None, parameter_set=0, kparameter_set=0):
    outfile = '%s%d.h5' % (fileprefix, int(number))
    
    prots = ['1' , '2', '3'] 

    protnum = len(prots)

    d = dataTables(proteins=prots)
    d.readInLigands(ligandfile, delimiter=';')
    d.readInTables(datafile)

    print outfile
    with open_file(outfile, 'w') as h5outhandle:
        ps = Paths(d, h5outhandle, seed)
        dr = dock_restraint()
        restraints = [dr]
        ratio = ps.sample_labels_monte_carlo(steps, protnum, protnum+1,
                                             restraints,
                                             parameter_set=parameter_set,
                                             kparameter_set=kparameter_set)

#setup random pathway sampling for convergence testing
def run_random_glycolysis(steps, prots, ligandfile, datafile, randomh5, randompicklefile):
    
    prots = ['1', '2', '3']  
    
    print 'Sampling %d random graphs' % steps

    d = dataTables(proteins=prots)
    d.readInLigands(ligandfile, delimiter=';')
    d.readInSelectedTables(datafile, ['dock'])

    with open_file(randomh5, mode='w') as h5outhandle:
        ps = Paths(d, h5outhandle)
        dr = dock_restraint()
        restraints = [dr]
        ps.random_graphs(steps, len(prots), len(prots)+1, restraints)
    
    scorelist = []
    print 'Getting stats'
    with open_file(randomh5, 'r') as h5file:
        table = h5file.root.paths.pathTable
        for row in table.iterrows():
            scorelist.append(row['obj'])
    sl = np.array(scorelist)
    pickle.dump(sl, open(randompicklefile, 'w'))
    print np.mean(sl), np.std(sl)

#if you know the right pathway

#analyze pathways to generate clusrer output files
def analyze(fileprefix, picklefile, randompicklefile, clusteroutfile):
    from analyze_test import cluster_pathways
    from analyze_test import pathens
    from analyze_test import convergence_test
    from analyze_test import print_pathways_in_cluster
    numstd = 1.5
    totalruns = 10

    best = None
    for i in range(1, totalruns+1):
        currfile = '%s%d.h5' % (fileprefix, i)
        with open_file(currfile, 'r') as h5:
            currtable = h5.root.paths.pathTable
             	     
            runbest = currtable.attrs.bestscore
            best = np.max([runbest, best])
	    print (str(best) + "printing best")
            print (str(runbest) + "printing runbest")
    sl = pickle.load(open(randompicklefile, 'r'))
    stdev = np.std(sl)
    cutoff = best - numstd*stdev
    print (cutoff) 
		     
 
    numiters = min(2, totalruns)
    convergence_test(cutoff=cutoff,
                     outpickle=picklefile,
                     fileprefix=fileprefix,
                     totalruns=totalruns,
                     maxgroupsize=totalruns,
                     numiters=numiters,
                     outfile=clusteroutfile)
    print_pathways_in_cluster(picklefile,1,0.5)
    print ("second cluster")
    print_pathways_in_cluster(picklefile,2,0.5)
    print ("third cluster")
    print_pathways_in_cluster(picklefile,3,0.5)

def main():
    n = sys.argv[1]
    s = sys.argv[2]
    glyco_set_up(ligandfile, datafile)
 
if __name__ == "__main__":

    ##############################################################
    # FILENAMES
    ligandfile = 'scores/smi_all_ID1.txt'

    datafile = 'test_docking.h5'
    fileprefix = 'test_output/test_'
    randomh5 = 'test_output/test_output.h5'
    randompicklefile = 'test_output/test_output.pickle'
    picklefile = 'test_output/test_dock.pickle'
    clusteroutfile = 'test_output/test_clusters_dock.txt'

    #dummydatafile = 'glyco/gly_data_dummy.h5'
    #dummyfileprefix = '/scrapp2/ganesans/pathways/glyco_dummy_'
    #dummyrandomh5 = '/scrapp2/ganesans/pathways/glyco_random_dummy.h5'
    #dummyrandompicklefile = '/scrapp2/ganesans/pathways/glyco_random_dummy.pickle'
    #dummypicklefile = 'output/gly_dummy.pickle'
    #dummyclusteroutfile = 'output/glycolysis_clusters_dummy.txt'

    ##############################################################

    proteins = ['1' , '2', '3'] 

    parser = OptionParser()
    parser.add_option('-e', action="store", dest="function_call",
                      help='Options: setup, run, random, analyze')
    parser.add_option('-d', '--dummy', action="store_true", default=False, dest="dummy")
    parser.add_option('-n', '--run_number', action='store', dest='run_number')
    parser.add_option('-s', '--steps', action="store", type='int', dest='num_steps')

    options, args = parser.parse_args()

    if options.dummy:
        datafile = dummydatafile
        fileprefix = dummyfileprefix
        randomh5 = dummyrandomh5
        randompicklefile = dummyrandompicklefile
        picklefile = dummypicklefile
        clusteroutfile = dummyclusteroutfile
        
        proteins = ['1', 'D', '3']
 
    if options.function_call == 'setup':
        if options.dummy:
            set_up_gly_dummy(ligandfile, seafile, datafile)
        else:
            glyco_set_up(ligandfile,datafile)
    elif options.function_call == 'run':
        n = options.run_number
        s = options.num_steps
        if (n is None) or (s is None):
            if len(args) > 1:
                n = args[0]
                s = int(args[1])
        run_standard_glycolysis(n, int(s), proteins, ligandfile, datafile, fileprefix,
                                kparameter_set=4, parameter_set=1)
    elif options.function_call == 'random':
        s = options.num_steps
        if s is None:
            if len(args) > 0:
                s = int(args[0])
        run_random_glycolysis(s, proteins, ligandfile, datafile, randomh5, randompicklefile)
    elif options.function_call == 'analyze':
        analyze(fileprefix, picklefile, randompicklefile, clusteroutfile)
    else:
        true_pathway(datafile)



