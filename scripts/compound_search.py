#!/bin/python

from openeye.oechem import *
from openeye.oegraphsim import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from read_ligand_data import ligandData
import time

def make_molecule(smi):
    mol = OEGraphMol()
    OEParseSmiles(mol, smi)
    OEAssignImplicitHydrogens(mol)
    OEAssignFormalCharges(mol)
    return mol

def check_molecule(m, mol):
    match = False
    for other in mol:
        match = match or OEExactGraphMatch(m, other)

    return match

def smilesToRDKitMol(smiles):
    rdkmol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(rdkmol)
    AllChem.AssignStereochemistry(rdkmol)
    return rdkmol

def OEMolToRDKitMol(mol):
    smiles = OECreateIsoSmiString(mol)
    return smilesToRDKitMol(smiles)

#def smilesToMorganFingerprint(smiles):
#    smiles = pybel.readstring('smi', smiles).write('can')
#    rdkmol = smilesToRDKitMol(smiles)
#    fp = AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True)
#    return fp

def standardizeMol(mol):
    OEAssignAromaticFlags(mol)
    OEAssignFormalCharges(mol)
    return mol

def OEMolToMorganFingerprint(mol):
    smiles = OECreateIsoSmiString(mol)
    rdkmol = Chem.MolFromSmiles(smiles)
    if rdkmol:
        Chem.SanitizeMol(rdkmol)
        AllChem.AssignStereochemistry(rdkmol)
        fp = AllChem.GetMorganFingerprintAsBitVect(rdkmol, 2, useChirality=True)
        return fp
    else:
        return False

def write_out_candidates(filename, candidates, keggdata):
    print 'Writing out %d compound candidates to %s' % (len(candidates), filename)
    with open(filename, 'w') as handle:
        for c in candidates:
            handle.write('%s\t%s\n' % (OECreateIsoSmiString(keggdata.molecules[c]), c))

class enzymeData(object):
    def __init__(self, dockingfile=None, reactionlist=None):
        self.dockingfile = dockingfile
        self.reactionlist = reactionlist
        self.molecules = []

    def get_molids(self, column=0, cutoff=100, delimiter=None):
        molids = []
        
        with open(self.dockingfile, 'r') as inhandle:
            lines = inhandle.readlines()
            cutoff = min(len(lines), cutoff)
            for i in range(cutoff):
                if not delimiter:
                    splitline = lines[i].split()
                else:
                    splitline = lines[i].split(delimiter)
                ligid = splitline[column].strip()
                molids.append(ligid)
        self.molecules = molids

    def find_similar_compounds(self, keggdata, tcmin=0.95):
        candidates = []
        for moleculeid in self.molecules:
            if moleculeid not in keggdata.molecules.keys():
                print 'Warning: %s not found' % moleculeid
            else:
                substrate = keggdata.molecules[moleculeid]
                fps = []
                for r in self.reactionlist:
                    lg = OELibraryGen(r)
                    lg.AddStartingMaterial(substrate, 0)
                    lg.SetValenceCorrection(True)
                    for p in lg.GetProducts():
                        numparts, partlist = OEDetermineComponents(p)
                        if numparts == 1:
                            #OEAssignFormalCharges(p)
                            fp = OEFingerPrint()
                            OEMakeFP(fp, p, OEFPType_Circular)
                            fps.append(fp)
                        elif numparts > 1:
                            pred = OEPartPredAtom(partlist)
                            for i in range(1, numparts + 1):
                                pred.SelectPart(i)
                                partmol = OEGraphMol()
                                OESubsetMol(partmol, p, pred)
                                #OEAssignFormalCharges(partmol)
                                fp = OEFingerPrint()
                                OEMakeFP(fp, partmol, OEFPType_Circular)
                                fps.append(fp)
                for fp in fps:
                    for kid in keggdata.fingerprints.keys():
                        tc = OETanimoto(fp, keggdata.fingerprints[kid])
                        if tc > tcmin:
                            #print kid, tc
                            #print keggdata.smiles[kid]
                            candidates.append(kid)
        candidates = list(set(candidates))
        return candidates

    def possible_reactants(self, keggdata):
        candidates = []
        for moleculeid in keggdata.molecules.keys():
            substrate = keggdata.molecules[moleculeid]
            for r in self.reactionlist:
                lg = OELibraryGen(r)
                m = lg.AddStartingMaterial(substrate, 0)
                if m:
                    candidates.append(moleculeid)
                    break
        return candidates

  
def verei():
    print 'Reading ligand data'
    keggsmifile = 'background_data/kegg_sep2007.smi'
    keggdata = ligandData(smilesfile=keggsmifile, delimiter=';')
    print 'Finished reading in ligand data'
    
    starttime = time.time()
    candidates = []
    rxn1 = ['[C:1][N:2][C](=[O])[N]>>[C:1][N:2][H]']
    rxn2 = ['[A:1][C:2]([H])([O:3][H])[C:4]([OH])[A:5]>>[A:1][C:2](=[O:3])[C:4]([H])[A:5]']
    rxn3 = ['[A:1][C:2]=[O:3]>>[A:1][C@:2]([H])[N]',
            '[A:1][C:2](N)([A:3])C(=O)O>>[A:1][C:2](=O)[A:3]',
            '[N:1][C:2]([H])[C:3](=[O:4])[A:5]>>[O:4]=[C:2][C@:3]([N:1])([H])[A:5]',
            '[N:1][C:2]([H])[C:3](=[O:4])[A:5]>>[O:4]=[C:2][C@@:3]([N:1])([H])[A:5]']
    rxn4 = ['[O:1]=[C:2]([H])[A:3]>>[O:1]=[C:2]([OH])[A:3]']
    dockfiles = ['verei/verei_gi121608226_rescore.txt',
                 'verei/verei_gi121608227_rescore.txt',
                 'verei/verei_gi121608228_rescore.txt',
                 'verei/verei_gi121608229_rescore.txt']
    #enzs = [enzymeData(dockingfile=dockfiles[0], reactionlist=rxn1),
    #        enzymeData(dockingfile=dockfiles[1], reactionlist=rxn2),
    #        enzymeData(dockingfile=dockfiles[2], reactionlist=rxn3),
    #        enzymeData(dockingfile=dockfiles[3], reactionlist=rxn4)]
    allrxns = rxn1 + rxn2 + rxn3 +rxn4
    enzs = [enzymeData(dockingfile='verei/all.smi', reactionlist=allrxns)]
    for e in enzs:
        e.get_molids(column=1, cutoff=10000)
        newc = e.find_similar_compounds(keggdata)
        candidates.extend(newc)
        print '> Number of candidates: %d' % len(newc)
    print time.time() - starttime
    candidates = list(set(candidates))
    write_out_candidates('verei_ligands.txt', candidates, keggdata)
    

def gmh():
    print 'Reading in kegg data with zinc ids'
    keggsmifile = 'background_data/ZINC_all.smi'
    keggdata = ligandData(smilesfile=keggsmifile)
    print 'Finished reading in ligand data'

    starttime = time.time()
    candidates = []
    isomerase=['[C:1]([H])[C:2](=[O:3])[C:4][C:5][C:6][O:7][H]>>[O:3]([H])[C:2]1([H])[C:1][O:7][C:6][C:5][C:4]1']
    transferase=["[C:1][O:2][H]>>[C:1][O:2][P]([O-])([O-])=O"]
    esterase=["[C:1][O:2][P]([O-])([O-])=O>>[C:1][O:2][H]"]
    guanylyltransferase=["[O:1][P:2]([O:3])([O:4])(=[O:5])>>C1=NC2=C(N1C3C(C(C(O3)COP(O)(=O)[O:3][P:2]([O:4])([O:1])=[O:5])O)O)NC(=NC2=O)N"]

    enzs = [enzymeData(dockingfile='gmh_data/combine.scores.1', reactionlist=isomerase),
            enzymeData(dockingfile='gmh_data/combine.scores.3', reactionlist=transferase),
            enzymeData(dockingfile='gmh_data/combine.scores.4', reactionlist=esterase),
            enzymeData(dockingfile='gmh_data/combine.scores.5', reactionlist=guanylyltransferase)]
    for e in enzs:
        e.get_molids(column=0, cutoff=1000)
        e.molecules = ['ZIN%s' % m for m in e.molecules]
        newc = e.find_similar_compounds(keggdata, tcmin=0.85)
        candidates.extend(newc)
        print '> Number of candidates: %d' % len(newc)
    print time.time() - starttime
    candidates = list(set(candidates))
    write_out_candidates('gmh_ligands.txt', candidates, keggdata)

def kdo():
    print 'Reading in kegg data with zinc ids'
    keggsmifile = 'background_data/ZINC_all_neutral.smi'
    keggdata = ligandData(smilesfile=keggsmifile)
    print 'Finished reading in ligand data'
    starttime = time.time()
    candidates = []
    
    kdsf=['[H][O:1][C:2]([H])[C:3]=[O:4]>>[O:1]=[C:2][C:3]([H])[O:4][H]']
    kdsa=['[C:1]=[O:2]>>[C:1]([O:2][H])CC(=O)C(=O)O']
    phosphatase=['[C:1][O:2][P]([O])([O])=O>>[C:1][O:2]']
    kdsb=['[H][O:1][C:2][A:5][A:6][A:7][C:3]=[O:4]>>[O:1]1[C:2][A:5][A:6][A:7][C:3]1[O:4]P(=O)(OC[C@H]1O[C@H]([C@@H]([C@@H]1O)O)n1c(=O)nc(cc1)N)O']
    
    #enzs = [enzymeData(dockingfile='kdo8p_data/combine.scores.1', reactionlist=kdsf),
    #        enzymeData(dockingfile='kdo8p_data/combine.scores.3', reactionlist=phosphatase),
    #        enzymeData(dockingfile='kdo8p_data/combine.scores.4', reactionlist=kdsa),
    #        enzymeData(dockingfile='kdo8p_data/combine.scores.6a', reactionlist=kdsb)]

    smilesfile='kdo8p_data/top1500_smiles.smi'
    enzs = [enzymeData(dockingfile=smilesfile, reactionlist=kdsf),
            enzymeData(dockingfile=smilesfile, reactionlist=phosphatase),
            enzymeData(dockingfile=smilesfile, reactionlist=kdsa),
            enzymeData(dockingfile=smilesfile, reactionlist=kdsb)]

    for e in enzs:
        e.get_molids(column=1, cutoff=5000)
        e.molecules = ['%s' % m for m in e.molecules]
        newc = e.find_similar_compounds(keggdata, tcmin=0.75)
        candidates.extend(newc)
        print '> Number of candidates: %d' % len(newc)
    print time.time() - starttime
    candidates = list(set(candidates))
    write_out_candidates('kdo_ligands.txt', candidates, keggdata)


def this_is_a_test():
    oxidoreductase=['[H][C:1][O:2][H]>>[C:1]=[O:2]']
    aminotransferase=['[C:1]=[O]>>[H][C:1][NH2]']
    phosphatase=['[C:1][O:2][P]([O])([O])=O>>[C:1][O:2][H]']
    acetyltransferase=['[C:1][O:2][H]>>[C:1][O:2][C]([CH3])=O']
    acetatelyase=['[C:1][O][C]([CH3])=O>>[C:1][S][H]']

    keggsmifile = 'background_data/ZINC_all_neutral.smi'
    keggdata = ligandData(smilesfile=keggsmifile)
    reactions = oxidoreductase + aminotransferase + phosphatase + acetyltransferase + acetatelyase
    e = enzymeData(dockingfile=None, reactionlist=reactions)
    candidates = e.possible_reactants(keggdata)
    write_out_candidates('./test_reactability.smi', candidates, keggdata)

if __name__ == "__main__":
    #gmh()
    #verei()
    #kdo()
    #gulonate()
    starttime = time.time()
    #this_is_a_test()
    print time.time() - starttime
