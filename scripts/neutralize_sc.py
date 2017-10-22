usage = """ contribution from Hans de Winter

    FROM RDKit Cookbook (http://www.rdkit.org/docs/Cookbook.html)

    Retrieved 11/18/2013


##########################################################

usage: python neutralize.py <compounds>.smi
	<compounds>.smi has the following format with no header:
	[smiles];[code]

output:
	neutralized.smi
"""
from optparse import OptionParser
from rdkit import Chem
from rdkit.Chem import AllChem

def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutralizeCharges(mol, code, badf, line_number, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions

    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        try:
            while mol.HasSubstructMatch(reactant):
#                print "In smiles: '%s' replace '%s' -> '%s'" % (smiles, reactant, product)
                replaced = True
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]
        except:
            print "failed to neutralize molecule with smiles '%s' for code '%s' on line '%s'" % (mol, code, line_number)
            break

    return Chem.MolToSmiles(mol,True)

def run_smi(smi_fname, bad_fname, out_smi_fname, delimiter):
    smi_f = open(smi_fname)
    out_smi_f = open(out_smi_fname, 'w')
    badf = open(bad_fname, 'w')
    line_number = 1
    for line in smi_f:
        try:
            smiles, code = line[:-1].split(delimiter)
        except:
            print "Unable to upack line %s: '%s' into <smiles>;<code>" % ( str(line_number), line[:-1])
            badf.write(line)
            continue
            

        try:
            mol = Chem.MolFromSmiles(smiles)
        
        except:
            badf.write("%s;%s\n" % (smiles, code))
            out_smi_f.write("%s;%s\n" % (smiles, code))
            continue

        if mol is None:
            print "failed to read smiles '%s' for code '%s' on line '%s'" % (
                smiles, code, line_number)
            badf.write("%s;%s\n" % (smiles, code))
            continue

        neutralized_smiles  = NeutralizeCharges(
            mol, code, badf, line_number) 
        out_smi_f.write("%s;%s\n" % (neutralized_smiles, code))
        line_number += 1

    smi_f.close()
    out_smi_f.close()

def run_sdf(sdf_fname, bad_fname, out_smi_fname):
    out_smi_f = open(out_smi_fname, 'w')
    badf = open(bad_fname, 'w')

    suppl = Chem.SDMolSupplier(sdf_fname)
    molecule_number = 1
    for mol in suppl:
        if mol is None:
            print "FAILED to parse molecule # %s " % molecule_number
            continue

        code = mol.GetProp("_Name")
        neutralized_smiles  = NeutralizeCharges(mol, code, badf, molecule_number) 
        out_smi_f.write("%s;%s\n" % (neutralized_smiles, code))
        molecule_number += 1

    out_smi_f.close()


if __name__ == '__main__':
    parser = OptionParser(usage)
    
    parser.add_option('-o','--outfile',dest='outfile',metavar='FILE',
                      help='Specify outfile file name (default %default)',
                      default="cleaned.smi")

    parser.add_option('-b','--badfile',dest='badfile',metavar='FILE',
                      help='Specify badfile file name (default <outfile>.bad)')

    parser.add_option('-d','--delimiter',dest='delimiter',metavar='FILE',
                      help='Specify delimiter in line (default %default)',
                      default=";")

    options,args = parser.parse_args()

    try:
        input_fname, = args
    except:
        parser.error("Please provide a single argument of the .smi file.\nYou provided: (" + ", ".join(args) + ")")
        print usage

    if options.badfile is None or options.badfile == "":
        badfile = options.outfile + ".bad"
    else:
        badfile = options.badfile

    if  input_fname[-4:] == ".smi":
        print 'Running'
        print options.outfile
        run_smi(input_fname, badfile, options.outfile, delimiter=options.delimiter)
    else:
        print "Unrecognized extension: '%s'" % input_fname[-4:]
        
