# Reformat scores so in the form of:
# ID;score


def reformat_file(infile, outfile, idcol, scorecol):
    scoredict = {}
    with open(infile, 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            if not line.startswith('#'):
                line = line.strip()
                stringlist = line.split()
                identifier = stringlist[idcol-1]
                score = stringlist[scorecol-1]
                scoredict[identifier] = score

    with open(outfile, 'w') as outhandle:
        for key in scoredict.keys():
            identifier = key
            score = scoredict[key]
            outhandle.write('%s,%s\n' % (identifier, score))

if __name__ == "__main__":
    infile = '6_gi1612_9733_rescore_duplicateyes.txt'
    outfile = 'gly6_scores.txt'
    filedict = {'1_gi1613_0320_h2o45_rescore_duplicatesyes.txt':'gly1_scores.txt',
                '2_gi1613_1851_rescore_duplicatesyes.txt':'gly2_scores.txt',
                '3_gi16131754_rescore_duplicateyes.txt':'gly3_scores.txt',
                '4_gi1613_0826_rescore_duplicateyes.txt':'gly4_scores.txt',
                '5_gi1613_1757_rescore_duplicateyes.rept':'gly5_scores.txt',
                '6_gi1612_9733_rescore_duplicateyes.txt':'gly6_scores.txt',
                '7_gi1613_0827_rescore_duplicateyes.txt':'gly7_scores.txt',
                '8_gi1613_1483_rescore_duplicateyes.txt':'gly8_scorse.txt',
                '9_gi1613_0686_rescore_duplicateyes.txt':'gly9_scores.txt',
                '10_gi1612_9632_rescore_duplicateyes.txt':'gly10_scores.txt'}
    for f in filedict.keys():
        print f
        reformat_file(f, filedict[f], idcol=7, scorecol=5)
