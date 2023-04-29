import numpy as np

from SQRNdbnseq import SQRNdbnseq


if __name__ == "__main__":

    from collections import Counter

    rst = None

    threads = 1

    #SAM riboswitch
    seq = "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAACA"
    dbn = "(((((((((....(((((...(((.[[[[)))......)))))(((..(((((...(((((....))))).)))..)).)))...(]]]](((((.......)))))..))))))))))." 
    #TRNA
    #seq = "GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA"
    #dbn = "(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))...."
    #Twister ribozyme
    #seq = "GAAAUAAUGCUACCAGACUGCACAAAAAGUCUGCCCGUUCCAAGUCGGGAUGAAAGCAGAGGAUUUC"
    #dbn = "((((...(((({{((((((........))))))(((..[[[..}}.))).....))))..]]]))))"
    #rst = "(x((xxx....xx..............,,,,,,..............................)).)"
    #rst = "..................................................................." 
    
    queue = []

    with open("CoRToise.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-2,3):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            db = lines[ii+2].strip()
            if len(sq) < 150:
                queue.append([nm, sq, db, rst])

    #seq = "AAACCACGAGGAAGAGAGGUAGCGUUUUCUCCUGAGCGUGAAGCCGGCUUUCUGGCGUUGCUUGGCUGCAACUGCCGUCAGCCAUUGAUGAUCGUUCUUCUCUCCGUAUUGGGGAGUGAGAGGGAGAGAACGCGGUCUGAGUGGU"
    #dbn = "..(((((......(.......((((((((((((....(...(...((((..(.((((((((......))))).))))..))))..)......)..((((((((((.....)))))).))))))))))))))))...)...)))))"

    seq = "GGUAAAGAAUGAAAAAACACGAUUCGGUUGGUAGUCCGGAUGCAUGAUUGAGAAUGUCAGUAACCUUCCCCUCCUCGGGAUGUCCAUCAUUCUUUAAUAUCUUUUAUGAGGAGGGAA"
    dbn = None
    
    #queue  = [["Bs_yjdf_GGARBS", seq, dbn, rst],]
    #queue += [["default", seq, dbn, rst],]

    """ TOP FROM FIRST MANUAL TUNING
    bpweights = {
                 'GU' : -1.95,
                 'AU' :  1.5,
                 'GC' :  4.25,
                 }


    subopt = 0.9
    toplim  = 5
    conslim = 1
    
    minlen = 2
    minbpscore = 8
    minfinscorefactor = 0.75
    distcoef = 0.18
    bracketweight = 19.0
    orderpenalty = -0.08
    fiveprime = 0.215
    gupen  = 0.0
    maxstemnum = 10**6"""

    """ Grid best on 230429 """
    bpweights = {
                 'GU' : -1,
                 'AU' :  2,
                 'GC' :  4,
                 }


    subopt = 0.85
    toplim  = 5
    conslim = 1
    
    minlen = 2
    minbpscore = 6
    minfinscorefactor = 1
    distcoef = 0.1
    bracketweight = -1
    orderpenalty = 1
    fiveprime = 0.0
    gupen  = 0.0
    maxstemnum = 10**6
    

    resultsB = []
    resultsC = []

    for obj in queue:

        name, seq, dbn, rst = obj

        result = SQRNdbnseq(seq, bpweights, rst, dbn,
                            subopt, minlen, minbpscore,
                            minfinscorefactor, bracketweight,
                            distcoef, orderpenalty, fiveprime,
                            maxstemnum, gupen, conslim, toplim,
                            threads)

        print(name)
        print(seq)
        if dbn:
            print(dbn)
        print('='*len(seq))
        print(result[0],result[2])
        print('_'*len(seq))
        for rank, pred in enumerate(result[1]):
            if rank == result[3][-1]-1:
                print(pred, result[3])
            else:
                print(pred)
        print("#"*len(seq))

        resultsC.append(result[2])
        resultsB.append(result[3])

    tpC = sum(x[0] for x in resultsC)
    fpC = sum(x[1] for x in resultsC)
    fnC = sum(x[2] for x in resultsC)
    fsC = [x[3] for x in resultsC]
    prC = [x[4] for x in resultsC]
    rcC = [x[5] for x in resultsC]

    m1 = round(2*tpC / (2*tpC + fpC + fnC), 3)
    m2 = round(np.mean(fsC), 3)

    print(round(2*m1*m2/(m1+m2),3), round(2*tpC / (2*tpC + fpC + fnC), 3),
          round(np.mean(fsC), 3), round(np.mean(prC), 3), round(np.mean(rcC), 3))

    tpB = sum(x[0] for x in resultsB)
    fpB = sum(x[1] for x in resultsB)
    fnB = sum(x[2] for x in resultsB)
    fsB = [x[3] for x in resultsB]
    prB = [x[4] for x in resultsB]
    rcB = [x[5] for x in resultsB]
    rkB = [x[6] for x in resultsB]

    m1 = round(2*tpB / (2*tpB + fpB + fnB), 3)
    m2 = round(np.mean(fsB), 3)

    
    
    print(round(2*m1*m2/(m1+m2),3), round(2*tpB / (2*tpB + fpB + fnB), 3),
          round(np.mean(fsB), 3), round(np.mean(prB), 3), round(np.mean(rcB), 3))
    print(Counter(rkB))















    
