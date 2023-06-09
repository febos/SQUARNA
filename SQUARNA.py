import numpy as np

from SQRNdbnseq import SQRNdbnseq, BPMatrix, DBNToPairs, AnnotateStems


if __name__ == "__main__":

    from collections import Counter

    rst = None

    #SAM riboswitch
    seq = "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAACA"
    dbn = "(((((((((....(((((...(((.[[[[)))......)))))(((..(((((...(((((....))))).)))..)).)))...(]]]](((((.......)))))..))))))))))." 
    #TRNA
    #seq = "GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA"
    #dbn = "(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))...."
    #Twister ribozyme
    #seq = "GAAAUAAUGCUACCAGACUGCACAAAAAGUCUGCCCGUUCCAAGUCGGGAUGAAAGCAGAGGAUUUC"
    #dbn = "((((...(((({{((((((........))))))(((..[[[..}}.))).....))))..]]]))))"
    #rst = "X_XX___....__..............,,,,,,.........................xx......x"
    #rst = "..................................................................." 
    
    queue = []

    with open("SRtrain150.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-2,3):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            db = lines[ii+2].strip()
            queue.append([nm, sq, db, rst])


    #seq = "AAACCACGAGGAAGAGAGGUAGCGUUUUCUCCUGAGCGUGAAGCCGGCUUUCUGGCGUUGCUUGGCUGCAACUGCCGUCAGCCAUUGAUGAUCGUUCUUCUCUCCGUAUUGGGGAGUGAGAGGGAGAGAACGCGGUCUGAGUGGU"
    #dbn = "..(((((......(.......((((((((((((....(...(...((((..(.((((((((......))))).))))..))))..)......)..((((((((((.....)))))).))))))))))))))))...)...)))))"

    #seq = "GGUAAAGAAUGAAAAAACACGAUUCGGUUGGUAGUCCGGAUGCAUGAUUGAGAAUGUCAGUAACCUUCCCCUCCUCGGGAUGUCCAUCAUUCUUUAAUAUCUUUUAUGAGGAGGGAA"
    #dbn = None

    #seq = "GGACUUAUAGAUGGCUAAAAUCUGAGUCCA"
    #dbn = "((((((..((((.......))))))))))."
    
    #seq  = "CGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGCGGAAACGCUGAUGUCGUAUACAGGGCU"
    #dbn  = "(((((((((...[[[[[..)))))))))((((((((((((....))))).)))))))....]]]]]"
    
    #seq = "GGGACCAGUUGAACCUGAACAGGGUAAUGCCUGCGCAGGGAGGGUGCUUGUUCACAGGCUGAGAAAGUCCCUGUGUC"
    #dbn = None
    #rst = "(..........................................................................)."
    
    #queue  = [["default", seq, dbn, rst],]
    #queue += [["default", seq, dbn, rst],]  

    poor = {6, 28, 54, 55, 61, 73, 81, 90, 129, 133, 144, 150, 155, 158, 160, 161, 162, 163, 165, 168, 173, 174, 175, 177, 179, 182, 187, 191, 200, 202, 203, 204, 205, 206, 209, 211, 215, 223, 224, 226, 228, 232, 238, 239, 244, 245, 247, 252, 254, 258, 259, 260, 265}

    #queue = [queue[i] for i in range(len(queue)) if i in poor]

    paramsets = []

    #NN = 264
    #queue = queue[NN:NN+1]

    """ TOP SINGLE 
    paramsets.append({"bpweights" : {'GU' : -1.25,
                           'AU' :  1.25,
                           'GC' :  3.25,},
            "suboptmax" : 0.9,
            "suboptmin" : 0.65,
            "suboptsteps": 1,
            "minlen" : 2,
            "minbpscore" : 4.5,
            "minfinscorefactor" : 1.0,
            "distcoef" : 0.09,
            "bracketweight" :  -2.0,##
            "orderpenalty"  : 1.00,
            "loopbonus": 0.125,
            "maxstemnum" : 10**6,
            "mode": "diffedge",##
           })  """ 

    """ TOP ONE  """
    paramsets.append({"bpweights" : {'GU' : -1.25,
                           'AU' :  1.25,
                           'GC' :  3.25,},
            "suboptmax" : 0.9,
            "suboptmin" : 0.65,
            "suboptsteps": 1,
            "minlen" : 2,
            "minbpscore" : 4.5,
            "minfinscorefactor" : 1.25,
            "distcoef" : 0.09,
            "bracketweight" :  -2.0,##
            "orderpenalty"  : 1.00,
            "loopbonus": 0.125,
            "maxstemnum" : 10**6,
            "mode": "diffedge",##
           })  
    

    """ TOP TWO """
    paramsets.append({"bpweights" : {'GU' : 1,
                           'AU' :  1,
                           'GC' :  2,},
            "suboptmax" : 0.9,
            "suboptmin" : 0.65,
            "suboptsteps": 1,
            "minlen" : 2,
            "minbpscore" : 3,
            "minfinscorefactor" : 0.99,
            "distcoef" : 0.1,
            "bracketweight" :  -2.0,##
            "orderpenalty"  : 1.35,
            "loopbonus": 0.125,
            "maxstemnum" : 10**6,
            "mode": "diffedge",##
           })    

    threads = 4
    
    toplim     = 5
    conslim    = 1
    hardrest   = False
    rankbydiff = False
    
    resultsB = []
    resultsC = []

    for obj in queue:

        name, seq, dbn, rst = obj

        result = SQRNdbnseq(seq, rst, dbn,
                            paramsets, conslim, toplim,
                            hardrest, rankbydiff,
                            threads)

        print(name)
        print(seq)
        if dbn:
            print(dbn)
        print('_'*len(seq))
        print(result[0],result[2])
        print('='*len(seq))
        for rank, pred in enumerate(result[1][:toplim]):
            if rank == result[3][-1]-1:
                print(' '.join([str(gg) for gg in pred]), result[3])
            else:
                print(' '.join([str(gg) for gg in pred]))
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

    #print([i for i in range(len(fsB)) if fsB[i]<0.7]) # poorly predicted sequences

    m1 = round(2*tpB / (2*tpB + fpB + fnB), 3)
    m2 = round(np.mean(fsB), 3)

    
    
    print(round(2*m1*m2/(m1+m2),3), round(2*tpB / (2*tpB + fpB + fnB), 3),
          round(np.mean(fsB), 3), round(np.mean(prB), 3), round(np.mean(rcB), 3))
    print(Counter(rkB))

    #for x in fsB:
    #    print(x)















    
