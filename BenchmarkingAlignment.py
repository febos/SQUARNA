import os, time, glob

from SQUARNA import ReadStockholm
from SQRNdbnseq import PairsToDBN, DBNToPairs


def PredictSQUARNAs1(dataset, fam):

    command = "python SQUARNA.py i={} a step3=1 > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs2(dataset, fam):

    command = "python SQUARNA.py i={} a step3=2 > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs3i(dataset, fam):

    command = "python SQUARNA.py i={} a step3=i > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictSQUARNAs3u(dataset, fam):

    command = "python SQUARNA.py i={} a step3=u > outp3.tmp".format("datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictRNAalifold(dataset, fam):

    command = "RNAalifold --noPS {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictCentroidAlifold(dataset, fam):

    command = "~/software/centroid-rna-package-master/build/src/centroid_alifold"+\
              " {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].split()[0]


def PredictIPknot(dataset, fam):

    command = "~/software/ipknot-master/build/ipknot"+\
              " {} > outp3.tmp".format("datasets/{}/aln/{}.aln".format(dataset,fam))
    os.system(command)
    with open("outp3.tmp") as file:
        lines = file.readlines()
        return lines[-1].strip()


def PredictRscapeNested(dataset, fam):

    command = "cd tmp; ~/software/rscape/rscape_v1.6.1/bin/R-scape"+\
              " --fold --covmin 4 --nofigures --rna {} > outp3.tmp".format("../datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("tmp/outp3.tmp") as file:

        names = []
        dbns  = {}
        flag = False
        
        for line in file:
            if line.startswith("# The predicted CaCoFold structure"):
                flag = True
            elif flag and line.startswith("# SS_cons"):
                _, name, dbn = line.strip().split()
                if name not in dbns:
                    names.append(name)
                    dbns[name] = ''
                dbns[name] += dbn
    return PairsToDBN(DBNToPairs(dbns["SS_cons"]),len(dbns["SS_cons"]))             
                 

def PredictRscapeTotal(dataset, fam):

    command = "cd tmp; ~/software/rscape/rscape_v1.6.1/bin/R-scape"+\
              " --fold --covmin 4 --nofigures --rna {} > outp3.tmp".format("../datasets/{}/sto/{}.sto".format(dataset,fam))
    os.system(command)
    with open("tmp/outp3.tmp") as file:

        names = []
        dbns  = {}
        flag = False
        
        for line in file:
            if line.startswith("# The predicted CaCoFold structure"):
                flag = True
            elif flag and line.startswith("# SS_cons"):
                _, name, dbn = line.strip().split()
                if name not in dbns:
                    names.append(name)
                    dbns[name] = ''
                dbns[name] += dbn

    pairs = DBNToPairs(dbns["SS_cons"])
    seen = set(pos for bp in pairs for pos in bp)

    for name in names:
        if name == "SS_cons":
            continue
        for v, w in DBNToPairs(dbns[name]):
            if v not in seen and w not in seen:
                seen.add(v)
                seen.add(w)
                pairs.append((v, w))
    return PairsToDBN(sorted(pairs), len(dbns["SS_cons"])) 


if __name__ == "__main__":

    #dataset = "SubAli" # RNAStralignExt / Rfam14.9 / RfamPDB / SubAli / SeqSim
    #tool    = "IPknot"

    for dataset, tool in (("RNAStralignExt","SQUARNAs1"),
                          ("RNAStralignExt","SQUARNAs2"),
                          ("RNAStralignExt","SQUARNAs3i"),
                          ("RNAStralignExt","SQUARNAs3u"),
                          ("Rfam14.9","SQUARNAs1"),
                          ("Rfam14.9","SQUARNAs2"),
                          ("Rfam14.9","SQUARNAs3i"),
                          ("Rfam14.9","SQUARNAs3u"),
                          ("RfamPDB","SQUARNAs1"),
                          ("RfamPDB","SQUARNAs2"),
                          ("RfamPDB","SQUARNAs3i"),
                          ("RfamPDB","SQUARNAs3u"),
                          ("SubAli","SQUARNAs1"),
                          ("SubAli","SQUARNAs2"),
                          ("SubAli","SQUARNAs3i"),
                          ("SubAli","SQUARNAs3u"),
                          ("SeqSim","SQUARNAs1"),
                          ("SeqSim","SQUARNAs2"),
                          ("SeqSim","SQUARNAs3i"),
                          ("SeqSim","SQUARNAs3u"),
                          ("S01AliCMclean","SQUARNAs1"),
                          ("S01AliCMclean","SQUARNAs2"),
                          ("S01AliCMclean","SQUARNAs3i"),
                          ("S01AliCMclean","SQUARNAs3u"),
                          ("S01AliUngapclean","SQUARNAs1"),
                          ("S01AliUngapclean","SQUARNAs2"),
                          ("S01AliUngapclean","SQUARNAs3i"),
                          ("S01AliUngapclean","SQUARNAs3u"),
                          ("S01Aliclean","SQUARNAs1"),
                          ("S01Aliclean","SQUARNAs2"),
                          ("S01Aliclean","SQUARNAs3i"),
                          ("S01Aliclean","SQUARNAs3u"),):
                
        outname = "{}_{}".format(dataset,tool)
        title = '\t'.join("NAME LEN DEPTH TIME TP FP FN PRC RCL FS DBN PRED".split())
        outp1 = open(outname+'.fas','w')
        outp2 = open(outname+'.tsv','w')
        outp2.write(title+'\n')

        t0 = time.time()

        famfiles = glob.glob("datasets/{}/sto/*".format(dataset))
        
        fams = []

        for famfile in famfiles:
            fam = os.path.basename(famfile).split('.')[0]
            headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm(famfile)
            LEN    = len(gcdict['SS_cons'])
            DEPTH  = len(seqnames)
            refdbn = PairsToDBN(DBNToPairs(gcdict['SS_cons']), LEN)
            fams.append((LEN, DEPTH, fam, refdbn))
        fams.sort()

        cnt = 0  
        for LEN, DEPTH, fam, refdbn in fams:
            cnt += 1
            name = '>'+fam
            print(name,end='')

            preddbn = {"SQUARNAs1":PredictSQUARNAs1,
                        "SQUARNAs2":PredictSQUARNAs2,
                        "SQUARNAs3i":PredictSQUARNAs3i,
                        "SQUARNAs3u":PredictSQUARNAs3u,
                        "RNAalifold":PredictRNAalifold,
                        "CentroidAlifold":PredictCentroidAlifold,
                        "IPknot": PredictIPknot,
                        "RscapeNested": PredictRscapeNested,
                        "RscapeTotal" : PredictRscapeTotal,
                        }[tool](dataset, fam)

            t1 = time.time()-t0

            print("...COMPLETE ({}sec) == {}/{}".format(round(t1,3), cnt, len(fams)))

            pairsr = set(DBNToPairs(refdbn))
            pairsq = set(DBNToPairs(preddbn))

            TP = len(pairsr & pairsq)
            FP = len(pairsq - pairsr)
            FN = len(pairsr - pairsq)
                    
            FS = 2*TP / (2*TP + FN + FP) if (TP + FN + FP) else 1
            PRC = (TP / (TP + FP)) if (TP+FP) else 1
            RCL = (TP / (TP + FN)) if (TP+FN) else 1

            outp1.write(name+'\n')
            outp1.write(refdbn+'\n')
            outp1.write(preddbn+'\n')
            outp1.write("LEN={} DEPTH={}, TIME={}sec TP={} FP={} FN={} PRC={} RCL={} FS={}\n"\
                        .format(LEN, DEPTH,round(t1,3),
                                TP,FP,FN,
                                round(PRC,3),
                                round(RCL,3),
                                round(FS,3)))
            res = [name[1:], LEN, DEPTH, round(t1,3), TP, FP, FN,
                    round(PRC,3), round(RCL,3), round(FS,3),
                    refdbn,preddbn]
            outp2.write('\t'.join([str(g) for g in res])+'\n')

        outp1.close()
        outp2.close()

