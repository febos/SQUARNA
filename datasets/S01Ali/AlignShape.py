import glob, os


files = glob.glob("afa/*")

for file in files:

    fam = os.path.basename(file)[:-4]
    print(fam)
    with open(file) as inp:
        seq = inp.readlines()[4].strip()

    shapefile = 'shape2.0_unaligned/{}.dat'.format(fam)

    with open(shapefile) as inp:
        reacts = [float(x.split()[1]) for x in inp.readlines() if x]

    alireacts = []
    print(seq)
    cur = 0
    for k,ch in enumerate(seq):
        if ch != '-':
            alireacts.append(reacts[cur])
            cur += 1
        else:
            alireacts.append(-999)


    newshapefile = shapefile.replace('unaligned','aligned')

    with open(newshapefile,'w') as outp:
        for k,v in enumerate(alireacts):
            outp.write("{} {}".format(k+1, round(v,4) if v!=int(v) else int(v))+'\n')
