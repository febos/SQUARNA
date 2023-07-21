import os

fams = '''RF00001
RF00002
RF00003
RF00004
RF00005
RF00007
RF00008
RF00009
RF00010
RF00011
RF00012
RF00013
RF00015
RF00017
RF00020
RF00023
RF00024
RF00025
RF00026
RF00027
RF00028
RF00029
RF00030
RF00031
RF00032
RF00037
RF00044
RF00050
RF00059
RF00061
RF00066
RF00075
RF00080
RF00083
RF00100
RF00102
RF00106
RF00114
RF00162
RF00164
RF00166
RF00167
RF00168
RF00169
RF00174
RF00175
RF00177
RF00207
RF00209
RF00210
RF00228
RF00233
RF00234
RF00240
RF00250
RF00254
RF00373
RF00374
RF00375
RF00379
RF00380
RF00382
RF00386
RF00390
RF00442
RF00458
RF00480
RF00488
RF00500
RF00504
RF00505
RF00507
RF00522
RF00525
RF00619
RF00622
RF00634
RF00957
RF01047
RF01051
RF01054
RF01073
RF01084
RF01330
RF01344
RF01357
RF01380
RF01381
RF01415
RF01689
RF01704
RF01725
RF01727
RF01734
RF01739
RF01750
RF01763
RF01786
RF01807
RF01826
RF01831
RF01846
RF01852
RF01854
RF01857
RF01959
RF01960
RF01998
RF02001
RF02012
RF02033
RF02095
RF02253
RF02340
RF02348
RF02359
RF02519
RF02540
RF02541
RF02542
RF02543
RF02545
RF02546
RF02553
RF02678
RF02679
RF02680
RF02681
RF02683
RF02796
RF03013
RF03054
RF04190
RF04222'''.split('\n')

print(len(fams))

for fam in fams:
    for fmt in ("afa","aln","sto"):
        file1 = "Rfam14.9/{}/{}.{}".format(fmt,fam,fmt)
        file2 = file1.replace("Rfam14.9","RfamPDB")
        os.system("cp {} {}".format(file1, file2))
        