import glob
files = glob.glob("SeqSim/afa/*")
len(files)

files.sort()
import os
cnt = 0
for k,file in enumerate(files):
    if k==0 or file==files[-1]:
        cnt += 1
        print(file)
    else:
        if os.path.basename(file).split('_')[:2] != os.path.basename(files[k+1]).split('_')[:2]:
            cnt += 1
            print(file)
            print(cnt)
            cnt = 0
            print(files[k+1])
        else:
            cnt += 1
print(cnt)
