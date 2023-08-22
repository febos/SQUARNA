import glob, os


files = glob.glob("afa/*")

for file in files:

    outpfile = file.replace("afa/","for_shapesorter/")

    with open(outpfile,'w') as outp:
        
        with open(file) as inp:
            lines = inp.readlines()
            seq = lines[4].strip()
            gaps = {i for i in range(len(seq)) if seq[i] == '-'}

            for line in lines[3:]:
                if line.startswith('>'):
                    outp.write(line)
                else:
                    outp.write(''.join([line[i]
                                        for i in range(len(line))
                                        if i not in gaps]))
