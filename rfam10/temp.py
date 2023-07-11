
def DBNToPairs(dbn):
    """Convert the dbn string into a sorted list of base pairs"""
    pairs = set()

    # keys == closing brackets, values == matching opening brackets
    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',}
    # 30 bp stacks for 30 allowed pseudoknot levels
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],}

    for i,v in enumerate(dbn):
        # if we observe an opening bracket
        # then add its index into the matching stack
        if v in stack: 
            stack[v].append(i)
        # else if we observe the closing bracket
        # take the opening index from the matching stack
        # and add the base pair to the pairs set
        elif v in closing:
            # this is to handle closing brackets with no
            # opening partner - they will be ignored
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(), i))

    return sorted(pairs)


with open("Rscape.txt") as inp:
    lines = inp.readlines()

    for i in range(0,len(lines)-2,3):
        name = lines[i].strip()
        ref = lines[i+1].strip()
        pred = lines[i+2].strip()

        rb = set(DBNToPairs(ref))
        pb = set(DBNToPairs(pred))

        TP = len(pb & rb)
        FP = len(pb - rb)
        FN = len(rb - pb)

        FS = 2*TP / (2*TP + FP + FN) if 2*TP + FP + FN else 1.0
        print(name)
        print(round(FS,3))
