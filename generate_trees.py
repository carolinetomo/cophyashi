import p4
import sys

def f5(seq, idfun=None):
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

p4.read(sys.argv[1])
t = p4.var.trees[0]
di = []
alt = True
alt2 = False
for i in range(5):
    d = t.dupe()
    if alt == True:
        d.nni()
        alt = False
    else:
        if alt2 == False:
            d.randomSpr()
            alt2 = True
        else:
            d.nni()
            alt2 = False
        alt = True
    di.append(d.writeNewick(toString=True,spaceAfterComma=False))
x = f5(di)

for i in range(0,len(x)):
    outfl = open(sys.argv[1]+str(i),"w")
    outfl.write(x[i].strip())
    outfl.close()

