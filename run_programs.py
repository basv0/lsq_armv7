# -*- coding: utf-8 -*-

import os
def avg(list):
    sum = 0
    for elm in list:
        sum += elm
    av=sum/len(list)
    return av
    
R=5 # number of repeatitions

#M = [2,3,4,5,6,7,8,9,10] # object order
M = [3]

nm = ['gsl','double','optim','neonh','neonv']

for nam in nm:
    exe='mnk_test_'+nam
    f=nam+'.txt'
    os.system('rm -f '+f)
    for m in M:
        cmd='./'+exe+' '+ str(m)+' ./ident3.dat >> '+f
        for r in range (0,R):
	    print cmd
            os.system(cmd)
    g=open(f,'r')
    h=open('plot_'+f,'w')
    rs=[]
    for m in M:
        tm=[]
        for r in range(0,R):
            lin=g.readline()
            vals=lin.split()
            tm.append(float(vals[2]))
        tm.remove(min(tm))
        tm.remove(max(tm))
        tim=avg(tm)
        rl=vals[0]+' '+vals[1]+' '+str(tim)+'\n'
        h.write(rl)
    h.close()
    g.close()
    