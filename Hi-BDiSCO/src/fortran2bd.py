f = open('out.fort','r')
content = [x.rstrip('\n') for x in f]
f.close()

data = [x.split() for x in content]

f = open('dim.in','r')
content = [x.rstrip('\n') for x in f]
f.close()

struc = [x.split() for x in content]

f = open('LHbound.0.out','r')
content = [x.rstrip('\n') for x in f]
f.close()

lhstruc = [x.split() for x in content[2:]]


f1 = open('restart_r.txt','w')
f2 = open('restart_a.txt','w')
f3 = open('restart_b.txt','w')
f4 = open('restart_c.txt','w')


nuc_num = int(struc[1][0])

dna_start = nuc_num*4

dna_tot = 0


k=0

for i in range(nuc_num):
    f1.write(str("%0.15f" % float(data[i*4][0]))+'\n')
    f1.write(str("%0.15f" % float(data[i*4][1]))+'\n')
    f1.write(str("%0.15f" % float(data[i*4][2]))+'\n')
    f2.write(str("%0.15f" % float(data[i*4+1][0]))+'\n')
    f2.write(str("%0.15f" % float(data[i*4+1][1]))+'\n')
    f2.write(str("%0.15f" % float(data[i*4+1][2]))+'\n')
    f3.write(str("%0.15f" % float(data[i*4+2][0]))+'\n')
    f3.write(str("%0.15f" % float(data[i*4+2][1]))+'\n')
    f3.write(str("%0.15f" % float(data[i*4+2][2]))+'\n')
    f4.write(str("%0.15f" % float(data[i*4+3][0]))+'\n')
    f4.write(str("%0.15f" % float(data[i*4+3][1]))+'\n')
    f4.write(str("%0.15f" % float(data[i*4+3][2]))+'\n')
    dna_tot += int(struc[5+i][0])
    for j in range(int(struc[5+i][0])):
        f1.write(str("%0.15f" % float(data[dna_start+k*4][0]))+'\n')
        f1.write(str("%0.15f" % float(data[dna_start+k*4][1]))+'\n')
        f1.write(str("%0.15f" % float(data[dna_start+k*4][2]))+'\n')
        f2.write(str("%0.15f" % float(data[dna_start+k*4+1][0]))+'\n')
        f2.write(str("%0.15f" % float(data[dna_start+k*4+1][1]))+'\n')
        f2.write(str("%0.15f" % float(data[dna_start+k*4+1][2]))+'\n')
        f3.write(str("%0.15f" % float(data[dna_start+k*4+2][0]))+'\n')
        f3.write(str("%0.15f" % float(data[dna_start+k*4+2][1]))+'\n')
        f3.write(str("%0.15f" % float(data[dna_start+k*4+2][2]))+'\n')
        f4.write(str("%0.15f" % float(data[dna_start+k*4+3][0]))+'\n')
        f4.write(str("%0.15f" % float(data[dna_start+k*4+3][1]))+'\n')
        f4.write(str("%0.15f" % float(data[dna_start+k*4+3][2]))+'\n') 
        k=k+1
    
    


f1.close()
f2.close()
f3.close()
f4.close()




