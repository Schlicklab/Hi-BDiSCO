f=open('out_f.xyz','r')
content=[x for x in f]
f.close()
data = [x.split() for x in content]

core_data=[]
dna_data=[]
abc_data=[]

flag=0

for i in range(len(data)):
    if data[i][0]=='OC':
        core_data.append(data[i][1:])
        flag=1
    elif data[i][0]=='CA':
        dna_data.append(data[i][1:])
    elif data[i][0]=='H1' and flag==1:
        abc_data.append(data[i][1:])
    elif data[i][0]=='H2' and flag==1:
        abc_data.append(data[i][1:])
    elif data[i][0]=='H3' and flag==1:
        abc_data.append(data[i][1:])
        flag=0


f=open('1.dat','w')

for i in dna_data:
    for j in range(3):
        f.write(str("%6.4e" % float(i[j]))+' ')
    f.write('\n')

for i in range(len(core_data)):
    for j in range(3):
        f.write(str("%6.4e" % float(core_data[i][j]))+' ')
    f.write('\n')
    for x in range(3):
        for y in range(3):
            f.write(str("%6.4e" % float(abc_data[i*3+x][y]))+' ')
        f.write('\n')

f.close()

