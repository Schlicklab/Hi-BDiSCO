import pyBigWig
import copy as copy
import math
import numpy as np
import os
import cooler
import random

def nuc_pos(filename, chr, start, end):
	f = open(filename,'r')
	content = [x for x in f]
	f.close()
	data = [x.split() for x in content]
	l = []
	for i in data:
		if (i[0]==chr and int(i[1])>=start and int(i[1])<=end):
			l.append(int(i[1])-start)
	return l

def nuc_choice(l_all, min_ll):
	l_all.sort()
	l=[l_all[0]]
	for i in l_all:
		if (i-l[-1]>min_ll+147):
			l.append(i)
	return l

def linker_length(l_all):
	l = []
	for i in range(1,len(l_all)):
		l.append(l_all[i]-l_all[i-1]-147)
	return l

def ac_peak(file_macs, chr, start, end):
	f=open(file_macs,'r')
	content = [x for x in f]
	f.close()
	data = [x.split() for x in content]
	data_ext = []
	for i in data:
		if (i[0]==chr and ((int(i[1])>=start and int(i[1]) < end) or (int(i[2])>start and int(i[2])<=end))):
			data_ext.append(i)
	return data_ext

def nuc_ac(pos, data_ext, start):
	for i in data_ext:
		if ((pos+start) > int(i[1]) and (pos+start) < int(i[2])):
			return 1
			break
	return 0

def LH_pos(rate, l_LH_values, l_choice_U):
	LH = []
	l_LH_copy = copy.copy(l_LH_values)
	l_LH_copy.sort()
	LH_sort_index = math.floor(rate*len(l_choice_U))
	med_LH = l_LH_copy[LH_sort_index]
	for i in l_LH_values:
		if (i< med_LH):
			LH.append(0)
		else:
			LH.append(1)
	return LH

prefix_f='output/ini_struct/'

f=open('input.txt','r')
c=[x for x in f]
f.close()

d=[x.split() for x in c]

for i in d:
	if i[0]=='chr':
		chr=i[-1]
	elif i[0]=='start':
		start=int(i[-1])
	elif i[0]=='end':
		end=int(i[-1])
	elif i[0]=='nucleosome_position':
		nuc_file=i[-1]
	elif i[0]=='tail_acetylation':
		tail_file=i[-1]
	elif i[0]=='LH':
		LH_file=i[-1]
	elif i[0]=='LH_ratio':
		LH_ratio=float(i[-1])
	elif i[0]=='Hi-C':
		HiC_file=i[-1]
	elif i[0]=='Hi-C-path':
		hicpath=i[-1]
	elif i[0]=='N_rep':
		ncopy=int(i[-1])
	elif i[0]=='N_sim':
		nsimulated=int(i[-1])	

if not os.path.isfile(nuc_file):
	print('Please provide a nucleosome peak file.')
	exit(0)


#####Nucleosome positions

print ('Assigning nucleosome positions...')

l_all = nuc_pos(nuc_file, chr, start, end)
l_choice = nuc_choice(l_all, 20)
ll = linker_length(l_choice)

bead_U = []

for i in ll:
	bead_U.append(round(i*0.34/3-1))


nuc_pos=[]

for i in bead_U:
	if (i>8 and i<36):
		nuc_pos.append(round(i/6)*6)
	elif (i>36):
		nuc_pos.append(36)
	elif (i<2):
		nuc_pos.append(2)
	else:
		nuc_pos.append(i)

nuc_pos.append(2)

f=open(prefix_f+'ll_elig.dat','w')
for i in nuc_pos:
	for j in range(i):
		f.write('0 ')
	f.write('\n')
f.close()

##### Acetylation Island

print ('Assigning histone tail acetylation...')

if not os.path.isfile(tail_file):
	print ('No tail acetylation file provided. Using wildtype tails for all.')
	l_ac=[0]*len(l_choice)
else:	
	data_ext = ac_peak(tail_file, chr, start, end)
	l_ac=[]
	for i in l_choice:
		l_ac.append(nuc_ac(i, data_ext, start))

f=open(prefix_f+'fold_elig.dat','w')
for i in l_ac:
	f.write(str(i)+'\n')

f.close()


##### LH

print ('Assigning linker histones...')

if not os.path.isfile(LH_file):
	print ('No linker histone Chip-seq data provided. Zero LH assigned.')
	LH=[0]*len(l_choice)
else:
	bw_LH = pyBigWig.open(LH_file)
	bw_LH_values = bw_LH.values(chr, start, end)
	l_LH_value = []
	for i in l_choice:
		l_LH_value.append(bw_LH_values[i])
	LH = LH_pos(LH_ratio, l_LH_value, l_choice)

f=open(prefix_f+'LHbound.0.in','w')
f.write('#\n')
f.write('#\n')
for i in LH:
        f.write(str(i)+'\n')
f.close()


################# Generate Input Files 

print ('Generating initial structures...')

S=np.array([[0,0,1],[0,1,0],[1,0,0],[1,0,1],[1,1,1]])

number_of_cores = len(nuc_pos)

z_core_increment_array=[]
core_distance=[]
vector=[]

for i in nuc_pos:
	if i==2:
		z_core_increment_array.append(4.2)
	elif i==3:
		z_core_increment_array.append(4.2)
	elif i==4:
		z_core_increment_array.append(4.2)
	elif i==5:
		z_core_increment_array.append(4.4)
	elif i==6:
		z_core_increment_array.append(5.0)
	elif i==7:
		z_core_increment_array.append(5.1)
	elif i==8:
		z_core_increment_array.append(5.2)
	else:
		z_core_increment_array.append(i*2.5+5)


for i in nuc_pos:
	if i==2:
		core_distance.append(4.7)
	elif i==3:
		core_distance.append(5.5)
	elif i==4:
		core_distance.append(6.5)
	elif i==5:
		core_distance.append(6.7)
	elif i==6:
		core_distance.append(8.5)
	elif i==7:
		core_distance.append(8.75)
	elif i==8:
		core_distance.append(9.0)
	else:
		core_distance.append(7.25)
	z_temp = sum(z_core_increment_array)+1
	vector.append([-core_distance[-1],core_distance[-1],z_temp])

alpha=np.pi*0/80
beta=-np.pi*0/180
gama=np.pi*205/180

totlink=sum(nuc_pos)

cores=np.zeros((number_of_cores,3))
linkers=np.zeros((totlink,3))
linker_ends=np.zeros(((1+number_of_cores)*2,3))
core_orientations=np.zeros((number_of_cores*3,3))
l1=np.zeros((3))
l2=np.zeros((3))

rot_x=np.array([[1,0,0],[0,np.cos(alpha),np.sin(alpha)],[0,-np.sin(alpha),np.cos(alpha)]])
rot_y=np.array([[np.cos(beta),0,-np.sin(beta)],[0,1,0],[np.sin(beta),0,np.cos(beta)]])
rot_z=np.array([[np.cos(gama),np.sin(gama),0],[-np.sin(gama),np.cos(gama),0],[0,0,1]])

#### First core ####
#### coordinates definitions ####
l1[0] = cores[0][0]+5.12
l1[1] = cores[0][1]-3.1
l1[2] = cores[0][2]-1.8
l2[0] = cores[0][0]+3.1
l2[1] = cores[0][1]-5.12
l2[2] = cores[0][2]+1.8

#### movement of linkers to the center
l1_m=l1-cores[0]
l2_m=l2-cores[0]

#### first core rotation
rot_m=np.dot(np.dot(rot_x,rot_y),np.eye(3))

#### first core position
cores[0][0]=vector[0][0]
cores[0][1]=vector[0][1]
cores[0][2]=vector[0][2]

#### first core orientation
core_orientations[0]=rot_m[0]
core_orientations[1]=rot_m[1]
core_orientations[2]=rot_m[2]

l1_m=np.dot(l1_m,rot_m)
l1_m=l1_m+cores[0]

l2_m=np.dot(l2_m,rot_m)
l2_m=l2_m+cores[0]

linker_ends[0]=l1_m
linker_ends[1]=l2_m

i=1
r_x=linker_ends[(i-1)*2+1][0]-linker_ends[(i-1)*2][0]
r_y=linker_ends[(i-1)*2+1][1]-linker_ends[(i-1)*2][1]
r_z=linker_ends[(i-1)*2+1][2]-linker_ends[(i-1)*2][2]

pom_x=r_x/(2*(nuc_pos[i-1]-1))
pom_y=r_y/(2*(nuc_pos[i-1]-1))
pom_z=r_z/(2*(nuc_pos[i-1]-1))

for l in range(nuc_pos[i-1]):
	linkers[l][0]=-pom_x*l+linker_ends[(i-1)*2+1][0]
	linkers[l][1]=-pom_y*l+linker_ends[(i-1)*2+1][0]
	linkers[l][2]=-pom_z*l+linker_ends[(i-1)*2+1][0]

rot_m=np.dot(rot_m,rot_z)

j=1

for i in range(1,number_of_cores):
	core_orientations[i*3]=rot_m[0]
	core_orientations[i*3+1]=rot_m[1]
	core_orientations[i*3+2]=rot_m[2]
	vector=np.dot(vector,rot_z)
	z_temp=z_temp-z_core_increment_array[i-1]
	cores[i][0]=vector[i][0]
	cores[i][1]=vector[i][1]
	cores[i][2]=z_temp
	l1[0]=3.1
	l1[1]=-5.12
	l1[2]=-1.8
	l2[0]=5.12
	l2[1]=-3.1
	l2[2]=1.8
	l1_m=np.dot(l1,rot_m)
	l1_m=l1_m+cores[i]
	l2_m=np.dot(l2,rot_m)
	l2_m=l2_m+cores[i]
	j+=1
	linker_ends[j][0]=l1_m[0]
	linker_ends[j][1]=l1_m[1]
	linker_ends[j][2]=l1_m[2]
	j+=1
	linker_ends[j][0]=l2_m[0]
	linker_ends[j][1]=l2_m[1]
	linker_ends[j][2]=l2_m[2]
	rot_m=np.dot(rot_m,rot_z)

for i in range(number_of_cores):
	r_x=linker_ends[i*2+3][0]-linker_ends[i*2][0]
	r_y=linker_ends[i*2+3][1]-linker_ends[i*2][1]
	r_z=linker_ends[i*2+3][2]-linker_ends[i*2][2]
	pom_x=r_x/(2*(nuc_pos[i]-1))
	pom_y=r_y/(2*(nuc_pos[i]-1))
	pom_z=r_z/(2*(nuc_pos[i]-1))
	for l in range(nuc_pos[i]):
		linkers[sum(nuc_pos[0:i])+l][0]=pom_x*2*l+linker_ends[i*2][0]
		linkers[sum(nuc_pos[0:i])+l][1]=pom_y*2*l+linker_ends[i*2][1]
		linkers[sum(nuc_pos[0:i])+l][2]=pom_z*2*l+linker_ends[i*2][2]


count=0
f=open(prefix_f+'data_mod','w')
for i in range(number_of_cores):
	f.write(str("%18.13f" % cores[i][0]).rjust(18)+str("%18.13f" % cores[i][1]).rjust(20)+str("%18.13f" % cores[i][2]).rjust(20)+'  1\n')
	for j in range(3):
		f.write(str("%18.13f" % core_orientations[i*3+j][0]).rjust(18)+str("%18.13f" % core_orientations[i*3+j][1]).rjust(20)+str("%18.13f" % core_orientations[i*3+j][2]).rjust(20)+'\n')	
	for k in range(nuc_pos[i]):
		f.write(str("%18.13f" % linkers[count][0]).rjust(18)+str("%18.13f" % linkers[count][1]).rjust(20)+str("%18.13f" % linkers[count][2]).rjust(20)+'  0\n')
		count+=1
		for j in range(3):
			f.write(str("%18.13f" % core_orientations[i*3+j][0]).rjust(18)+str("%18.13f" % core_orientations[i*3+j][1]).rjust(20)+str("%18.13f" % core_orientations[i*3+j][2]).rjust(20)+'\n')  
f.close()

f=open(prefix_f+'1.dat','w')
for i in linkers:
	f.write(str("%6.4e" % i[0])+'  '+str("%6.4e" % i[1])+'  '+str("%6.4e" % i[2])+'\n')
for i in range(number_of_cores):
	f.write(str("%6.4e" % cores[i][0])+'  '+str("%6.4e" % cores[i][1])+'  '+str("%6.4e" % cores[i][2])+'\n')
	for j in range(3):
		f.write(str("%6.4e" % core_orientations[i*3+j][0])+'  '+str("%6.4e" % core_orientations[i*3+j][1])+'  '+str("%6.4e" % core_orientations[i*3+j][2])+'\n')

f.close()


f=open(prefix_f+'LH.in','w')
for i in LH:
	f.write(str(i)+'\n')

f.close()

f=open(prefix_f+'dim.in','w')
f.write('#cores\n')
f.write(str(number_of_cores)+'\n')
f.write('#LH\n')
f.write('0 6 22\n')
f.write('#DNA\n')

for i in nuc_pos:
	f.write(str(i)+'\n')	
f.close()

#### Analysis Hi-C data
print ('Analyzing Hi-C data...')

if not os.path.isfile(HiC_file):
	print('Please provide a Hi-C file.')
	exit(0)

beads_num=sum(nuc_pos)+len(nuc_pos)*16

l=[]
count=1

for i in nuc_pos:
	l.extend([count]*16)
	count+=1
	for j in range(i):
		l.append(count)
		count+=1

c=cooler.Cooler(HiC_file+'::'+hicpath)

A=c.matrix().fetch((chr,start,end))

tmax=max([max(x) for x in A])
tmin=min([min(x) for x in A])

bond_selection=[[] for i in range(nsimulated)]

print ('Distributing restraints...')

for i in range(0,len(A)-1):
	for j in range(i+1,len(A)):
		if np.isnan(A[i][j])==False and abs(l[int(i*beads_num/len(A))]-l[int(j*beads_num/len(A))])>1:
			for x in random.sample(range(ncopy),int(ncopy*(1+(ncopy-1)*(A[i][j]-tmin)/(tmax-tmin))/ncopy)):
				if x<nsimulated:
					bond_selection[x].append([l[int(i*beads_num/len(A))],l[int(j*beads_num/len(A))]])

for x in range(1,len(bond_selection)+1):
	os.mkdir(prefix_f+'copy_'+str(x))
	f=open(prefix_f+'copy_'+str(x)+'/ex_force.txt','w')
	for i in bond_selection[x-1]:
		f.write(str(i[0])+' '+str(i[1])+'\n')
	f.close()

print ('Done')
