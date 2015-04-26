#!/usr/bin/python

import os

for p in range(1, 9):
	for R in range(1,5):
		for Z in range(1,6):
			mat_name = 'm%d%d%d00' %(R, Z, p)
			mat_comp=[]
			nb_isotope =0
			with open('mk1.txt', 'r+') as f:
				for line in f:
					if mat_name in line:
						#mat_comp.append(line)
						ne = f.next()
						while 'm' not in ne:
							mat_comp.append(ne.split()[0]+' %.8e' %(float(ne.split()[1])+float(ne.split()[3]) )+ '\n')
							ne = f.next()
							nb_isotope = nb_isotope + 1
						print nb_isotope
				with open(mat_name, 'w+') as outputf:
					outputf.write(''.join(mat_comp))
				outputf.close()
f.close()
