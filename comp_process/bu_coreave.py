
import mocup

import re
import os

#'''search for volume of 20 zones from input file inp1 and make a list of volume list_vol'''
searchfile = open('inp1', 'r+')
lines = searchfile.readlines() 
output_str = []
list_vol = [[None for _ in range(5)] for _ in range(4)] 
print list_vol
for i in range (1, 5):
	for j in range(1, 6):
		print j
		linenb =0
		string_to_find = 'c Depletion Zones %d%d' %(i, j)
		for k in range(linenb, len(lines)):
			linenb = linenb +1
			if string_to_find in lines[k]:
				list_vol[i-1][j-1]=lines[k+3].split('=')[1]
				output_str.append(str(i)+str(j) + ' ' + lines[k+3]+ '\n')
	outf=open('volume', 'w+')
	outf.write(''.join(output_str))
	print list_vol
	outf.close()
	searchfile.close()

# ''' calculate volume averaged fuel composition over the core, for 8 different depletion passes'''
for B in range(1,9):
	for R in range(1,5):
    		for Z in range(1,6):
        		node_mat = mocup.material()
        		print dir(node_mat)
            		# mat_loc = 'moi_files/moi.%d%d%d00.eq.pch' % (R,Z,B)
            		mat_loc = 'fuel_comp/m%d%d%d00' % (R,Z,B)
            		mat = mocup.material()
            		mat.import_ocf(mat_loc)
            		node_mat = node_mat + mat * float(list_vol[R-1][Z-1])
	mat_loc = 'bu_files/vol_ave_from_file/bu%d.eq' % (B)
	node_mat.make_ocf(mat_loc,lib=node_mat.lib)
