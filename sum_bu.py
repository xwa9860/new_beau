
import mocup

tot_mat = mocup.material()
for R in range(1,5):
    for Z in range(1,6):
        print dir(tot_mat)

        mat_loc = 'bu_files/bu9.%d%d.eq' % (R,Z)
        mat = mocup.material()
        mat.import_ocf(mat_loc)
        tot_mat = tot_mat + mat

mat_loc='bu_files/tot.eq' 
tot_mat.make_ocf(mat_loc,lib=tot_mat.lib)
