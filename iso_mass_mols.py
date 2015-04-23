import mocup

mpo = mocup.mpo()
mpo.populate(timestep=1)

mass = 0
moles = 0
for cell in mpo.cell:
        moles += mpo.volume[cell]*mpo.concentration[cell]['551370']/.6022
        mass += mpo.volume[cell]*mpo.concentration[cell]['551370']/.6022*137.

print 'mass (g): ', mass
print 'moles: ', moles
