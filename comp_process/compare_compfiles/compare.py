#!/usr/bin/python

import mocup
import material

mat=material.material()
mat.import_ocf('moi.11100.eq.pch')
mat.mcf(1, mcf_loc = 'moi.11100.eq.mcnpformat')

