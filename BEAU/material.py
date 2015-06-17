#!/usr/bin/env python
import depletion


class material:

    def __init__(self):
        self.comp = {}
        self.lib = 'pwru50'
        self.tmp = 0
        self.scat = ''
        # in this dictionary the keys are a string of the nuclide ID and the
        # values are the compositions in moles

    def __add__(self, other):
        # this method adds to materials together

        new = material()
        for nuclide in self.comp.keys():
            if nuclide in other.comp:
                new.comp[nuclide] = self.comp[nuclide] + other.comp[nuclide]
            else:
                new.comp[nuclide] = self.comp[nuclide]
        for nuclide in other.comp.keys():
            if nuclide in self.comp:
                # do nothing these nuclides have already been summed
                pass
            else:
                new.comp[nuclide] = other.comp[nuclide]

        new.tmp = self.tmp*self.moles() + other.tmp*other.moles()

        return new

    def __mul__(self, other):
        new = material()
        new.tmp = self.tmp
        for nuclide in self.comp.keys():
            new.comp[nuclide] = other*self.comp[nuclide]
        return new

    def make_ocf(self, ocf_loc, lib='None'):
        # This method makes and origen compsition file from the material

        LIB = lib
        self.ORIGEN_nucl_list(lib=LIB)

        STRUCTURAL = ''
        ACTINIDES = ''
        FISSIONPRODUCTS = ''

        for nuclide in self.comp.keys():
            if nuclide in self.structonly_nucl:
                STRUCTURAL += '1 %s %1.4E 0 0 0 0 0 0 \n' % (nuclide,
                                                             self.comp[nuclide])
            elif nuclide in self.actin_nucl:
                ACTINIDES += '2 %s %1.4E 0 0 0 0 0 0 \n' % (nuclide,
                                                            self.comp[nuclide])
            elif nuclide in self.fisprod_nucl:
                FISSIONPRODUCTS += '3 %s %1.4E 0 0 0 0 0 0 \n' % (
                    nuclide, self.comp[nuclide])
        open(ocf_loc, 'w').write(STRUCTURAL+ACTINIDES+FISSIONPRODUCTS+'0 0 0 0')

    def import_ocf(self, ocf_loc):
        # This method appends to a material object based on the contents of the
        # file defined

        import re

        self.comp = {}

        pairs = re.compile(
            '\s\d{5,6}\s+\d[.]\d{4,4}[eE+-]{2,2}\d{2,2}').findall(open(ocf_loc).read())

        for pair in pairs:

            self.comp[pair.split()[0]] = float(pair.split()[1])

    def import_mcf(self, inp_loc, mat_ID):
        # This method imports a material vector from MCNP and places it in the
        # materials composition vector

        import re

        inp = open(inp_loc).read()
        token = 'm%d ' % int(mat_ID)

        k = inp.find(token)

        if k == -1:
            print token
            print 'this material was not in the MCNP input deck: %s' % inp_loc

        else:
            # find end of material vector
            l = len(inp) - k
            for card in ['\nm', '\nM', '\nf', '\nF', '\nk', '\nK', '\ns', '\nS']:
                m = inp[k:].find(card)
                if m != -1 and m < l:
                    l = m

            self.mat_card = inp[k:k+l]

            if '$ tmp=' in self.mat_card:
                k = self.mat_card.find('$ tmp=')
                self.tmp = float(
                    self.mat_card[
                        k:].split()[1].replace(
                        'tmp=',
                        '').replace(
                        'K',
                        ''))

            mcnp_nucls = re.compile(
                '\s+\d{4,5}[.]\d{2,2}[a-zA-Z]\s+').findall(self.mat_card)
            mcnp_nucls += re.compile('\s+\d{4,5}\s+').findall(self.mat_card)

            self.comp = {}

            for mcnp_nucl in mcnp_nucls:

                # find the composition

                k = self.mat_card.find(mcnp_nucl)

                comp = float(self.mat_card[k:].split()[1])

                # transform mcnp nuclide to origen nuclide

                mcnp_nucl = mcnp_nucl.replace('.', ' ').split()[0]

                if mcnp_nucl == '6000':
                    nucl = '60120'
                elif mcnp_nucl == '14000':
                    nucl = '140280'
                elif mcnp_nucl == '95601':
                    nucl = '952421'
                elif mcnp_nucl == '52601':
                    nucl = '521251'
                elif mcnp_nucl == '61601':
                    nucl = '621471'
                elif mcnp_nucl[-3:] == '601':
                    print mcnp_nucl
                    print 'add this nuclide to mocup.py import_ method'
                    poop
                elif int(mcnp_nucl[-3:]) > 299:
                    # this is a metastable isotope
                    aw = int(mcnp_nucl[-3:]) - 400
                    nucl = mcnp_nucl[:-3] + '0'*(len(str(aw))-3) + str(aw) + '1'
                else:
                    # groundstate isotope
                    nucl = mcnp_nucl + '0'

                if nucl in self.comp.keys():
                    self.comp[nucl] += comp
                else:
                    self.comp[nucl] = comp

    def import_scf(self, scf_loc):

        self.comp = {}

        mat_card = open(scf_loc).read().split()

        a = depletion.Depletion()
        nct = a.nct

        while len(mat_card) >= 2:
            if mat_card[0] in nct.keys():
                # This value is a nuclide ID

                # enter this value into the material composition vecor
                self.comp[nct[mat_card[0]]] = float(mat_card[1])

                # advance the card
                mat_card = mat_card[2:]
            else:
                mat_card = mat_card[1:]

    def normalize_at(self):
        self.normal_at_comp = {}
        for nuclide in self.comp.keys():
            self.normal_at_comp[nuclide] = self.comp[nuclide]/self.moles()

    def normalize_mass(self):
        self.normal_mass_comp = {}
        for nuclide in self.comp.keys():
            self.normal_mass_comp[nuclide] = self.comp[
                nuclide]*float(nuclide[-4:-1])/self.mass()

    def heavy_mass(self):
        mass = 0
        for nuclide in self.comp.keys():
            AW = float(nuclide[-4:-1])
            Z = float(nuclide[:-4])
            if Z > 80.:
                mass += AW*self.comp[nuclide]
        return mass

    def heavy_moles(self):
        moles = 0
        for nuclide in self.comp.keys():
            Z = float(nuclide[:-4])
            if Z > 80.:
                moles += self.comp[nuclide]
        return moles

    def mass(self):
        mass = 0
        for nuclide in self.comp.keys():
            AW = float(nuclide[-4:-1])
            mass += AW*self.comp[nuclide]
        return mass

    def moles(self):
        moles = 0
        for nuclide in self.comp.keys():
            moles += self.comp[nuclide]
        return moles

    def AW(self):
        AW = self.mass()/self.moles()
        return AW

    def mcf(self, m, xs='72c', mcf_loc='', tmp=0, scat=''):

        # This method makes an string of an MCNP5 formated material definition
        # m is the material number and should be an integer
        # xs is the cross section library and is set to 72c by default
        # mcf_loc is the location where this material card is to be written, if
        # it is left blank it will not be written

        mcnp5_mat = 'm%d $' % int(m)

        if tmp == 0 and self.tmp != 0:
            tmp = self.tmp

        if tmp != 0:
            xstmp = {
                293.4469961: '70c',
                599.6752494: '71c',
                899.512874: '72c',
                1199.373695: '73c',
                2498.608212: '74c'}
            mcnp5_mat += ' tmp=%1.5EK' % tmp
            scattmp = {}
            scattmp['grph'] = {
                293.4469961: '10t',
                399.7912317: '11t',
                499.7332405: '12t',
                599.6752494: '13t',
                699.6288564: '14t',
                799.5708652: '15t',
                999.4664811: '16t',
                1199.373695: '17t',
                1599.164927: '18t',
                1998.956159: '19t'}

            if scat == '':
                scat = self.scat

        self.ORIGEN_nucl_list()

        a = depletion.Depletion()
        nct = a.nct

        metastable = {}

        for nucl in nct.keys():
            if nct[nucl][-1] == '1' and (nct[nucl] not in metastable):
                if nct[nucl] == '952421':
                    # recall in mcnp metastable naming convention is reversed
                    # for am242
                    metastable['952420'] = nct[nucl][:-4] + \
                        str(int(nct[nucl][-4:-1])+400) + '.00c'
                else:
                    metastable[nct[nucl]] = nct[nucl][
                        :-4] + str(int(nct[nucl][-4:-1])+400) + '.00c'

        not_in = []

        for nuclide in self.comp.keys():
            if nuclide not in self.structonly_nucl + self.actin_nucl + self.fisprod_nucl:
                not_in.append(nuclide)

        not_in = sorted(not_in)

        if scat != '' and tmp == 0:
            print '''
                error! mcnp5 composition file generator will crash if you use
                scattering modification without defining a temperature\n
                please modify the material object by either\n\t
                A: defining a realistic temperature for this material (preferred), or \n\t
                B: removing the scat kernel flad (mat.scar) or mcf(m,scat='flag)
                '''

        if tmp != 0:

            # find upper and lower limits
            above = {}
            below = {}
            for t in xstmp.keys():
                if (t - tmp) > 0:
                    above[t - tmp] = [xstmp[t], t]
                elif (tmp - t) > 0:
                    below[tmp - t] = [xstmp[t], t]

            if len(above.keys()) == 0:
                # there are no cross section libraries at higher temperatures
                # then the tmp
                tmp = 0
                xs = xstmp[sorted(xstmp.keys())[-1]]
            elif len(below.keys()) == 0:
                # there are no cross section libraries at lower temperatures
                # then the tmp
                tmp = 0
                xs = xstmp[sorted(xstmp.keys())[0]]
            else:
                # there are xs libraries that bound the temperature
                xs_a = above[sorted(above.keys())[0]][0]
                tmp_a = above[sorted(above.keys())[0]][1]
                xs_b = below[sorted(below.keys())[0]][0]
                tmp_b = below[sorted(below.keys())[0]][1]
                f_a = (tmp_b**0.5 - tmp**0.5)/(tmp_b**0.5 - tmp_a**0.5)
                f_b = 1 - f_a

            if scat != '':
                import math
                scat_err = 2000
                scat_lib = ''
                for t in scattmp[scat].keys():
                    if math.fabs(t - tmp) < scat_err:
                        scat_lib = scattmp[scat][t]
                        scat_err = math.fabs(t - tmp)

            if tmp != 0:
                for nuclide in sorted(not_in):
                    if nuclide == '60120':
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (
                                '6000', xs_b, (f_b*self.comp[nuclide]), '6000', xs_a, (f_a*self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (
                                '6000', xs_b, 1e-30, '6000', xs_a, 1e-30)
                    elif nuclide not in metastable:
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (nuclide[
                                :-1],
                                xs_b,
                                (f_b*self.comp[nuclide]),
                                nuclide[
                                    :-1],
                                xs_a,
                                (f_a*self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (nuclide[
                                :-1],
                                xs_b,
                                1.e-30,
                                nuclide[
                                    :-1],
                                xs_a,
                                1.e-30,
                                )
                    else:
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (metastable[nuclide][
                                :-4],
                                xs_b,
                                (f_b*self.comp[nuclide]),
                                metastable[nuclide][
                                    :-4],
                                xs_a,
                                (f_a*self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (metastable[nuclide][
                                :-4],
                                xs_b,
                                1.e-30,
                                metastable[nuclide][
                                    :-4],
                                xs_a,
                                1.e-30,
                                )
                for nuclide in sorted(self.structonly_nucl) + sorted(self.actin_nucl) + sorted(self.fisprod_nucl):
                    if nuclide in self.comp.keys():
                        if nuclide == '60120':
                            if self.comp[nuclide] > 1.e-30:
                                mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (
                                    '6000', xs_b, (f_b*self.comp[nuclide]), '6000', xs_a, (f_a*self.comp[nuclide]))
                            else:
                                mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (
                                    '6000', xs_b, 1e-30, '6000', xs_a, 1e-30)
                        elif nuclide not in metastable:
                            if nuclide in self.comp.keys():
                                if self.comp[nuclide] > 1.e-30:
                                    mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (nuclide[
                                        :-1],
                                        xs_b,
                                        (f_b*self.comp[nuclide]),
                                        nuclide[
                                            :-1],
                                        xs_a,
                                        (f_a*self.comp[nuclide]))
                                else:
                                    mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (nuclide[
                                        :-1],
                                        xs_b,
                                        1.e-30,
                                        nuclide[
                                            :-1],
                                        xs_a,
                                        1.e-30)
                        else:
                            if self.comp[nuclide] > 1.e-30:
                                mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (metastable[nuclide][
                                    :-4],
                                    xs_b,
                                    (f_b*self.comp[nuclide]),
                                    metastable[nuclide][
                                        :-4],
                                    xs_a,
                                    (f_a*self.comp[nuclide]))
                            else:
                                mcnp5_mat += '\n          %s.%s %1.5E %s.%s %1.5E' % (metastable[nuclide][
                                    :-4],
                                    xs_b,
                                    1.e-30,
                                    metastable[nuclide][
                                        :-4],
                                    xs_a,
                                    1.e-30)
        if tmp == 0:
            for nuclide in sorted(not_in):
                if nuclide == '60120':
                    if self.comp[nuclide] > 1.e-30:
                        mcnp5_mat += '\n          %s.%s %1.5E' % ('6000',
                                                                  xs,
                                                                  (self.comp[nuclide]))
                    else:
                        mcnp5_mat += '\n          %s.%s %1.5E' % ('6000',
                                                                  xs,
                                                                  1e-30)
                elif nuclide not in metastable:
                    if self.comp[nuclide] > 1.e-30:
                        mcnp5_mat += '\n          %s.%s %1.5E' % (
                            nuclide[:-1], xs, (self.comp[nuclide]))
                    else:
                        mcnp5_mat += '\n          %s.%s %1.5E' % (
                            nuclide[: -1], xs, 1.e-30)
                else:
                    if self.comp[nuclide] > 1.e-30:
                        mcnp5_mat += '\n          %s.%s %1.5E' % (
                            metastable[nuclide][: -4], xs,
                            (self.comp[nuclide]))
                    else:
                        mcnp5_mat += '\n          %s.%s %1.5E' % (
                            metastable
                            [nuclide]
                            [: -4], xs,
                            1.e-30)

            for nuclide in sorted(self.structonly_nucl) + sorted(self.actin_nucl) + sorted(self.fisprod_nucl):
                if nuclide in self.comp.keys():
                    if nuclide == '60120':
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E' % ('6000',
                                                                      xs,
                                                                      (self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E' % (
                                '6000', xs, 1e-30)
                    elif nuclide not in metastable:
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E' % (
                                nuclide[:-1], xs, (self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E' % (
                                nuclide[: -1], xs, 1.e-30)
                    else:
                        if self.comp[nuclide] > 1.e-30:
                            mcnp5_mat += '\n          %s.%s %1.5E' % (
                                metastable[nuclide][: -4], xs,
                                (self.comp[nuclide]))
                        else:
                            mcnp5_mat += '\n          %s.%s %1.5E' % (
                                metastable
                                [nuclide]
                                [: -4], xs,
                                1.e-30)

        if scat != '':
            mcnp5_mat += '\nmt%d %s.%s' % (m, scat, scat_lib)

        self.mat_card = mcnp5_mat

        if mcf_loc != '':
            open(mcf_loc, 'w').write(mcnp5_mat)

        return mcnp5_mat

    def mcf_mcnp_nuclides(self, m, mcf_loc='', inp_loc='inp1'):

        import re

        a = depletion.Depletion()
        nct = a.nct

        nuclides = {}
        not_in = {}

        mat2 = material()
        mat2.import_mcf(inp_loc, m)

        mcnp_nuclides = re.compile(
            '\s\d{4,5}[.]\d{2,2}[cC]\s').findall(mat2.mat_card)

        mcnp5_mat = 'm%d $' % int(m)

        for mcnp_nucl in mcnp_nuclides:
            mcnp_nucl = mcnp_nucl.replace('\n', '').replace(' ', '')
            nucl = nct[mcnp_nucl]
            if nucl in self.comp.keys():
                mcnp5_mat += '\n' + ' '*10 + \
                    mcnp_nucl + ' %1.5E' % self.comp[nucl]
            else:
                mcnp5_mat += '\n' + ' '*10 + mcnp_nucl + ' %1.5E' % (1e-30)
        if mcf_loc != '':
            open(mcf_loc, 'w').write(mcnp5_mat)

        self.mat_card = mcnp5_mat
        return mcnp5_mat

    def mocup_strings(self, n=101, lib='72c'):
        self.ORIGEN_nucl_list()

        a = depletion.Depletion()
        nct = a.nct

        metastable = {}

        poop

        for nucl in nct.keys():
            if nct[nucl][-1] == '1' and (nct[nucl] not in metastable):
                if nct[nucl] == '952421':
                    # recall in mcnp metastable naming convention is reversed
                    # for am242
                    metastable['952420'] = nct[nucl][:-4] + \
                        str(int(nct[nucl][-4:-1])+400) + '.00c'
                else:
                    print nct[nucl]
                    print nct[nucl][:-4] + str(int(nct[nucl][-4:-1])+400) + '.00c'
                    metastable[nct[nucl]] = nct[nucl][
                        :-4] + str(int(nct[nucl][-4:-1])+400) + '.00c'

        print metastable
        poop

        self.single_mat = 'c    Single isotopes for depletion'
        self.tally = 'c    begin_mocup_reaction_rate_tallies\nc    time dependent reaction rates\nfc104  Reaction rates\nf4:n\n          _cells_\nfm104  (1)'

        for nucl in sorted(self.actin_nucl) + sorted(self.fisprod_nucl):
            if nucl in self.comp.keys():
                if nucl in self.actin_nucl:
                    s = '(-6)'
                else:
                    s = '    '

                if nucl not in metastable:
                    self.single_mat += '\nm%d %s.%s 1.0' % (n, nucl[:-1], lib)
                else:
                    self.single_mat += '\nm%d %s.%s 1.0' % (n,
                                                            metastable[nucl][
                                                                :-4],
                                                            lib)
                self.tally += '\n          (1.0 %d (16) (17) %s (102))' % (n, s)
                n += 1

    def scf(
            self,
            title,
            scf_loc='',
            volume='',
            temperature=0,
            density='mole/cc',
            burn='no',
            options='',
            library='03c'):

        # This method prints a serpent composition file, scf
        # options
        # 	title:       this is the title of the material
        #	scf_loc:     this is file location
        #		if this is left blank -- '' -- then the scf will not be printed
        #	volume:      this is the volume in cm
        # 	temperature: this is the temperatures in Kelvin
        #	density:     this is the units of the composition vector
        #		mole/cc
        #		at/bn-cm
        #		if this is a number, this will be the normalization factor
        #	burn:        this indicates when this is a burnable materials or not
        #	options:     this is a string to add additional options in the scf
        #	library:     this is the cross section library (ie. XXc in ZZAAAM.XXc)

        if str(density)[0] == 'M' or str(density)[0] == 'm':
            multiplier = 1.
        elif str(density)[0] == 'A' or str(density)[0] == 'a':
            multiplier = (1./.602213129)
        else:
            multiplier = density/self.mass()

        self.mat_card = '\nmat ' + title + ' '

        self.mat_card += '-%1.5E ' % (self.mass()*multiplier)

        if volume != '':
            self.mat_card += 'vol %1.5E ' % volume

        if temperature != 0:
            self.mat_card += 'tmp %1.5E ' % temperature

        if burn == 'yes' or burn == '1':
            self.mat_card += 'burn 1 '

        if options != '':
            self.mat_card += options

        for nucl in self.comp.keys():
            Z = nucl[:-4]

            if len(Z) < 2:
                Z = ' ' + Z

            if Z in [' 6', '23']:
                A = '000'
            else:
                A = nucl[-4:-1]

            M = nucl[-1]

            if M == '1':
                # this isotope is metastable
                # update the A
                A = '3' + A[1:]

            self.mat_card += '\n%s%s.%s %1.5E' % (Z,
                                                  A,
                                                  library,
                                                  self.comp[nucl])

        self.mat_card += '\n'
        if scf_loc != '':
            open(scf_loc, 'w').write(self.mat_card)

    def kcf(
            self,
            mat_num,
            title='',
            kcf_loc='',
            temperature=0,
            density='mole/cc'):

        # This method prints a serpent composition file, scf
        # options
        #   mat_num:     this is the material number
        #   title:       this is the title of the material
        #   scf_loc:     this is file location
        #       if this is left blank -- '' -- then the scf will not be printed
        #   temperature: this is the temperatures in Kelvin
        #   density:     this is the units of the composition vector
        #       mole/cc
        #       at/bn-cm
        #       if this is a number, this will be the normalization factor

        if str(density)[0] == 'M' or str(density)[0] == 'm':
            multiplier = (1/.602213129)
        elif str(density)[0] == 'A' or str(density)[0] == 'a':
            multiplier = (.602213129)
        else:
            multiplier = density/self.mass()*(.602213129)

        # convert mat_num into an int
        mat_num = int(mat_num)

        if title == '':
            title = 'material %d' % mat_num

        self.mat_card = '\n\'' + title + ' '

        ZtoSymbol = {
            1: 'h', 2: 'he', 3: 'li', 4: 'be', 5: 'b', 6: 'c', 7: 'n', 8: 'o', 9            : 'f', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p', 16            : 's', 17: 'cl', 18: 'ar', 19: 'k', 20: 'ca', 21: 'sc', 22: 'ti', 23            : 'v', 24: 'cr', 25: 'mn', 26: 'fe', 27: 'co', 28: 'ni', 29: 'cu',
            30: 'zn', 31: 'ga', 32: 'ge', 33: 'as', 34: 'se', 35: 'br', 36:
            'kr', 37: 'rb', 38: 'sr', 39: 'y', 40: 'zr', 41: 'nb', 42: 'mo', 43
            : 'tc', 44: 'ru', 45: 'rh', 46: 'pd', 47: 'ag', 48: 'cd', 49: 'in',
            50: 'sn', 51: 'sb', 52: 'te ', 53: 'i', 54: 'xe', 55: 'cs', 56:
             'ba', 57: 'la', 58: 'ce', 59: 'pr', 60: 'nd', 61: 'pm', 62: 'sm',
            63: 'eu', 64: 'gd', 65: 'tb', 66: 'dy', 67: 'ho', 68: 'er', 69:
             'tm', 70: 'yb', 71: 'lu', 72: 'hf', 73: 'ta', 74: 'w', 75: 're', 76            : 'os', 77: 'ir', 78: 'pt', 79: 'au', 80: 'hg', 81: 'tl', 82: 'pb',
            83: 'bi', 84: 'po', 85: 'at', 86: 'rn', 87: 'fr', 88: 'ra', 89:
             'ac', 90: 'th', 91: 'pa', 92: 'u', 93: 'np', 94: 'pu', 95: 'am', 96            : 'cm', 97: 'bk', 98: 'cf', 99: 'es', 100: 'fm'}

        for nucl in self.comp.keys():
            Z = int(nucl[:-4])
            A = int(nucl[-4:-1])
            M = nucl[-1]

            if Z in [6, 8]:
                isotope = ZtoSymbol[Z]
            else:
                isotope = ZtoSymbol[Z] + '-' + str(A)

            self.mat_card += '\n  %s %d 0 %1.5E %1.5E end' % (isotope,
                                                              mat_num,
                                                              (self.comp[nucl]*multiplier),
                                                              temperature)

        self.mat_card += '\n'
        if kcf_loc != '':
            open(kcf_loc, 'w').write(self.mat_card)

    def addnux(self):

        self.comp = {}
        self.comp['902320'] = 0
        self.comp['902330'] = 0
        self.comp['912330'] = 0
        self.comp['922330'] = 0
        self.comp['922340'] = 0
        self.comp['922350'] = 0
        self.comp['922360'] = 0
        self.comp['922370'] = 0
        self.comp['922380'] = 0
        self.comp['922390'] = 0
        self.comp['932360'] = 0
        self.comp['932370'] = 0
        self.comp['932380'] = 0
        self.comp['932390'] = 0
        self.comp['942360'] = 0
        self.comp['942370'] = 0
        self.comp['942380'] = 0
        self.comp['942390'] = 0
        self.comp['942400'] = 0
        self.comp['942410'] = 0
        self.comp['942420'] = 0
        self.comp['942430'] = 0
        self.comp['942440'] = 0
        self.comp['952410'] = 0
        self.comp['952421'] = 0
        self.comp['952420'] = 0
        self.comp['952430'] = 0
        self.comp['952440'] = 0
        self.comp['962420'] = 0
        self.comp['962430'] = 0
        self.comp['962440'] = 0
        self.comp['962450'] = 0
        self.comp['962460'] = 0
        self.comp['962470'] = 0
        self.comp['962480'] = 0
        self.comp['972490'] = 0
        self.comp['982490'] = 0
        self.comp['982500'] = 0
        self.comp['350810'] = 0
        self.comp['360830'] = 0
        self.comp['360840'] = 0
        self.comp['360860'] = 0
        self.comp['370850'] = 0
        self.comp['370870'] = 0
        self.comp['380880'] = 0
        self.comp['380890'] = 0
        self.comp['380900'] = 0
        self.comp['390890'] = 0
        self.comp['390910'] = 0
        self.comp['400910'] = 0
        self.comp['400920'] = 0
        self.comp['400930'] = 0
        self.comp['400940'] = 0
        self.comp['400950'] = 0
        self.comp['400960'] = 0
        self.comp['420950'] = 0
        self.comp['420970'] = 0
        self.comp['420980'] = 0
        self.comp['420990'] = 0
        self.comp['421000'] = 0
        self.comp['430990'] = 0
        self.comp['441010'] = 0
        self.comp['441020'] = 0
        self.comp['441030'] = 0
        self.comp['441040'] = 0
        self.comp['441050'] = 0
        self.comp['441060'] = 0
        self.comp['451030'] = 0
        self.comp['451050'] = 0
        self.comp['461040'] = 0
        self.comp['461050'] = 0
        self.comp['461060'] = 0
        self.comp['461070'] = 0
        self.comp['461080'] = 0
        self.comp['461100'] = 0
        self.comp['471090'] = 0
        self.comp['481100'] = 0
        self.comp['481110'] = 0
        self.comp['481130'] = 0
        self.comp['481140'] = 0
        self.comp['491150'] = 0
        # self.comp['521200'] = 0
        self.comp['531270'] = 0
        self.comp['531290'] = 0
        self.comp['541310'] = 0
        self.comp['541320'] = 0
        self.comp['541340'] = 0
        self.comp['541350'] = 0
        self.comp['541360'] = 0
        self.comp['551330'] = 0
        self.comp['551340'] = 0
        self.comp['551350'] = 0
        self.comp['551370'] = 0
        self.comp['561380'] = 0
        self.comp['561400'] = 0
        self.comp['571390'] = 0
        self.comp['581410'] = 0
        self.comp['581420'] = 0
        self.comp['581430'] = 0
        self.comp['591410'] = 0
        self.comp['591430'] = 0
        self.comp['601430'] = 0
        self.comp['601440'] = 0
        self.comp['601450'] = 0
        self.comp['601460'] = 0
        self.comp['601470'] = 0
        self.comp['601480'] = 0
        self.comp['601500'] = 0
        self.comp['611470'] = 0
        self.comp['611480'] = 0
        self.comp['611481'] = 0
        self.comp['611490'] = 0
        self.comp['621470'] = 0
        self.comp['621490'] = 0
        self.comp['621500'] = 0
        self.comp['621510'] = 0
        self.comp['621520'] = 0
        self.comp['621530'] = 0
        self.comp['621540'] = 0
        self.comp['631530'] = 0
        self.comp['631540'] = 0
        self.comp['631550'] = 0
        self.comp['631560'] = 0
        self.comp['641550'] = 0
        self.comp['641560'] = 0
        self.comp['641570'] = 0
        self.comp['641580'] = 0

    def addnux_1(self):

        # This object sets the composition vector to a blank vector with all
        # isotopes in TRITON addnux=1

        self.comp = {}

        self.comp['922340'] = 0.
        self.comp['922350'] = 0.
        self.comp['922360'] = 0.
        self.comp['922380'] = 0.
        self.comp['932370'] = 0.
        self.comp['942380'] = 0.
        self.comp['942390'] = 0.
        self.comp['942400'] = 0.
        self.comp['942410'] = 0.
        self.comp['942420'] = 0.
        self.comp['952410'] = 0.
        self.comp['952420'] = 0.
        self.comp['952430'] = 0.
        self.comp['962420'] = 0.
        self.comp['962430'] = 0.
        self.comp['541350'] = 0.
        self.comp['621490'] = 0.

    def addnux_2(self):

        # This object sets the composition vector to a blank vector with all
        # isotopes in TRITON addnux=2

        self.addnux_1()

        self.comp['360830'] = 0.
        self.comp['410930'] = 0.
        self.comp['400940'] = 0.
        self.comp['420950'] = 0.
        self.comp['430990'] = 0.
        self.comp['451030'] = 0.
        self.comp['451050'] = 0.
        self.comp['441060'] = 0.
        self.comp['471090'] = 0.
        self.comp['501260'] = 0.

        self.comp['531350'] = 0.
        self.comp['541310'] = 0.
        self.comp['551330'] = 0.
        self.comp['551340'] = 0.
        self.comp['551350'] = 0.
        self.comp['551370'] = 0.
        self.comp['591430'] = 0.
        self.comp['581440'] = 0.
        self.comp['601430'] = 0.

        self.comp['601450'] = 0.
        self.comp['601460'] = 0.
        self.comp['601470'] = 0.
        self.comp['611470'] = 0.
        self.comp['611480'] = 0.
        self.comp['611490'] = 0.
        self.comp['601480'] = 0.
        self.comp['621470'] = 0.
        self.comp['621500'] = 0.

        self.comp['621510'] = 0.
        self.comp['621520'] = 0.
        self.comp['631510'] = 0.
        self.comp['631530'] = 0.
        self.comp['631540'] = 0.
        self.comp['631550'] = 0.
        self.comp['641520'] = 0.
        self.comp['641540'] = 0.
        self.comp['641550'] = 0.
        self.comp['641560'] = 0.

        self.comp['641570'] = 0.
        self.comp['641580'] = 0.
        self.comp['641600'] = 0.
        self.comp['962440'] = 0.

    def addnux_3(self):

        # This object sets the composition vector to a blank vector with all
        # isotopes in TRITON addnux=3

        self.addnux_2()

        self.comp['320720'] = 0.
        self.comp['320730'] = 0.
        self.comp['320740'] = 0.
        self.comp['320760'] = 0.
        self.comp['330750'] = 0.
        self.comp['350790'] = 0.
        self.comp['340760'] = 0.
        self.comp['340770'] = 0.
        self.comp['340780'] = 0.
        self.comp['340800'] = 0.

        self.comp['340820'] = 0.
        self.comp['350810'] = 0.
        self.comp['360800'] = 0.
        self.comp['360820'] = 0.
        self.comp['360840'] = 0.
        self.comp['360850'] = 0.
        self.comp['360860'] = 0.
        self.comp['370850'] = 0.
        self.comp['370860'] = 0.
        self.comp['370870'] = 0.

        self.comp['380890'] = 0.
        self.comp['380900'] = 0.
        self.comp['390890'] = 0.
        self.comp['390900'] = 0.
        self.comp['390910'] = 0.
        self.comp['400900'] = 0.
        self.comp['400910'] = 0.
        self.comp['400920'] = 0.
        self.comp['400930'] = 0.
        self.comp['400950'] = 0.

        self.comp['400960'] = 0.
        #self.comp['420920'] = 0.
        #self.comp['420940'] = 0.
        self.comp['420960'] = 0.
        self.comp['420970'] = 0.
        self.comp['420980'] = 0.
        self.comp['420990'] = 0.
        self.comp['421000'] = 0.
        self.comp['410940'] = 0.
        self.comp['410950'] = 0.

        #self.comp['440960'] = 0.
        #self.comp['440980'] = 0.
        self.comp['440990'] = 0.
        self.comp['441000'] = 0.
        self.comp['441010'] = 0.
        self.comp['441020'] = 0.
        self.comp['441030'] = 0.
        self.comp['441040'] = 0.
        self.comp['441050'] = 0.
        #self.comp['461020'] = 0.

        self.comp['461040'] = 0.
        self.comp['461050'] = 0.
        self.comp['461060'] = 0.
        self.comp['461070'] = 0.
        self.comp['461080'] = 0.
        self.comp['461100'] = 0.
        self.comp['471070'] = 0.
        self.comp['471110'] = 0.
        #self.comp['481060'] = 0.
        self.comp['481080'] = 0.

        self.comp['481100'] = 0.
        self.comp['481110'] = 0.
        self.comp['481120'] = 0.
        self.comp['481130'] = 0.
        self.comp['481140'] = 0.
        self.comp['481151'] = 0.
        self.comp['481160'] = 0.
        self.comp['491130'] = 0.
        self.comp['491150'] = 0.
        #self.comp['501120'] = 0.

        self.comp['501140'] = 0.
        self.comp['501150'] = 0.
        self.comp['501160'] = 0.
        self.comp['501170'] = 0.
        self.comp['501180'] = 0.
        self.comp['501190'] = 0.
        self.comp['501200'] = 0.
        self.comp['501220'] = 0.
        self.comp['501230'] = 0.
        self.comp['501240'] = 0.

        self.comp['501250'] = 0.
        self.comp['511210'] = 0.
        self.comp['511230'] = 0.
        self.comp['511240'] = 0.
        self.comp['511250'] = 0.
        self.comp['511260'] = 0.
        #self.comp['521200'] = 0.
        self.comp['521220'] = 0.
        self.comp['521230'] = 0.
        #self.comp['521240'] = 0.

        self.comp['521250'] = 0.
        #self.comp['521260'] = 0.
        self.comp['521271'] = 0.
        self.comp['521280'] = 0.
        self.comp['521291'] = 0.
        self.comp['521300'] = 0.
        self.comp['521320'] = 0.
        self.comp['531270'] = 0.
        self.comp['531290'] = 0.
        self.comp['531300'] = 0.

        self.comp['531310'] = 0.
        #self.comp['541240'] = 0.
        #self.comp['541260'] = 0.
        self.comp['541280'] = 0.
        self.comp['541290'] = 0.
        self.comp['541300'] = 0.
        self.comp['541320'] = 0.
        self.comp['541330'] = 0.
        self.comp['541340'] = 0.
        self.comp['541360'] = 0.

        self.comp['561340'] = 0.
        self.comp['561350'] = 0.
        self.comp['561360'] = 0.
        self.comp['561370'] = 0.
        self.comp['561380'] = 0.
        self.comp['561400'] = 0.
        self.comp['551360'] = 0.
        self.comp['571390'] = 0.
        self.comp['591410'] = 0.
        self.comp['591420'] = 0.

        self.comp['571400'] = 0.
        self.comp['601420'] = 0.
        self.comp['601440'] = 0.
        self.comp['601500'] = 0.
        self.comp['581400'] = 0.
        self.comp['581410'] = 0.
        self.comp['581420'] = 0.
        self.comp['581430'] = 0.
        self.comp['611510'] = 0.
        #self.comp['621440'] = 0.

        self.comp['621480'] = 0.
        self.comp['621530'] = 0.
        self.comp['621540'] = 0.
        self.comp['631520'] = 0.
        self.comp['631560'] = 0.
        self.comp['631570'] = 0.
        self.comp['651590'] = 0.
        self.comp['651600'] = 0.
        self.comp['661600'] = 0.
        self.comp['661610'] = 0.
        self.comp['661620'] = 0.
        self.comp['661630'] = 0.

        self.comp['661640'] = 0.
        self.comp['671650'] = 0.
        self.comp['681660'] = 0.
        self.comp['681670'] = 0.
        #self.comp['711750'] = 0.
        #self.comp['711760'] = 0.
        #self.comp['731810'] = 0.
        #self.comp['741830'] = 0.
        #self.comp['741840'] = 0.

        #self.comp['741860'] = 0.
        #self.comp['751850'] = 0.
        #self.comp['751870'] = 0.
        #self.comp['791970'] = 0.
        self.comp['912310'] = 0.
        self.comp['912330'] = 0.
        self.comp['902300'] = 0.
        self.comp['902320'] = 0.
        self.comp['922320'] = 0.
        self.comp['922330'] = 0.
        self.comp['952421'] = 0.

    def ORIGEN_nucl_list(self, lib=None):

        #import xslibrary

        import re

        a = depletion.Depletion()
        nct = a.nct

        if lib is None:
            decay_lib = a.dir + 'decay.lib'
        else:
            decay_lib = a.dir + lib + '.lib'
        print lib

        xslib = open(decay_lib).read()

        # split library into structural, actinide and fission product libraries

        structlib = xslib.split('\n  -1')[0]
        actinlib = xslib.split('\n  -1')[1]
        fisprodlib = xslib.split('\n  -1')[2]

        # find the library numbers of each library

        self.structlib_num = structlib.split()[0]
        self.actinlib_num = actinlib.split()[0]
        self.fisprodlib_num = fisprodlib.split()[0]

        # determine the structural nuclides

        self.actin_nucl = []
        _compile_ = '\\s%s\s+\d+' % self.actinlib_num
        for token in re.compile(_compile_).findall(actinlib):
            self.actin_nucl.append(token.split()[1])

        self.fisprod_nucl = []
        _compile_ = '\\s%s\s{2,3}\d+' % self.fisprodlib_num
        for token in re.compile(_compile_).findall(fisprodlib):
            self.fisprod_nucl.append(token.split()[1])

        self.struct_nucl = []
        self.structonly_nucl = []
        _compile_ = '\\s%s\s{2,3}\d{5,6}' % self.structlib_num
        for token in re.compile(_compile_).findall(structlib):
            self.struct_nucl.append(token.split()[1])
            if token.split()[1] not in (self.actin_nucl + self.fisprod_nucl):
                self.structonly_nucl.append(token.split()[1])

    def mocup_strings(self, n=101, lib='72c'):
        self.ORIGEN_nucl_list()

        a = depletion.Depletion()
        nct = a.nct

        metastable = {}

        for nucl in nct.keys():
            if nct[nucl][-1] == '1' and (nct[nucl] not in metastable):
                if nct[nucl] == '952421':
                    # recall in mcnp metastable naming convention is reversed
                    # for am242
                    metastable['952420'] = nct[nucl][:-4] + \
                        str(int(nct[nucl][-4:-1])+400) + '.00c'
                else:
                    metastable[nct[nucl]] = nct[nucl][
                        :-4] + str(int(nct[nucl][-4:-1])+400) + '.00c'

        self.single_mat = 'c    Single isotopes for depletion'
        self.tally = 'c    begin_mocup_reaction_rate_tallies\nc    time dependent reaction rates\nfc4  Reaction rates\nf4:n\n          _cells_\nfm4  (1)'

        for nucl in sorted(self.actin_nucl) + sorted(self.fisprod_nucl):
            if nucl in self.comp.keys():
                if nucl in self.actin_nucl:
                    s = '(-6)'
                else:
                    s = '    '

                if nucl not in metastable:
                    self.single_mat += '\nm%d %s.%s 1.0' % (n, nucl[:-1], lib)
                else:
                    self.single_mat += '\nm%d %s.%s 1.0' % (n,
                                                            metastable[nucl][
                                                                :-4],
                                                            lib)
                self.tally += '\n          (1.0 %d (16) (17) %s (102))' % (n, s)
                n += 1

    def decay_heat(self, lib=''):

        # This method uses the materials composition list to generate an origen
        # input deck for determining the decay heat evolution, execute ORIGEN
        # and then return a vector of decay heat v. time.

        import os
        import re

        self.ORIGEN_nucl_list(lib=self.lib)

        # determine if there are any activation products
        structonly = 'no'
        actin = 'no'
        fissprod = 'no'
        for nucl in self.comp.keys():
            if nucl in self.structonly_nucl:
                structonly = 'yes'
            elif nucl in self.actin_nucl:
                actin = 'yes'
            elif nucl in self.fisprod_nucl:
                fissprod = 'yes'
        skele = '  -1                                                                            \n  -1                                                                            \n  -1                                                                            \n  TIT   Decay                                  \n  BAS   DECAY HEAT CALCULATION                                                       \n  RDA   A. Cisneros Generated by mocup.py material.decay_heat function\n  LIP   0 0 0                                                                   \n  LIB   0 1 2 3    _STRUCT_LIB_ _ACTIN_LIB_  _FISSPROD_LIB_ 9 50 0 4 0                                     \n  PHO   0 102 103 10                                                            \n  _OPTL_ \n  _OPTA_ \n  _OPTF_ \n  CUT   3 1.0E-24  28 1.0E-75          -1                                       \n  INP   1 -2  -1  -1  1  1                                                      \n  DEC      1 1 2 1 0\n  DEC      2 2 3 1 0\n  DEC      4 3 4 1 0\n  DEC      8 4 5 1 0\n  DEC      16 5 6 1 0\n  DEC      24 6 7 1 0\n  DEC      32 7 8 1 0\n  DEC      40 8 9 1 0\n  DEC      48 9 10 1 0\n  DEC      60 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  DEC      120 1 2 1 0\n  DEC      240 2 3 1 0\n  DEC      480 3 4 1 0\n  DEC      960 4 5 1 0\n  DEC      1440 5 6 1 0\n  DEC      1920 6 7 1 0\n  DEC      2400 7 8 1 0\n  DEC      2880 8 9 1 0\n  DEC      3360 9 10 1 0\n  DEC      3600 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  DEC      12240 1 2 1 0\n  DEC      20880 2 3 1 0\n  DEC      29520 3 4 1 0\n  DEC      38160 4 5 1 0\n  DEC      46800 5 6 1 0\n  DEC      55440 6 7 1 0\n  DEC      64080 7 8 1 0\n  DEC      72720 8 9 1 0\n  DEC      81360 9 10 1 0\n  DEC      86400 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  DEC      172800 1 2 1 0\n  DEC      259200 2 3 1 0\n  DEC      345600 3 4 1 0\n  DEC      432000 4 5 1 0\n  DEC      518400 5 6 1 0\n  DEC      604800 6 7 1 0\n  DEC      691200 7 8 1 0\n  DEC      777600 8 9 1 0\n  DEC      864000 9 10 1 0\n  DEC      950400 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  DEC      1036800 1 2 1 0\n  DEC      1123200 2 3 1 0\n  DEC      1209600 3 4 1 0\n  DEC      1296000 4 5 1 0\n  DEC      1382400 5 6 1 0\n  DEC      1468800 6 7 1 0\n  DEC      1555200 7 8 1 0\n  DEC      1641600 8 9 1 0\n  DEC      1728000 9 10 1 0\n  DEC      1814400 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  DEC      1900800 1 2 1 0\n  DEC      1987200 2 3 1 0\n  DEC      2073600 3 4 1 0\n  DEC      2160000 4 5 1 0\n  DEC      2246400 5 6 1 0\n  DEC      2332800 6 7 1 0\n  DEC      2419200 7 8 1 0\n  DEC      2505600 8 9 1 0\n  DEC      2592000 9 10 1 0\n  DEC      2678400 10 11 1 0\n  OUT 11 1 0 0\n  PCH 11 11 11\n  MOV 11 1 0 1\n  STP   4                                                                       \n0                                                                               \n\n'

        if structonly == 'yes':
            _STRUCT_LIB_ = '-%s' % self.structlib_num
            _OPTL_ = 'OPTL   8 8 8 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
        else:
            _STRUCT_LIB_ = '0'
            _OPTL_ = 'OPTL   8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'

        skele = skele.replace('_STRUCT_LIB_', _STRUCT_LIB_)
        skele = skele.replace('_OPTL_', _OPTL_)

        if actin == 'yes':
            _ACTIN_LIB_ = '-%s' % self.actinlib_num
            _OPTA_ = 'OPTA   8 8 8 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
        else:
            _ACTIN_LIB_ = '0'
            _OPTA_ = 'OPTA   8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'

        skele = skele.replace('_ACTIN_LIB_', _ACTIN_LIB_)
        skele = skele.replace('_OPTA_', _OPTA_)

        if fissprod == 'yes':
            _FISSPROD_LIB_ = '-%s' % self.fisprodlib_num
            _OPTF_ = 'OPTF   8 8 8 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
        else:
            _FISSPROD_LIB_ = '0'
            _OPTF_ = 'OPTF   8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'

        skele = skele.replace('_FISSPROD_LIB_', _FISSPROD_LIB_)
        skele = skele.replace('_OPTF_', _OPTF_)

        open('TAPE5.INP', 'w').write(skele)

        self.make_ocf('TAPE4.INP', lib=self.lib)

        a = depletion.Depletion()

        # TAPE3 - blanket
        open('TAPE3.INP', 'w').write('')

        # TAPE9 - Neutron Cross Sections
        TAPE9 = 'cat %sdecay.lib %s%s.lib > TAPE9.INP' % (
            a.dir, a.dir, self.lib)
        print TAPE9
        os.system(TAPE9)

        # TAPE10 - Gamma Cross Sectinos
        TAPE10 = 'ln -s %sgxuo2brm.lib TAPE10.INP' % (a.dir)
        print TAPE10
        os.system(TAPE10)

        if self.lib[0:3] in ['amo', 'emo', 'fft']:
            # origen uses a fast reactor set of cross sections, thus origen will
            # be executed with o2_fast
            origen = 'o2_fast'
        else:
            # origen uses a fast reactor set of cross sections, thus origen will
            # be executed with o2_therm
            origen = 'o2_therm'

        print origen
        os.system(origen)

        TAPE6 = open('TAPE6.OUT').read()
        TAPE12 = open('TAPE12.OUT').read()
        tokens = re.compile('\s\d+\s{13,13}THERMAL').findall(TAPE12)

        print tokens

        sub_times = []
        sub_times.append([1.0000E+00,
                          2.0000E+00,
                          4.0000E+00,
                          8.0000E+00,
                          1.6000E+01,
                          2.4000E+01,
                          3.2000E+01,
                          4.0000E+01,
                          4.8000E+01,
                          6.0000E+01,
                          ])
        sub_times.append([1.2000E+02,
                          2.4000E+02,
                          4.8000E+02,
                          9.6000E+02,
                          1.4400E+03,
                          1.9200E+03,
                          2.4000E+03,
                          2.8800E+03,
                          3.3600E+03,
                          3.6000E+03,
                          ])
        sub_times.append([1.2240E+04,
                          2.0880E+04,
                          2.9520E+04,
                          3.8160E+04,
                          4.6800E+04,
                          5.5440E+04,
                          6.4080E+04,
                          7.2720E+04,
                          8.1360E+04,
                          8.6400E+04,
                          ])
        sub_times.append([1.7280E+05,
                          2.5920E+05,
                          3.4560E+05,
                          4.3200E+05,
                          5.1840E+05,
                          6.0480E+05,
                          6.9120E+05,
                          7.7760E+05,
                          8.6400E+05,
                          9.5040E+05,
                          ])
        sub_times.append([1.0368E+06,
                          1.1232E+06,
                          1.2096E+06,
                          1.2960E+06,
                          1.3824E+06,
                          1.4688E+06,
                          1.5552E+06,
                          1.6416E+06,
                          1.7280E+06,
                          1.8144E+06,
                          ])
        sub_times.append([1.9008E+06,
                          1.9872E+06,
                          2.0736E+06,
                          2.1600E+06,
                          2.2464E+06,
                          2.3328E+06,
                          2.4192E+06,
                          2.5056E+06,
                          2.5920E+06,
                          2.6784E+06,
                          ])

        decay_heat = {}
        substep_decay_heat = {}
        i = 0

        for time_list in sub_times:
            if structonly == 'yes':
                skeleton_token = 'PAGE     '
                token1 = skeleton_token[
                    :-len(tokens[i].split()[0])] + tokens[i].split()[0]
                k = TAPE6.find(token1)
                l = TAPE6[k:].find('0TOTAL')
                heat = TAPE6[k+l:].split()[2:12]
                for j in range(len(time_list)):
                    decay_heat[time_list[j]] = float(heat[j])*1e-6
                # advance in the decay heat table token vector
                i += 1

            if actin == 'yes':
                skeleton_token = 'PAGE     '
                token1 = skeleton_token[
                    :-len(tokens[i].split()[0])] + tokens[i].split()[0]
                k = TAPE6.find(token1)
                l = TAPE6[k:].find('0TOTAL')
                heat = TAPE6[k+l:].split()[2:12]
                for j in range(len(time_list)):
                    if time_list[j] in decay_heat.keys():
                        decay_heat[time_list[j]] += float(heat[j])*1e-6
                    else:
                        decay_heat[time_list[j]] = float(heat[j])*1e-6
                # advance in the decay heat table token vector
                i += 1

            if fissprod == 'yes':
                skeleton_token = 'PAGE     '
                token1 = skeleton_token[
                    :-len(tokens[i].split()[0])] + tokens[i].split()[0]
                k = TAPE6.find(token1)
                l = TAPE6[k:].find('0TOTAL')
                heat = TAPE6[k+l:].split()[2:12]
                for j in range(len(time_list)):
                    if time_list[j] in decay_heat.keys():
                        decay_heat[time_list[j]] += float(heat[j])*1e-6
                    else:
                        decay_heat[time_list[j]] = float(heat[j])*1e-6
                # advance in the decay heat table token vector
                i += 1

        time_out = 'Time (s),'
        decay_heat_out = '\nDecay Heat (MW),'
        for i in range(len(decay_heat)):
            time_out += '%1.5E,' % sorted(decay_heat.keys())[i]
            decay_heat_out += '%1.5E,' % decay_heat[sorted
                                                    (decay_heat.keys())[i]]
        print time_out + decay_heat_out
        return decay_heat
