#!/usr/bin/env python


class mpo:
    # This class holds all the information mocup can receive from an mcnp output file
    # cell numbers of depletion zones
    # flux tallies
    # libraries of cross sections
    #

    def __init__(self):
        self.cell = []
        self.ORIGEN_power = {}
        self.ORIGEN_fissionpower = {}
        self.ORIGEN_fluxes = {}
        self.temperature = {}

    def populate(self, timestep=None, print_mpo=None, ):

        import math
        from os import listdir
        import re

        self.transport_module = 'mcnp'

        self.timestep = timestep or ''

        a = Depletion()
        nct = a.nct
        dacay = a.decay

        # import OUTP file
        outp_loc = 'outp%d' % timestep
        OUTP = open(outp_loc).read()

        # search for cells

        OUTP1 = OUTP
        self.mcnp_nuclides = {}
        if 'c    time dependent reaction rates' in OUTP1:
            k = OUTP1.find('c    time dependent reaction rates')
            self.RR_tally_ID = re.compile('[fF]\d*4[:][nN]').search(
                OUTP1[
                    k:]).group().replace(
                'f',
                '').replace(
                'F',
                '').replace(
                ':n',
                '').replace(
                ':N',
                '')
            token = re.compile('[fF][mM]\d+').search(OUTP[k:]).group()
            l = OUTP1[k:].find(token) + len(token)

            m = re.compile('\s[mMfFkKsSpP]').search(OUTP1[k+l:]).start()

            tally = OUTP1[k+l:k+l+m]

            multipliers = re.compile(
                '\s\(\s*1[.]?\d*\s*\d+\s+\(').findall(tally)

            for multiplier in multipliers:
                self.mcnp_nuclides[multiplier.split()[-2]] = 0
        else:
            print 'There is no MOCUP tally in this input deck!'
            exit()

        z = OUTP.find('1problem')
        OUTP1 = OUTP[z:]

        # find depletion cells
        token = '1tally' + ' '*(4 - len(self.RR_tally_ID)) + self.RR_tally_ID

        # search for cells
        # slice out the cells defined in the reaction rates tally
        self.cell = []
        k = OUTP1.find(token)
        l = OUTP1[k:].find('cell:')
        m = OUTP1[k+l:].find('cell ')

        cells = OUTP1[k+l:k+l+m]
        while 'cell:' in cells:
            k = cells.find('cell:')
            l = cells[k:].find('\n')
            line = cells[k:k+l].split()
            self.cell = self.cell + line[1:]
            cells = cells[k+l:]

        # populate the volume table
        self.volume = {}
        self.density = {}
        self.molar_density = {}
        self.mat_ID = {}

        OUTP1 = OUTP
        k = OUTP1.find('1cells ')
        l = OUTP1[k:].find('importance')
        m = OUTP1[k+l:].find('total')
        table = OUTP[k+l+len('importance'):k+l+m]
        k = table.find('\n')
        for line in table.split('\n'):
            if len(line.split()) > 3 and line[7:12].replace(' ', '') in self.cell:
                self.volume[
                    line[
                        7:12].replace(
                        ' ',
                        '')] = float(
                    line[
                        43:54].replace(
                        ' ',
                        ''))
                self.density[
                    line[
                        7:12].replace(
                        ' ',
                        '')] = float(
                    line[
                        31:42].replace(
                        ' ',
                        ''))
                self.molar_density[
                    line[
                        7:12].replace(
                        ' ',
                        '')] = float(
                    line[
                        19:30].replace(
                        ' ',
                        ''))
                self.mat_ID[
                    line[
                        7:12].replace(
                        ' ',
                        '')] = line[
                    12:17].replace(
                    ' ',
                    '')

        # populate material correspondence dictionary
        OUTP1 = OUTP
        for mat in self.mcnp_nuclides.keys():
            mat_ID = 'm%s ' % mat
            k = OUTP1.find(mat_ID)
            self.mcnp_nuclides[mat] = OUTP1[k:].split()[1]

        # determine the nuclide concentration vector for each depletion material

        self.concentration = {}

        moi_files = listdir('moi_files')

        # split table 40 into individual material cards

        k = OUTP1.find('print table 40')
        k += OUTP1[k:].find('component nuclide, atom fraction') + \
            len('component nuclide, atom fraction')
        l = OUTP1[k:].find('1material')

        table40 = OUTP1[k:k+l]

        # split table 40 into individual material cards

        # find the material tokens
        tokens = re.compile('\s+\d+\s+').findall(table40)
        material_cards = {}
        for i in range(len(tokens)-1):
            k = table40.find(tokens[i])+len(tokens[i])
            l = table40.find(tokens[i+1])

            material_cards[tokens[i].split()[0]] = table40[k:l]

        material_cards[tokens[i+1].split()[0]] = table40[l+len(tokens[i+1]):]

        for cell in self.cell:

            self.concentration[cell] = {}

            # populate the concentration vector based on mcnp output file
            OUTP1 = OUTP

            # grep material vector from material card in table 40

            pairs = re.compile(
                '\d+,\s\d+[.]\d+[eE+-]+\d+').findall(material_cards[self.mat_ID[cell]])

            for pair in pairs:

                nuclide = pair.split(',')[0]

                # generate nuclide ID

                if nuclide == '6000':
                    nuclide = '6012'
                Z = int(nuclide[:-3])

                A = int(nuclide[-3:])

                # check to see if this is
                if A > 300:
                    if A == '601':
                        print "please use endf7 formatted XS's!"
                        exit()

                    A += -400
                    M = 1
                    print '%d' % (M + 10*A + 10000*Z)
                else:
                    M = 0

                nucl = '%d' % (M + 10*A + 10000*Z)

                # update concentration

                if nucl in self.concentration[cell].keys():
                    self.concentration[cell][
                        nucl] += self.molar_density[cell]*float(pair.split()[1])
                else:
                    self.concentration[cell][nucl] = self.molar_density[
                        cell]*float(pair.split()[1])

            # update the concentration vector if there is a ocf file

            ocf_loc = 'moi.%s.%d.pch' % (cell, (self.timestep - 1))
            if ocf_loc in moi_files:
                # There is a ocf file for the material vector
                # populate the concentration vector based on this material

                mat = material()
                mat.import_ocf('moi_files/'+ocf_loc)

                for nuclide in self.mcnp_nuclides.values():
                    if nct[nuclide] in mat.comp.keys():
                        self.concentration[cell][
                            nct[nuclide]] = mat.comp[
                            nct[nuclide]]/self.volume[cell]*.602214078

        # define the temperature distribution for depletion materials
        if 'print table 72' in OUTP:
            OUTP1 = OUTP
            k = OUTP1.find('print table 72') + len('print table 72')
            l = OUTP1[k:].find('\n1')
            table72 = OUTP1[k:k+l]
            pairs = re.compile(
                '\d{1,5}\s+\d[.]\d{4,4}[eE+-]{2,2}\d{2,2}').findall(table72)
            for pair in pairs:
                if pair.split()[0] in self.cell:
                    self.temperature[
                        pair.split()[0]] = (
                        float(
                            pair.split()[1]) /
                        8.6173324e-11)

        OUTP1 = OUTP
        k = OUTP1.find('1problem summary')
        OUTP1 = OUTP1[k:]

        self.flux_tally = {}

        # isolate the flux tally

        token = '1tally'+' '*(4-len(self.RR_tally_ID))+self.RR_tally_ID
        k = OUTP1.find(token)
        l = OUTP1[k:].find('1status')
        flux_tallies = OUTP1[k:k+l]

        f4_tallies = re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+[-]*\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}').findall(flux_tallies)
        for string in f4_tallies:
            self.flux_tally[string.split()[1]] = float(string.split()[5])

        def safedivide(numerator, denominator):
            try:
                return numerator/denominator
            except ZeroDivisionError:
                return 0.

        self.NG = {}
        self.NF = {}
        self.N2N = {}
        self.N3N = {}
        self.NG_sig = {}
        self.NF_sig = {}
        self.N2N_sig = {}
        self.N3N_sig = {}

        for cell in self.cell:
            self.NG[cell] = {}
            self.NF[cell] = {}
            self.N2N[cell] = {}
            self.N3N[cell] = {}
            self.NG_sig[cell] = {}
            self.NF_sig[cell] = {}
            self.N2N_sig[cell] = {}
            self.N3N_sig[cell] = {}

        # populate (n,gamma) library mt=102

        f4_tallies = re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d+\s+102\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2} \d[.]\d{4,4}').findall(flux_tallies)
        for string in f4_tallies:
            self.NG[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = safedivide(
                float(string.split()[-2]), self.flux_tally
                [string.split()[1]])
            self.NG_sig[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = float(
                string.split()[-1])

        # populate (n,fission) library mt=18 or -6

        f4_tallies = re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d+\s+18\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2} \d[.]\d{4,4}').findall(flux_tallies)
        f4_tallies += re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d+\s+-6\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2} \d[.]\d{4,4}').findall(flux_tallies)

        for string in f4_tallies:
            self.NF[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = safedivide(
                float(string.split()[-2]), self.flux_tally
                [string.split()[1]])
            self.NF_sig[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = float(
                string.split()[-1])

        # populate (n,2n) library mt=16

        f4_tallies = re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d+\s+16\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2} \d[.]\d{4,4}').findall(flux_tallies)
        for string in f4_tallies:
            self.N2N[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = safedivide(
                float(string.split()[-2]), self.flux_tally
                [string.split()[1]])
            self.N2N_sig[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = float(
                string.split()[-1])

        # populate (n,3n) library mt=17

        f4_tallies = re.compile(
            r'cell\s+\d+\s+multiplier bin:\s{3,3}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d+\s+17\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2} \d[.]\d{4,4}').findall(flux_tallies)
        for string in f4_tallies:
            self.N3N[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = safedivide(
                float(string.split()[-2]), self.flux_tally
                [string.split()[1]])
            self.N3N_sig[
                string.split()[1]][
                nct[self.mcnp_nuclides[string.split()[5]]]] = float(
                string.split()[-1])

        # determine power distribution in each cell
        self.power = {}
        self.power_density = {}
        self.decay_heat = {}
        sum_power = 0

        for cell in self.cell:
            power = 0
            power_density = 0
            power_decay = 0
            for nuclide in self.concentration[cell].keys():
                if nuclide in self.NF[cell].keys():
                    Z = float(nuclide[-6:-4])
                    A = float(nuclide[-4:-1])
                    Q = 1.29927E-3*(Z**2)*(A**.5)+33.12
                    power += Q*self.NF[cell][nuclide]*self.flux_tally[cell] * \
                        self.concentration[cell][nuclide]*self.volume[cell]
                    power_density += Q * \
                        self.NF[cell][nuclide]*self.flux_tally[cell]*self.concentration[cell][nuclide]
                # if nuclide in self.decy.keys():
                #	power_decay += self.concentration[cell][nuclide]*self.volume[cell]*1.e24*self.decay[nuclide]

            # add power in cell to power distribution
            self.power[cell] = power
            #self.decay_heat[cell] = power_decay
            self.power_density[cell] = power_density
            sum_power += power

        self.global_power_density = sum_power

        # normalize the power distribution
        if sum_power != 0.:
            for cell in self.cell:
                self.power[cell] = self.power[cell]/sum_power
        else:
            for cell in self.cell:
                self.power[cell] = 1./len(self.cell)

        if print_mpo == 'yes':
            self.write_mpo()

    def populate_serpent(self, timestep=None, print_mpo=None, ):

        # This method collects all the information required to perform the point
        # depletion calculation with origen from SERPENT output files generated
        # in externally coupled model

        import math
        import re

        self.transport_module = 'serpent'

        self.timestep = timestep

        outp_loc = 'inp%d_res.m' % self.timestep
        outp = open(outp_loc).read().split()

        # find neutron fission source
        k = outp.index('TOT_GENRATE')
        source = float(outp[k+6])

        k = outp.index('TOT_POWER')
        power = float(outp[k+6])

        burn_loc = 'inp%d.burn' % self.timestep
        burn_card = open(burn_loc).read()
        burn_cards = {}

        # cells
        # divide the burn card into segments indexed by their cell
        self.cell = []
        self.mat_ID = {}
        self.burn_cards = {}
        self.volume = {}
        self.flux_tally = {}
        self.power = {}
        self.density = {}
        self.concentration = {}
        self.power_density = {}
        self.mcnp_nuclides = {}

        self.NG = {}
        self.NF = {}
        self.N2N = {}
        self.N3N = {}

        self.NF_sig = {}
        self.NG_sig = {}
        self.N2N_sig = {}
        self.N3N_sig = {}

        # idenfify the depletion cells and split the burn_card into depletion
        # specific burn cards

        # identify the tokens for Material specific burn_cards
        tokens = re.compile(
            '% Material +\d+[/]\d+ [\(].+[\)]').findall(burn_card)

        # this definition is inserted for the case that there is only one
        # depletion cell
        i = -1

        for i in range(len(tokens)-1):
            self.cell.append(
                tokens[i].split()[-1].replace("('", '').replace("')", ''))
            self.mat_ID[self.cell[-1]] = self.cell[-1]
            k = burn_card.find(tokens[i])
            l = burn_card.find(tokens[i+1])
            self.burn_cards[self.cell[-1]] = burn_card[k:l]

        # the last burn card has different index definitions
        self.cell.append(
            tokens[i+1].split()[-1].replace("('", '').replace("')", ''))
        self.mat_ID[self.cell[-1]] = self.cell[-1]
        k = burn_card.find(tokens[i+1])
        self.burn_cards[self.cell[-1]] = burn_card[k:]

        ng_match = re.compile('102 \d \d[.]\d{5,5}E[-+]\d{2,2} \d[.]\d{5,5}')
        n2n_match = re.compile('16 \d \d[.]\d{5,5}E[-+]\d{2,2} \d[.]\d{5,5}')
        n3n_match = re.compile('17 \d \d[.]\d{5,5}E[-+]\d{2,2} \d[.]\d{5,5}')
        nf_match = re.compile('18 \d \d[.]\d{5,5}E[-+]\d{2,2} \d[.]\d{5,5}')

        for cell in self.cell:

            # nuclides concentration density volume

            burn_card = self.burn_cards[cell]

            k = burn_card.split().index('VOL')

            self.volume[cell] = float(burn_card.split()[k+1])

            k = burn_card.split().index('ADENS')

            self.density[cell] = float(burn_card.split()[k+1])

            k = burn_card.split().index('FLUX')

            self.flux_tally[cell] = float(
                burn_card.split()[k+1])/(source*self.volume[cell])

            k = burn_card.split().index('PDE')

            self.power[cell] = float(
                burn_card.split()[
                    k+1])*self.volume[cell]/power

            # 1group cross sections
            # N2N
            # N3N
            # NG
            # NF

            # isolate cross sections

            k = burn_card.find('\nXS ')
            l = burn_card.find(burn_card[k:].split()[2])
            XS = burn_card[l:]

            self.NG[cell] = {}
            self.NF[cell] = {}
            self.N2N[cell] = {}
            self.N3N[cell] = {}
            self.NG_sig[cell] = {}
            self.NF_sig[cell] = {}
            self.N2N_sig[cell] = {}
            self.N3N_sig[cell] = {}
            self.concentration[cell] = {}
            self.mcnp_nuclides = {}
            self.power_density[cell] = 0

            token = re.compile(
                '\d{4,5}[.]\d{2,2}[a-zA-Z] {2,3}\d{4,5} \d {2,4}\d{1,3}[.]\d{5,5} {3,3}\d[.]\d{5,5}E[+-]\d{2,2}').search(self.burn_cards[cell]).group()
            k = self.burn_cards[cell].find(token)
            l = self.burn_cards[cell][k:].find('\n\n')
            nucl_cards = self.burn_cards[cell][k:k+l].rsplit('\n')
            for nucl_card in nucl_cards:
                list = nucl_card.split()
                nucl = list[1]+list[2]

                if nucl == '60000':
                    nucl = '60120'

                self.concentration[cell][nucl] = float(
                    list[4])*self.density[cell]

                self.NG[cell][nucl] = 0.
                self.NG_sig[cell][nucl] = 0.
                for string in ng_match.findall(nucl_card):
                    self.NG[cell][nucl] += float(string.split()[2])
                    self.NG_sig[cell][nucl] = float(string.split()[3])

                for string in nf_match.findall(nucl_card):
                    self.NF[cell][nucl] = float(string.split()[2])
                    self.NF_sig[cell][nucl] = float(string.split()[3])
                    Z = float(nucl[-6:-4])
                    A = float(nucl[-4:-1])
                    Q = 1.29927E-3*(Z**2)*(A**.5)+33.12
                    self.power_density[cell] += Q*self.NF[cell][nucl] * \
                        self.flux_tally[cell]*self.concentration[cell][nucl]

                for string in n2n_match.findall(nucl_card):
                    self.N2N[cell][nucl] = float(string.split()[2])
                    self.N2N_sig[cell][nucl] = float(string.split()[3])

                for string in n3n_match.findall(nucl_card):
                    self.N3N[cell][nucl] = float(string.split()[2])
                    self.N3N_sig[cell][nucl] = float(string.split()[3])

        self.global_power_density = 0
        for cell in self.cell:
            self.global_power_density += self.power_density[cell] * \
                self.volume[cell]

        for nucl_card in nucl_cards:
            self.mcnp_nuclides[nucl_card.split()[0]] = nucl_card.split()[0]

        if print_mpo == 'yes':
            self.write_mpo()

    def write_mpo(self):
        # this method returns mpo files just like the real mocup

        a = Depletion()
        nct = a.nct

        for cell in self.cell:
            output = ''

            # add timestep
            output += '.%d\n1' % self.timestep

            #  add cell, material number, flux, atomic density, volume
            output += '\n %s  %s  %1.6E  %1.6E  %1.6E' % (cell,
                                                          self.mat_ID[cell],
                                                          (self.flux_tally[cell]),
                                                          self.density[cell],
                                                          self.volume[cell])

            # write number of nuclides
            output += '\n  %d' % len(self.concentration[cell].keys())

            # write nuclide information

            for nuclide in self.mcnp_nuclides.values():

                if nct[nuclide] in self.concentration[cell].keys():
                    if nct[nuclide] in self.NF[cell].keys() or nct[nuclide] in self.NG[cell].keys() or nct[nuclide] in self.N2N[cell].keys() or nct[nuclide] in self.N3N[cell].keys():
                        # this nuclide is in the material composition for this depletion cell and has at least one cross section
                        # print nuclide information mcnp ID, origen ID and
                        # concentration
                        output += '\n  %s %s %1.6E' % (nuclide,
                                                       nct[nuclide],
                                                       self.concentration[cell][
                                                           nct[nuclide]])

                        # determine the number of cross sections
                        XSn = 0
                        XS = ''

                        if nct[nuclide] in self.N2N[cell].keys():
                            XS += '\n   SN2N %1.6E %1.4f' % (self.N2N
                                                             [cell]
                                                             [nct[nuclide]],
                                                             self.N2N_sig
                                                             [cell]
                                                             [nct[nuclide]])
                            XSn += 1
                        if nct[nuclide] in self.N3N[cell].keys():
                            XS += '\n   SN3N %1.6E %1.4f' % (self.N3N
                                                             [cell]
                                                             [nct[nuclide]],
                                                             self.N3N_sig
                                                             [cell]
                                                             [nct[nuclide]])
                            XSn += 1
                        if nct[nuclide] in self.NF[cell].keys():
                            XS += '\n   SNF %1.6E %1.4f' % (self.NF
                                                            [cell]
                                                            [nct[nuclide]],
                                                            self.NF_sig
                                                            [cell]
                                                            [nct[nuclide]])
                            XSn += 1
                        if nct[nuclide] in self.NG[cell].keys():
                            XS += '\n   SNG %1.6E %1.4f' % (self.NG
                                                            [cell]
                                                            [nct[nuclide]],
                                                            self.NG_sig
                                                            [cell]
                                                            [nct[nuclide]])
                            XSn += 1

                        output += '\n   %d 1%s' % (XSn, XS)
            # write output file
            mpo_loc = 'mpo.%s.%d' % (cell, self.timestep)
            open(mpo_loc, 'w').write(output)

    def update(self):

        import math

        a = Depletion()
        nct = a.nct

        self.power = {}
        self.power_density = {}
        sum_power = 0

        for cell in self.cell:
            # import fuel composition from OCF file
            ocf_loc = 'moi_files/moi.%s.%d.pch' % (cell, self.timestep)
            mat = material()
            mat.import_ocf(ocf_loc)
            mat = mat*(.602214078/self.volume[cell])
            for nucl in self.concentration[cell].keys():
                if nucl in mat.comp.keys():
                    self.concentration[cell][nucl] = mat.comp[nucl]

            # ensure all mcnp nuclides are in composition vector
            # calculate power in cell
            power = 0
            power_density = 0

            for nuclide in self.mcnp_nuclides.values():
                if nct[nuclide] in self.concentration[cell].keys():
                    if nct[nuclide] in self.NF[cell].keys():
                        Z = float(nct[nuclide][-6:-4])
                        A = float(nct[nuclide][-4:-1])
                        Q = 1.29927E-3*math.pow(Z, 2)*math.pow(A, .5)+33.12
                        power += Q*self.NF[cell][nct[nuclide]]*self.flux_tally[cell] * \
                            self.concentration[cell][nct[nuclide]]*self.volume[cell]
                        power_density += Q * \
                            self.NF[cell][nct[nuclide]]*self.flux_tally[cell]*self.concentration[cell][nct[nuclide]]

            self.power[cell] = power
            self.power_density[cell] = power_density
            sum_power += power

        if sum_power != 0.:
            for cell in self.cell:
                self.power[cell] = self.power[cell]/sum_power
        else:
            for cell in self.cell:
                self.power[cell] = 1./len(self.cell)

    def update_mat(self, progression, eq_mat):
        # This method accepts a progression and eq_material objects from BEAU
        # and updates the mpo via the eq_material rather than an OCF file

        import math
        a = Depletion()
        nct = a.nct

        for cell in progression.cells:
            # import fuel composition from mat object

            self.concentration[cell] = eq_mat.Eq_comp[cell].comp

            # ensure all mcnp nuclides are in composition vector
            # calcuale power in cell
            power = 0
            power_density = 0

            for nuclide in self.mcnp_nuclides.values():
                if nct[nuclide] in self.concentration[cell].keys():
                    if nct[nuclide] in self.NF[cell].keys():
                        Z = float(nct[nuclide][-6:-4])
                        A = float(nct[nuclide][-4:-1])
                        Q = 1.29927E-3*math.pow(Z, 2)*math.pow(A, .5)+33.12
                        power += Q*self.NF[cell][nct[nuclide]]*self.flux_tally[cell] * \
                            self.concentration[cell][nct[nuclide]]*self.volume[cell]
                        power_density += Q * \
                            self.NF[cell][nct[nuclide]]*self.flux_tally[cell]*self.concentration[cell][nct[nuclide]]
                else:
                    self.concentration[cell][nct[nuclide]] = 0
            self.power[cell] = power
            self.power_density[cell] = power_density

        sum_power = 0
        for power in self.power.values():
            sum_power += power

        if sum_power != 0.:
            for cell in self.cell:
                self.power[cell] = self.power[cell]/sum_power
        else:
            for cell in self.cell:
                self.power[cell] = 1./len(self.cell)

    def split(self):

        # This method slits the mpo up into a list of mpos that only contain the
        # depletion information for a single cell

        # generate several new mpos by copying original mpo
        mpos = {}
        for cell in self.cell:
            mpo1 = mpo()
            mpos[cell] = mpo1

            # replace full system information with cell specific information

            mpos[cell[:]].timestep = self.timestep
            mpos[cell[:]].transport_module = self.transport_module[:]
            mpos[cell[:]].cell = [cell[:]]
            mpos[
                cell[:]].concentration = {
                cell[:]: self.concentration[cell].copy()}
            mpos[cell[:]].density = {cell[:]: self.density[cell]}
            mpos[cell[:]].volume = {cell[:]: self.volume[cell]}

            mpos[cell[:]].N2N = {cell[:]: self.N2N[cell].copy()}
            mpos[cell[:]].N3N = {cell[:]: self.N3N[cell].copy()}
            mpos[cell[:]].NG = {cell[:]: self.NG[cell].copy()}
            mpos[cell[:]].NF = {cell[:]: self.NF[cell].copy()}

            mpos[cell[:]].N2N_sig = {cell[:]: self.N2N_sig[cell].copy()}
            mpos[cell[:]].N3N_sig = {cell[:]: self.N3N_sig[cell].copy()}
            mpos[cell[:]].NG_sig = {cell[:]: self.NG_sig[cell].copy()}
            mpos[cell[:]].NF_sig = {cell[:]: self.NF_sig[cell].copy()}

            mpos[cell[:]].power = {cell[:]: self.power[cell]}
            mpos[cell[:]].power_density = {cell[:]: self.power_density[cell]}
            mpos[cell[:]].flux_tally = {cell[:]: self.flux_tally[cell]}
            mpos[cell[:]].mat_ID = {cell[:]: self.mat_ID[cell][:]}
            mpos[cell[:]].mcnp_nuclides = self.mcnp_nuclides.copy()
            if cell in self.ORIGEN_power.keys():
                mpos[cell[:]].ORIGEN_power[cell[:]] = self.ORIGEN_power[cell]
            if cell in self.ORIGEN_fissionpower.keys():
                mpos[
                    cell[:]].ORIGEN_fissionpower[
                    cell[:]] = self.ORIGEN_fissionpower[cell]
            if cell in self.ORIGEN_fluxes.keys():
                mpos[cell[:]].ORIGEN_fluxes[cell[:]] = self.ORIGEN_fluxes[cell]
            if cell in self.temperature.keys():
                mpos[cell[:]].temperature[cell[:]] = self.temperature[cell]

            # Note mcnp_nuclides list, timestep are left alone because they
            # apply to all cells

        return mpos

    def twin(self):
        mpo1 = mpo()

        mpo1.timestep = self.timestep
        mpo1.transport_module = self.transport_module[:]
        mpo1.cell = self.cell[:]
        mpo1.concentration = self.concentration.copy()
        mpo1.density = self.density.copy()
        mpo1.volume = self.volume.copy()
        mpo1.N2N = self.N2N.copy()
        mpo1.N3N = self.N3N.copy()
        mpo1.NG = self.NG.copy()
        mpo1.NF = self.NF.copy()

        mpo1.N2N_sig = self.N2N_sig.copy()
        mpo1.N3N_sig = self.N3N_sig.copy()
        mpo1.NG_sig = self.NG_sig.copy()
        mpo1.NF_sig = self.NF_sig.copy()

        mpo1.power = self.power.copy()
        mpo1.power_density = self.power_density.copy()
        mpo1.flux_tally = self.flux_tally.copy()
        mpo1.mat_ID = self.mat_ID.copy()
        mpo1.mcnp_nuclides = self.mcnp_nuclides.copy()
        mpo1.ORIGEN_power = self.ORIGEN_power.copy()
        mpo1.ORIGEN_fissionpower = self.ORIGEN_fissionpower.copy()
        mpo1.ORIGEN_fluxes = self.ORIGEN_fluxes.copy()
        mpo1.temperature = self.temperature.copy()

        return mpo1

    def __add__(mpo1, mpo2):
        # this method combines mpos in a single MPO
        mpo_new = mpo1.twin()

        mpo_new.cell += mpo2.cell[:]

        for cell in mpo2.cell:
            mpo_new.concentration[cell] = mpo2.concentration[cell].copy()
            mpo_new.density[cell] = mpo2.density[cell]
            mpo_new.volume[cell] = mpo2.volume[cell]

            mpo_new.N2N[cell] = mpo2.N2N[cell].copy()
            mpo_new.N3N[cell] = mpo2.N3N[cell].copy()
            mpo_new.NG[cell] = mpo2.NG[cell].copy()
            mpo_new.NF[cell] = mpo2.NF[cell].copy()

            mpo_new.N2N_sig[cell] = mpo2.N2N_sig[cell].copy()
            mpo_new.N3N_sig[cell] = mpo2.N3N_sig[cell].copy()
            mpo_new.NG_sig[cell] = mpo2.NG_sig[cell].copy()
            mpo_new.NF_sig[cell] = mpo2.NG_sig[cell].copy()

            mpo_new.power[cell] = mpo2.power[cell]
            mpo_new.power_density[cell] = mpo2.power_density[cell]
            mpo_new.mat_ID[cell] = mpo2.mat_ID[cell][:]
            mpo_new.flux_tally[cell] = mpo2.flux_tally[cell]
            if cell in mpo2.ORIGEN_power.keys():
                mpo_new.ORIGEN_power[cell] = mpo2.ORIGEN_power[cell]
            if cell in mpo2.ORIGEN_fissionpower.keys():
                mpo_new.ORIGEN_fissionpower[
                    cell] = mpo2.ORIGEN_fissionpower[cell]
            if cell in mpo2.ORIGEN_fluxes.keys():
                mpo_new.ORIGEN_fluxes[cell] = mpo2.ORIGEN_fluxes[cell]
            if cell in mpo2.temperature.keys():
                mpo_new.temperature[cell] = mpo2.temperature[cell]
        return mpo_new


def gofiss(depletion):

    import os
    mkdir = 'mkdir moi_files'
    print mkdir
    os.system(mkdir)

    if depletion.method[0] == 'b' or depletion.method[0] == 'B':
        for iteration in range(len(depletion.time)):
            mocup(depletion, (iteration))
    elif depletion.method[0] == 'p' or depletion.method[0] == 'P':

        # generate predictor timesteps, corrector timesteps and add additional
        # flux or power timestep at BOC

        depletion.predictor_time = [depletion.time[0]/2.]
        depletion.corrector_time = [0.]
        depletion.delta_time = [0, depletion.time[0]/2.]
        for i in range(len(depletion.time) - 1):
            depletion.predictor_time.append(
                depletion.time[i] +
                0.5 *
                depletion.time[
                    i +
                    1])
            depletion.corrector_time.append(depletion.time[i])
            depletion.delta_time.append(
                0.5 *
                depletion.time[i] +
                0.5 *
                depletion.time[
                    i +
                    1])
        depletion.predictor_time.append(depletion.time[-1])
        depletion.corrector_time.append(depletion.time[-1])
        depletion.delta_time.append(depletion.time[-1]*0.5)

        if len(depletion.power) > 0:
            depletion.power = [depletion.power[0]] + depletion.power
        if len(depletion.source) > 0:
            depletion.source = [depletion.source[0]] + depletion.source
        if len(depletion.power) == 0 and len(depletion.source) == 0:
            print 'Error! Power and Flux vectors are not defined!!!'

        for iteration in range(len(depletion.predictor_time)):

            print 'starting %d iteration' % iteration

            mocup(depletion, (iteration))

    # determine EOEC keffective
    if depletion.method[0] == 'b' or depletion.method[0] == 'B':
        iteration = len(depletion.time)
    elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
        iteration = len(depletion.predictor_time)

    transportBRO(depletion, iteration)

    analysis(depletion)


def mocup(depletion, iteration):
    # This method runs MCNP and the modules of the MOCUP to produce the next MCNP input deck
    # it accepts all the depletion parameters and the iteration number

    import os

    jteration = iteration + 1

    # -------------------------------------------------------------------------
    # run neutron transport
    # ------------------------------------------------------------------------

    transportBRO(depletion, iteration)

    # -------------------------------------------------------------------------
    # post process neutron transport output
# -------------------------------------------------------------------------

    processingBRO(depletion, iteration)

# -------------------------------------------------------------------------
    # convert cross sections from MPO to ORIGEN inputdeck and execute ORIGEN
# -------------------------------------------------------------------------

    origenBRO(depletion, iteration)

# -------------------------------------------------------------------------
    # compPRO convert discharge from ORIGEN into the next MCNP inputdeck
# -------------------------------------------------------------------------

    compBRO(depletion.mpos[iteration])


def transportBRO(depletion, iteration):

    import os

    # -------------------------------------------------------------------------
    # Run MCNP
    # -------------------------------------------------------------------------

    jteration = iteration + 1

    if depletion.runTRANS == 'yes':

        if depletion.transport_module[0] == 'm' or depletion.transport_module[0] == 'M':
            mcnp = '\\rm outp%d mctal%d runtp%d' % (
                jteration, jteration, jteration)
            print mcnp
            os.system(mcnp)
            mcnp = '\cp source srctp'
            print mcnp
            os.system(mcnp)
            mcnp = 'srun mcnp5.mpi i=inp%d o=outp%d mc=mctal%d runtpe=runtp%d >> mcnplog' % (
                jteration, jteration, jteration, jteration)
            print mcnp
            os.system(mcnp)
            mcnp = '\cp srctp source'
            print mcnp
            os.system(mcnp)
            mcnp = '\\rm runtp%d srct?' % jteration
            print mcnp
            os.system(mcnp)

        elif depletion.transport_module[0] == 's' or depletion.transport_module[0] == 'S':
            serpent = '\\rm inp%d_res.m inp%d.burn' % (jteration, jteration)
            print serpent
            os.system(serpent)
            serpent = 'mpiexec sss117 inp%d' % jteration
            print serpent
            os.system(serpent)

        else:
            print 'unknown transport module'
            exit()


def processingBRO(depletion, iteration):

    # This module post processes the neutron tranposrt output files

    jteration = iteration + 1

    if depletion.runTRANS == 'yes' or len(depletion.mpos) <= iteration or len(depletion.mpos) == 0:
        # generate new mpo object based on new mcnp output
        mpo1 = mpo()
        if depletion.transport_module[0] == 'm' or depletion.transport_module[0] == 'M':
            mpo1.populate(jteration, 'yes')
        elif depletion.transport_module[0] == 's' or depletion.transport_module[0] == 'S':
            mpo1.populate_serpent(jteration, 'yes')
        depletion.mpos.append(mpo1)
    else:
        # update the concentration vector in the mpo file, but do not grab
        # information you already know (cell numbers, volumes, XSs, etc.)
        depletion.mpos[iteration].update()


def origenBRO(depletion, iteration):

    if depletion.method[0] == 'b' or depletion.method[0] == 'B':
        origenB(depletion, iteration)
    elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
        origenPC(depletion, iteration)
    else:
        print 'error depletion method must be either predictor corrector or beginning of timestep'


def origenB(depletion, iteration):

    import math
    import re
    import os

    jteration = iteration + 1
    nct = depletion.nct

    #import libraries
    dir = depletion.dir

    # import cross section library
    lib_loc = dir + depletion.library + '.lib'
    NXS_lib = open(lib_loc).read()

    # Link TAPE files to be used by origen

    # TAPE3 - blanket
    open('TAPE3.INP', 'w').write('')

    # TAPE9 - Neutron Cross Sections
    TAPE9 = 'cat %sdecay.lib %s%s.lib > TAPE9.INP' % (
        dir, dir, depletion.library)
    print TAPE9
    os.system(TAPE9)

    # TAPE10 - Gamma Cross Sectinos
    TAPE10 = 'ln -s %sgxuo2brm.lib TAPE10.INP' % (dir)
    print TAPE10
    os.system(TAPE10)

    nuclides2 = []

    for nucl in depletion.mpos[iteration].mcnp_nuclides.values():
        nuclides2.append(depletion.nct[nucl])

    power_convergence = 'no'
    power_density = {}
    fission_power = {}
    for cell in depletion.mpos[iteration].cell:
        power_density[cell] = depletion.mpos[iteration].flux_tally[
            cell]/(depletion.mpos[iteration].global_power_density*1.60217464e-19)

    while power_convergence == 'no':
        if depletion.decayheat != 'only fission' and depletion.mode != 'source':
            # update fission power source or fission neutron source based on decay heat
            # update the fission power
            print 'fission power update'
            print fission_power
            for i in range(10):
                fission_power[i] = depletion.power[
                    iteration] - depletion.global_decayheat[i]
            print fission_power
        else:
            if depletion.mode != 'source':
                fission_power = depletion.power[iteration]
            # decay heat doesn't influence convergence therefore power is
            # already converged
            power_convergence = 'yes'

        print fission_power
        print power_convergence
        # reinitialize the fullcore decay heat vector
        depletion.global_decayheat = {}

        for cell in depletion.mpos[iteration].cell:

            mpo = depletion.mpos[iteration]

            # Generate a list of nuclides in the same order as the ORIGEN cross
            # section library
            nuclides = []
            mat = material()
            mat.ORIGEN_nucl_list(lib=depletion.library)
            for nucl in (mat.structonly_nucl + mat.actin_nucl + mat.fisprod_nucl):
                if nucl in nuclides2:
                    nuclides.append(nucl)
                    if nucl in mpo.concentration[cell].keys():
                        # There is a molar mass in the concentration vector, do
                        # nothing
                        pass
                    else:
                        mpo.concentration[cell][depletion.nct[nucl]] = 0

            # identify which isotopes are not depleted so these can be added
            # back into the material when the material vectors are updated
            mpo.structural_nuclides = []
            for nucl in nuclides2:
                if nucl not in nuclides:
                    mpo.structural_nuclides.append(nucl)

            # update skeleton

            # initialize text for cross sections
            _SCROSS_SECTIONS_ = ''
            _ACROSS_SECTIONS_ = ''
            _FCROSS_SECTIONS_ = ''
            _STRUCTURAL_ = ''
            _ACTINIDES_ = ''
            _FISSIONPRODUCTS_ = ''
            _SINVENTORY_ = ''
            _AINVENTORY_ = ''
            _FINVENTORY_ = ''

            # determine nuclides followed in MOCUP based on nuclides in MPO file
            # generate list of actinides and fission products to put in origen
            # input deck

            for nuclide in nuclides:
                moles = mpo.concentration[cell][
                    nuclide]*mpo.volume[cell]/.602214078
                if nuclide in mat.actin_nucl:
                    _ACTINIDES_ += '        %s\n' % nuclide
                    _AINVENTORY_ += '2 %s %1.6E 0 0.0\n' % (nuclide, moles)
                elif nuclide in mat.fisprod_nucl:
                    _FISSIONPRODUCTS_ += '        %s\n' % nuclide
                    _FINVENTORY_ += '3 %s %1.6E 0 0.0\n' % (nuclide, moles)
                elif nuclide in mat.structonly_nucl:
                    _STRUCTURAL_ += '        %s\n' % nuclide
                    _SINVENTORY_ += '1 %s %1.6E 0 0.0\n' % (nuclide, moles)

            # combine inventory from actinides and fission products
            _INVENTORY_ = _SINVENTORY_ + _AINVENTORY_ + _FINVENTORY_

            for nuclide in nuclides:
                # determine if fuel or fission product and which library to
                # append to and cross section format
                k = NXS_lib.split().index(nuclide)
                if nuclide in mat.actin_nucl:
                    lib_num = mat.actinlib_num
                elif nuclide in mat.fisprod_nucl:
                    lib_num = mat.fisprodlib_num
                elif nuclide in mat.structonly_nucl:
                    lib_num = mat. structlib_num

                # determine metastable production ratio
                token = '      '
                token = token[:-len(nuclide)] + nuclide
                k = NXS_lib.find(token)

                # initialize fraction of NG XS and N2N XS that end in ground
                # state and metastable states
                fNG = 1
                fNGex = 0
                fN2N = 1
                fN2Nex = 0

                # intialize cross sections
                NG = 0
                N2N = 0
                NGex = 0
                N2Nex = 0

                # grab cross sections from XS library to allocate the NG  and
                # N2N cross section to reactions that end in ground state and
                # metastable states
                if NXS_lib[k+14:k+16] != '  ':
                    NG = float(NXS_lib[k+7:k+12]) * \
                        math.pow(10, int(NXS_lib[k+13:k+16]))
                if NXS_lib[k+24:k+26] != '  ':
                    N2N = float(NXS_lib[k+17:k+22]) * \
                        math.pow(10, int(NXS_lib[k+23:k+26]))
                if NXS_lib[k+54:k+56] != '  ':
                    NGex = float(NXS_lib[k+47:k+52]) * \
                        math.pow(10, int(NXS_lib[k+53:k+56]))
                if NXS_lib[k+64:k+66] != '  ':
                    N2Nex = float(NXS_lib[k+57:k+62]) * \
                        math.pow(10, int(NXS_lib[k+63:k+66]))

                # determime fraction of reaction that end in ground state and
                # metastable state
                if NG + NGex != 0:
                    fNG = NG/(NG + NGex)
                    fNGex = NGex/(NG + NGex)
                if N2N + N2Nex != 0:
                    fN2N = N2N/(N2N+N2Nex)
                    fN2Nex = N2Nex/(N2N+N2Nex)

                # Initialize Cross Sections at 0
                SN3N = 0
                SN2N = 0
                SNF = 0
                SNG = 0

                # define cross sections
                if nuclide in mpo.N3N[cell].keys():
                    SN3N = mpo.N3N[cell][nuclide]
                if nuclide in mpo.N2N[cell].keys():
                    SN2N = mpo.N2N[cell][nuclide]
                if nuclide in mpo.NF[cell].keys():
                    SNF = mpo.NF[cell][nuclide]
                if nuclide in mpo.NG[cell].keys():
                    SNG = mpo.NG[cell][nuclide]

                # write cross sections for a given nuclide

                # determine if a nuclide is actinide or fission product
                if nuclide in mat.actin_nucl:
                    # actinide
                    # write line append to skeleton
                    _ACROSS_SECTIONS_ += '%s %s %1.3E %1.3E %1.3E %1.3E %1.3E %1.3E -1.000000\n' % (
                        lib_num, nuclide, (fNG*SNG), (fN2N*SN2N), SN3N, SNF, (fNGex*SNG), (fN2Nex*SN2N))
                elif nuclide in mat.fisprod_nucl:
                    token = '%s  %s' % (lib_num, nuclide)
                    k = NXS_lib.find(token)
                    l = NXS_lib[k:].find('\n')
                    m = NXS_lib[k+l+1:].find('\n')
                    _XS_ = NXS_lib[k:(k+l+m+2)]
                    # replace (N,G)
                    SNGs = '%1.3E' % (SNG*fNG)
                    _XS_ = _XS_[:12]+SNGs+_XS_[21:]
                    SNGs = '%1.3E' % (SN2N*fN2N)
                    _XS_ = _XS_[:22]+SNGs+_XS_[31:]
                    SNGs = '%1.3E' % (SNG*fNGex)
                    _XS_ = _XS_[:52]+SNGs+_XS_[61:]
                    SNGs = '%1.3E' % (SN2N*fN2Nex)
                    _XS_ = _XS_[:62]+SNGs+_XS_[71:]
                    _FCROSS_SECTIONS_ += _XS_
                elif nuclide in mat.structonly_nucl:
                    token = '%s        ' % (lib_num)
                    token = token[:-len(nuclide)] + nuclide
                    k = NXS_lib.find(token)
                    l = NXS_lib[k:].find('\n')
                    _XS_ = NXS_lib[k:(k+l+1)]
                    SNGs = '%1.3E' % (SNG*fNG)
                    _XS_ = _XS_[:12]+SNGs+_XS_[21:]
                    SNGs = '%1.3E' % (SN2N*fN2N)
                    _XS_ = _XS_[:22]+SNGs+_XS_[31:]
                    SNGs = '%1.3E' % (SNG*fNGex)
                    _XS_ = _XS_[:52]+SNGs+_XS_[61:]
                    SNGs = '%1.3E' % (SN2N*fN2Nex)
                    _XS_ = _XS_[:62]+SNGs+_XS_[71:]
                    _SCROSS_SECTIONS_ += _XS_

            actinide_num = mat.actinlib_num
            fissionproduct_num = mat.fisprodlib_num
            structural_num = mat.structlib_num

            # combine cross sections from actinides and fission products
            _CROSS_SECTIONS_ = _SCROSS_SECTIONS_ + \
                _ACROSS_SECTIONS_+_FCROSS_SECTIONS_

            # write skeleton file as TAPE5 and moi_files/moi.cell.iteration
            skele = ' -1'
            skele += '\n -1'
            skele += '\n -1'
            skele += '\n TIT    %s' % depletion.title
            skele += '\n BAS    BWR bundle'
            skele += '\n LIP    1 1 0 '
            if len(_STRUCTURAL_) > 0:
                skele += '\n LPU \n' + _STRUCTURAL_ + '        -1'
            skele += '\n LPU \n'
            skele += _ACTINIDES_
            skele += '        -1'
            skele += '\n LPU \n' + _FISSIONPRODUCTS_ + '        -1'
            if len(_STRUCTURAL_) > 0:
                skele += '\n LIB    0   1 2 3    -%s -%s -%s   9 50 0  4   0' % (
                    structural_num, actinide_num, fissionproduct_num)
            else:
                skele += '\n LIB    0   0 2 3    0 -%s -%s   9 50 0  4   0' % (
                    actinide_num, fissionproduct_num)
            skele += '\n OPTL   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n OPTA   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n OPTF   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n CUT    3 1.0E-24  28 1.0E-75   -1'
            skele += '\n INP    1   -2  -1  -1  1  1'
            skele += '\n BUP'
            skele += '\n _MODE_      _TIME1_   _POWER0_    1   2   4   2'
            skele += '\n _MODE_      _TIME2_   _POWER1_    2   3   4   0'
            skele += '\n _MODE_      _TIME3_   _POWER2_    3   4   4   0'
            skele += '\n _MODE_      _TIME4_   _POWER3_    4   5   4   0'
            skele += '\n _MODE_      _TIME5_   _POWER4_    5   6   4   0'
            skele += '\n _MODE_      _TIME6_   _POWER5_    6   7   4   0'
            skele += '\n _MODE_      _TIME7_   _POWER6_    7   8   4   0'
            skele += '\n _MODE_      _TIME8_   _POWER7_    8   9   4   0'
            skele += '\n _MODE_      _TIME9_   _POWER8_    9   10   4   0'
            skele += '\n _MODE_      _TIME10_   _POWER9_   10   11   4   0'
            skele += '\n BUP'
            skele += '\n OUT    11  1  0  0'
            skele += '\n PCH    11  11  11'
            skele += '\n STP    4'
            if len(_STRUCTURAL_) > 0:
                skele += '\n' + _CROSS_SECTIONS_+_INVENTORY_+'0'
            else:
                skele += '\n' + _ACROSS_SECTIONS_+_FCROSS_SECTIONS_ + \
                    _AINVENTORY_ + _FINVENTORY_ + '0'

            # write the mode in the origen input file
            # write power or flux in origen input file depending on mode
            if depletion.mode == 'power':
                # origen must be run in power mode
                skele = skele.replace('_MODE_', 'IRP')
                for i in range(10):
                    token = '_POWER%d_' % i
                    power_string = '%1.5E' % (fission_power[i])
                    skele = skele.replace(token, power_string)

            elif depletion.mode == 'source' or depletion.mode == 'flux':
                # origen must be run in flux mode
                skele = skele.replace('_MODE_', 'IRF')
                if depletion.mode == 'source':
                    for i in range(10):
                        token = '_POWER%d_' % i
                        power_string = '%1.5E' % (
                            depletion.source[iteration]*mpo.flux_tally[cell])
                        skele = skele.replace(token, power_string)
                else:
                    for i in range(10):
                        token = '_POWER%d_' % i
                        power_string = '%1.5E' % (
                            fission_power[i]*power_density[cell])
                        skele = skele.replace('_POWER_', power_string)
            else:
                print "ERROR! depletion mode must either be 'power' or 'flux'"
                print 'depletion mode is %s' % depletion.mode
                exit()

            # write the subtime steps in ORIGEN input file
            time_sum = 0
            for j in range(10):
                time_sum += depletion.time[iteration]/10
                dummy = '_TIME%d_' % (j + 1)
                time_string = '%1.5E' % time_sum
                skele = skele.replace(dummy, time_string)

            if depletion.start_ocf == 1 and jteration == 1:
                # MAKE OCF file from composition file
                mat = material()
                mat.comp = mpo.concentration[cell]
                mat = mat*(mpo.volume[cell]/.602214078)
                ocf_loc = 'moi_files/moi.%s.0.pch' % cell
                mat.make_ocf(ocf_loc, lib=depletion.library)
            # copy discharge composition from i - 1 cycle to charge of i cycle
            copy = '\cp moi_files/moi.%s.%d.pch TAPE4.INP' % (cell, iteration)
            print copy
            os.system(copy)

            open('TAPE5.INP', 'w').write(skele)
            # execute origen

            if depletion.library[0:3] in ['amo', 'emo', 'fft']:
                # origen uses a fast reactor set of cross sections, thus origen
                # will be executed with o2_fast
                origen = 'o2_fast'
            else:
                # origen uses a fast reactor set of cross sections, thus origen
                # will be executed with o2_therm
                origen = 'o2_therm'

            print origen
            os.system(origen)

            TAPE6 = open('TAPE6.OUT').read()
            TAPE12 = open('TAPE12.OUT').read()
            tokens = re.compile('\s\d+\s{13,13}THERMAL').findall(TAPE12)

            k = TAPE6.find('NEUT. FLUX')
            fluxes = TAPE6[k:].split()[3:13]

            k = TAPE6.find('SP POW,MW')
            powers = TAPE6[k:].split()[3:13]
            depletion.mpos[iteration].ORIGEN_power[cell] = 0
            depletion.mpos[iteration].ORIGEN_fissionpower[cell] = 0
            depletion.mpos[iteration].ORIGEN_fluxes[cell] = 0

            decay_heat = []
            substep_decay_heat = {}
            for token in tokens:
                skeleton_token = 'PAGE     '
                token1 = skeleton_token[
                    :-len(token.split()[0])] + token.split()[0]
                k = TAPE6.find(token1)
                l = TAPE6[k:].find('0TOTAL')
                decay_heat.append(TAPE6[k+l:].split()[1:12])

            ave_decay_heat = 0
            for heat_list in decay_heat:
                for i in range(len(heat_list)):
                    if (i*depletion.time[iteration]) in substep_decay_heat.keys():
                        substep_decay_heat[i] += float(heat_list[i])*1e-6
                    else:
                        substep_decay_heat[i] = float(heat_list[i])*1e-6
                    if i == 0 or i == range(len(heat_list))[-1]:
                        ave_decay_heat += 1. / \
                            (2.*len(powers))*float(heat_list[i])*1e-6
                    else:
                        ave_decay_heat += 1. / \
                            (len(powers))*float(heat_list[i])*1e-6

            timepower = 0
            for i in range(len(powers)):
                # transform substep_deacay_heat from substep to EFPD*MW
                timepower += (
                    (substep_decay_heat[i] + substep_decay_heat[i + 1]) / 2. +
                    float(powers[i])) * depletion.time[iteration] / len(powers)
                substep_decay_heat[timepower] = substep_decay_heat[i]
                del substep_decay_heat[i]
                # add to power and flux attributes
                depletion.mpos[iteration].ORIGEN_fissionpower[
                    cell] += float(powers[i])/len(powers)
                depletion.mpos[iteration].ORIGEN_fluxes[
                    cell] += float(fluxes[i])/len(fluxes)
                print depletion.mpos[iteration].ORIGEN_fluxes[cell]
            depletion.mpos[iteration].ORIGEN_power[cell] = depletion.mpos[
                iteration].ORIGEN_fissionpower[cell] + ave_decay_heat

            for i in range(10):
                if i in depletion.global_decayheat.keys():
                    depletion.global_decayheat[
                        i] += 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6
                else:
                    depletion.global_decayheat[
                        i] = 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6

            # move outputs in moi_files

            move = 'cat TAPE12.OUT TAPE6.OUT > moi_files/moi.%s.%d.out' % (
                cell, jteration)
            print move
            os.system(move)

            # move origen input deck, TAPE5.INP
            move = '\mv TAPE5.INP moi_files/moi.%s.%d.inp' % (cell, jteration)
            print move
            os.system(move)

            # move discharge composition, TAPE7.OUT
            move = '\mv TAPE7.OUT moi_files/moi.%s.%d.pch' % (cell, jteration)
            print move
            os.system(move)

        # update the estimate of the decay heat fraction
        sum_power = 0
        sum_fissionpower = 0
        for cell in depletion.mpos[iteration].cell:
            sum_power = depletion.mpos[iteration].ORIGEN_power[cell]
            sum_fissionpower = depletion.mpos[
                iteration].ORIGEN_fissionpower[cell]

        # check for convergence of power
        if power_convergence == 'no':

            print '\n\n\n'
            print sum_power
            print depletion.power[iteration]
            print '\n\n\n'

            if math.fabs(1. - sum_power / depletion.power[iteration]) <= depletion.power_err:
                # the calculated power is within error of the target power
                print '\n\n\n\npower converged!\n\n\n\n'
                power_convergence = 'yes'

                # clean workspace
                remove = '\\rm TAPE*'
                print remove
                os.system(remove)

                # update global decayheat vector
                for i in range(10):
                    depletion.global_decayheat[i] = sum_power - sum_fissionpower
            else:
                # the calculated power is not within error of the target
                print '\n\n\n\npower not converged\n\n\n\n'
                power_convergence = 'no'

        else:
            print '\n\n\n\npower converged!\n\n\n\n'
            # clean workspace
            print remove
            os.system(remove)


def origenPC(depletion, iteration):

    # This is the origen driver for mocup in predictor corrector mode

    import math
    import os
    import re

    jteration = iteration + 1
    nct = depletion.nct

    #import libraries
    dir = depletion.dir

    # import cross section library
    lib_loc = dir + depletion.library + '.lib'
    NXS_lib = open(lib_loc).read()

    # Link TAPE files to be used by origen

    # TAPE3 - blanket
    open('TAPE3.INP', 'w').write('')

    # TAPE9 - Neutron Cross Sections
    TAPE9 = 'cat %sdecay.lib %s%s.lib > TAPE9.INP' % (
        dir, dir, depletion.library)
    print TAPE9
    os.system(TAPE9)

    # TAPE10 - Gamma Cross Sectinos
    TAPE10 = 'ln -s %sgxuo2brm.lib TAPE10.INP' % (dir)
    print TAPE10
    os.system(TAPE10)

    nuclides2 = []
    for nucl in depletion.mpos[iteration].mcnp_nuclides.values():
        nuclides2.append(depletion.nct[nucl])

    power_convergence = 'no'
    fission_power = {}
    power_density = {}
    for cell in depletion.mpos[iteration].cell:
        power_density[cell] = depletion.mpos[iteration].flux_tally[
            cell]/(depletion.mpos[iteration].global_power_density*1.60217464e-19)

    while power_convergence == 'no':
        if depletion.decayheat != 'only fission' and depletion.mode != 'source':
            # update fission power source or fission neutron source based on decay heat
            # update the fission power
            for i in depletion.global_decayheat.keys():
                fission_power[i] = depletion.power[
                    iteration] - depletion.global_decayheat[i]
            print fission_power

        else:
            if depletion.mode != 'source':
                fission_power = depletion.power[iteration]
            # decay heat doesn't influence convergence therefore power is
            # already converged
            power_convergence = 'yes'

        # intiialize new estimate of the global decay heat
        depletion.global_decayheat = {}

        for cell in depletion.mpos[iteration].cell:

            mpo = depletion.mpos[iteration]

            # Generate a list of nuclides in the same order as the ORIGEN cross
            # section library
            nuclides = []
            mat = material()
            mat.ORIGEN_nucl_list(lib=depletion.library)
            for nucl in (mat.structonly_nucl + mat.actin_nucl + mat.fisprod_nucl):
                if nucl in nuclides2:
                    nuclides.append(nucl)
                    if nucl in mpo.concentration[cell].keys():
                        # There is a molar mass in the concentration vector, do
                        # nothing
                        pass
                    else:
                        mpo.concentration[cell][nucl] = 0

            # update skeleton

            # initialize text for cross sections
            _SCROSS_SECTIONS_ = ''
            _ACROSS_SECTIONS_ = ''
            _FCROSS_SECTIONS_ = ''
            _STRUCTURAL_ = ''
            _ACTINIDES_ = ''
            _FISSIONPRODUCTS_ = ''
            _SINVENTORY_ = ''
            _AINVENTORY_ = ''
            _FINVENTORY_ = ''

            # determine nuclides followed in MOCUP based on nuclides in MPO file
            # generate list of actinides and fission products to put in origen
            # input deck

            for nuclide in nuclides:
                moles = mpo.concentration[cell][
                    nuclide]*mpo.volume[cell]/.602214078
                if nuclide in mat.actin_nucl:
                    _ACTINIDES_ += '        %s\n' % nuclide
                    _AINVENTORY_ += '2 %s %1.6E 0 0.0\n' % (nuclide, moles)
                elif nuclide in mat.fisprod_nucl:
                    _FISSIONPRODUCTS_ += '        %s\n' % nuclide
                    _FINVENTORY_ += '3 %s %1.6E 0 0.0\n' % (nuclide, moles)
                elif nuclide in mat.structonly_nucl:
                    _STRUCTURAL_ += '        %s\n' % nuclide
                    _SINVENTORY_ += '1 %s %1.6E 0 0.0\n' % (nuclide, moles)

            # combine inventory from actinides and fission products
            _INVENTORY_ = _SINVENTORY_ + _AINVENTORY_ + _FINVENTORY_

            for nuclide in nuclides:
                # determine if fuel or fission product and which library to
                # append to and cross section format
                k = NXS_lib.split().index(nuclide)
                if nuclide in mat.actin_nucl:
                    lib_num = mat.actinlib_num
                elif nuclide in mat.fisprod_nucl:
                    lib_num = mat.fisprodlib_num
                elif nuclide in mat.structonly_nucl:
                    lib_num = mat. structlib_num

                # determine metastable production ratio
                token = '      '
                token = token[:-len(nuclide)] + nuclide
                k = NXS_lib.find(token)

                # initialize fraction of NG XS and N2N XS that end in ground
                # state and metastable states
                fNG = 1
                fNGex = 0
                fN2N = 1
                fN2Nex = 0

                # intialize cross sections
                NG = 0
                N2N = 0
                NGex = 0
                N2Nex = 0

                # grab cross sections from XS library to allocate the NG  and
                # N2N cross section to reactions that end in ground state and
                # metastable states
                if NXS_lib[k+14:k+16] != '  ':
                    NG = float(NXS_lib[k+7:k+12]) * \
                        math.pow(10, int(NXS_lib[k+13:k+16]))
                if NXS_lib[k+24:k+26] != '  ':
                    N2N = float(NXS_lib[k+17:k+22]) * \
                        math.pow(10, int(NXS_lib[k+23:k+26]))
                if NXS_lib[k+54:k+56] != '  ':
                    NGex = float(NXS_lib[k+47:k+52]) * \
                        math.pow(10, int(NXS_lib[k+53:k+56]))
                if NXS_lib[k+64:k+66] != '  ':
                    N2Nex = float(NXS_lib[k+57:k+62]) * \
                        math.pow(10, int(NXS_lib[k+63:k+66]))

                # determime fraction of reaction that end in ground state and
                # metastable state
                if NG + NGex != 0:
                    fNG = NG/(NG + NGex)
                    fNGex = NGex/(NG + NGex)
                if N2N + N2Nex != 0:
                    fN2N = N2N/(N2N+N2Nex)
                    fN2Nex = N2Nex/(N2N+N2Nex)

                # Initialize Cross Sections at 0
                SN3N = 0
                SN2N = 0
                SNF = 0
                SNG = 0

                # define cross sections
                if nuclide in mpo.N3N[cell].keys():
                    SN3N = mpo.N3N[cell][nuclide]
                if nuclide in mpo.N2N[cell].keys():
                    SN2N = mpo.N2N[cell][nuclide]
                if nuclide in mpo.NF[cell].keys():
                    SNF = mpo.NF[cell][nuclide]
                if nuclide in mpo.NG[cell].keys():
                    SNG = mpo.NG[cell][nuclide]

                # write cross sections for a given nuclide

                # determine if a nuclide is actinide or fission product
                if nuclide in mat.actin_nucl:
                    # actinide
                    # write line append to skeleton
                    _ACROSS_SECTIONS_ += '%s %s %1.3E %1.3E %1.3E %1.3E %1.3E %1.3E -1.000000\n' % (
                        lib_num, nuclide, (fNG*SNG), (fN2N*SN2N), SN3N, SNF, (fNGex*SNG), (fN2Nex*SN2N))
                    actinide_num = lib_num
                elif nuclide in mat.fisprod_nucl:
                    token = '%s  %s' % (lib_num, nuclide)
                    k = NXS_lib.find(token)
                    l = NXS_lib[k:].find('\n')
                    m = NXS_lib[k+l+1:].find('\n')
                    _XS_ = NXS_lib[k:(k+l+m+2)]
                    # replace (N,G)
                    SNGs = '%1.3E' % (SNG*fNG)
                    _XS_ = _XS_[:12]+SNGs+_XS_[21:]
                    SNGs = '%1.3E' % (SN2N*fN2N)
                    _XS_ = _XS_[:22]+SNGs+_XS_[31:]
                    SNGs = '%1.3E' % (SNG*fNGex)
                    _XS_ = _XS_[:52]+SNGs+_XS_[61:]
                    SNGs = '%1.3E' % (SN2N*fN2Nex)
                    _XS_ = _XS_[:62]+SNGs+_XS_[71:]
                    _FCROSS_SECTIONS_ += _XS_
                elif nuclide in mat.structonly_nucl:
                    token = '%s        ' % (lib_num)
                    token = token[:-len(nuclide)] + nuclide
                    k = NXS_lib.find(token)
                    l = NXS_lib[k:].find('\n')
                    _XS_ = NXS_lib[k:(k+l+1)]
                    SNGs = '%1.3E' % (SNG*fNG)
                    _XS_ = _XS_[:12]+SNGs+_XS_[21:]
                    SNGs = '%1.3E' % (SN2N*fN2N)
                    _XS_ = _XS_[:22]+SNGs+_XS_[31:]
                    SNGs = '%1.3E' % (SNG*fNGex)
                    _XS_ = _XS_[:52]+SNGs+_XS_[61:]
                    SNGs = '%1.3E' % (SN2N*fN2Nex)
                    _XS_ = _XS_[:62]+SNGs+_XS_[71:]
                    _SCROSS_SECTIONS_ += _XS_

            actinide_num = mat.actinlib_num
            fissionproduct_num = mat.fisprodlib_num
            structural_num = mat.structlib_num

            # combine cross sections from actinides and fission products
            _CROSS_SECTIONS_ = _SCROSS_SECTIONS_ + \
                _ACROSS_SECTIONS_+_FCROSS_SECTIONS_

            # write skeleton file as TAPE5 and moi_files/moi.cell.iteration
            skele = ' -1'
            skele += '\n -1'
            skele += '\n -1'
            skele += '\n TIT    %s' % depletion.title
            skele += '\n BAS    BWR bundle'
            skele += '\n LIP    1 1 0 '
            if len(_STRUCTURAL_) > 0:
                skele += '\n LPU \n' + _STRUCTURAL_ + '        -1'
            skele += '\n LPU \n'
            skele += _ACTINIDES_
            skele += '        -1'
            skele += '\n LPU \n' + _FISSIONPRODUCTS_ + '        -1'
            if len(_STRUCTURAL_) > 0:
                skele += '\n LIB    0   1 2 3    -%s -%s -%s   9 50 0  4   0' % (
                    structural_num, actinide_num, fissionproduct_num)
            else:
                skele += '\n LIB    0   0 2 3    0 -%s -%s   9 50 0  4   0' % (
                    actinide_num, fissionproduct_num)
            skele += '\n OPTL   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n OPTA   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n OPTF   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
            skele += '\n CUT    3 1.0E-24  28 1.0E-75   -1'
            skele += '\n INP    1   -2  -1  -1  1  1'
            skele += '\n BUP'
            skele += '\n _MODE_      _TIME1_   _POWER0_    1   2   4   2'
            skele += '\n _MODE_      _TIME2_   _POWER1_    2   3   4   0'
            skele += '\n _MODE_      _TIME3_   _POWER2_    3   4   4   0'
            skele += '\n _MODE_      _TIME4_   _POWER3_    4   5   4   0'
            skele += '\n _MODE_      _TIME5_   _POWER4_    5   6   4   0'
            skele += '\n _MODE_      _TIME6_   _POWER5_    6   7   4   0'
            skele += '\n _MODE_      _TIME7_   _POWER6_    7   8   4   0'
            skele += '\n _MODE_      _TIME8_   _POWER7_    8   9   4   0'
            skele += '\n _MODE_      _TIME9_   _POWER8_    9   10   4   0'
            skele += '\n _MODE_      _TIME10_   _POWER9_   10   11   4   0'
            skele += '\n BUP'
            skele += '\n OUT    11  1  0  0'
            skele += '\n PCH    11  11  11'
            skele += '\n STP    4'
            if len(_STRUCTURAL_) > 0:
                skele += '\n' + _CROSS_SECTIONS_+_INVENTORY_+'0'
            else:
                skele += '\n' + _ACROSS_SECTIONS_ + \
                    _FCROSS_SECTIONS_ + _AINVENTORY_ + _FINVENTORY_

            # write the mode in the origen input file
            # write power or flux in origen input file depending on mode
            if depletion.mode == 'power':
                # origen must be run in power mode
                for i in range(10):
                    skele = skele.replace('_MODE_', 'IRP')
                    power_string = '%1.5E' % (fission_power[i]*mpo.power[cell])
                    token = '_POWER%d_' % (i)
                    skele = skele.replace(token, power_string)

            elif depletion.mode == 'source' or depletion.mode == 'flux':
                # origen must be run in flux mode
                skele = skele.replace('_MODE_', 'IRF')
                for i in range(10):
                    if depletion.mode == 'source':
                        power_string = '%1.5E' % (
                            depletion.source[iteration]*power_density[cell])
                    else:
                        power_string = '%1.5E' % (
                            fission_power[i]*mpo.flux_tally[cell]/(mpo.global_power_density*1.60217646e-19))
                    token = '_POWER%d_' % (i)
                    skele = skele.replace(token, power_string)
            else:
                print "ERROR! depletion mode must either be 'power' or 'flux'"
                print 'depletion mode is %s' % depletion.mode
                exit()

            if depletion.start_ocf == 1 and jteration == 1:
                # MAKE OCF file from composition file
                mat = material()
                mat.comp = mpo.concentration[cell]
                mat = mat*(mpo.volume[cell]/.602214078)
                ocf_loc = 'moi_files/moi.%s.0.pch' % cell
                mat.make_ocf(ocf_loc, lib=depletion.library)
                ocf_loc = 'moi_files/cor.%s.0.pch' % cell
                mat.make_ocf(ocf_loc, lib=depletion.library)

            elif depletion.start_ocf == 0 and jteration == 1:
                # ensure there is a BOC composition in the moi_files director
                # for the predictor and corrector step
                moi_files_dir = os.listdir('moi_files')
                ocf_loc = 'moi.%s.0.pch' % cell
                cor_loc = 'cor.%s.0.pch' % cell
                if ocf_loc in moi_files_dir:
                    # there is a BOC ocf file in the moi_files directory
                    if cor_loc not in moi_files_dir:
                        copy = 'cp moi_files/%s moi_files/%s' % (ocf_loc,
                                                                 cor_loc)
                        print copy
                        os.system(copy)
                else:
                    # exit analysis because there is not BOC compsition
                    print 'there is no BOC compositoin ocf file, %s, in the moi_files directory' % ocf_loc

            # copy discharge composition from i - 1 cycle to charge of i cycle
            copy = '\cp moi_files/cor.%s.%d.pch TAPE4.INP' % (cell, iteration)
            print copy
            os.system(copy)

            # write the subtime steps in ORIGEN input file
            p_time_sum = 0
            c_time_sum = 0
            predictor = skele[:]
            corrector = skele[:]
            for j in range(10):
                p_time_sum += depletion.predictor_time[iteration]/10.
                c_time_sum += depletion.corrector_time[iteration]/10.
                dummy = '_TIME%d_' % (j + 1)
                p_time_string = '%1.5E' % p_time_sum
                c_time_string = '%1.5E' % c_time_sum
                predictor = predictor.replace(dummy, p_time_string)
                corrector = corrector.replace(dummy, c_time_string)

            open('TAPE5.INP', 'w').write(predictor)
            # execute origen

            if depletion.library[0:3] in ['amo', 'emo', 'fft']:
                # origen uses a fast reactor set of cross sections, thus origen
                # will be executed with o2_fast
                origen = 'o2_fast'
            else:
                # origen uses a fast reactor set of cross sections, thus origen
                # will be executed with o2_therm
                origen = 'o2_therm'

            print origen
            os.system(origen)

            # move outputs in moi_files PREDICTOR

            move = 'cat TAPE12.OUT TAPE6.OUT > moi_files/moi.%s.%d.out' % (
                cell, jteration)
            print move
            os.system(move)

            # move origen input deck, TAPE5.INP
            move = '\mv TAPE5.INP moi_files/moi.%s.%d.inp' % (cell, jteration)
            print move
            os.system(move)

            # move discharge composition, TAPE7.OUT
            move = '\mv TAPE7.OUT moi_files/moi.%s.%d.pch' % (cell, jteration)
            print move
            os.system(move)

            # execute origen for corrector step
            if depletion.corrector_time[iteration] != 0:
                open('TAPE5.INP', 'w').write(corrector)

                print origen
                os.system(origen)
                TAPE6 = open('TAPE6.OUT').read()

                k = TAPE6.find('NEUT. FLUX')
                fluxes = TAPE6[k:].split()[3:13]

                k = TAPE6.find('SP POW,MW')
                powers = TAPE6[k:].split()[3:13]
                depletion.mpos[iteration].ORIGEN_power[cell] = 0
                depletion.mpos[iteration].ORIGEN_fissionpower[cell] = 0
                depletion.mpos[iteration].ORIGEN_fluxes[cell] = 0

                TAPE12 = open('TAPE12.OUT').read()
                tokens = re.compile('\s\d+\s{13,13}THERMAL').findall(TAPE12)

                decay_heat = []

                substep_decay_heat = {}

                for token in tokens:

                    skeleton_token = 'PAGE     '
                    token1 = skeleton_token[
                        :-len(token.split()[0])] + token.split()[0]

                    k = TAPE6.find(token1)
                    l = TAPE6[k:].find('0TOTAL')

                    decay_heat.append(TAPE6[k+l:].split()[1:12])

                ave_decay_heat = 0
                for heat_list in decay_heat:
                    for i in range(len(heat_list)):
                        if (i*depletion.corrector_time[iteration]) in substep_decay_heat.keys():
                            substep_decay_heat[i] += float(heat_list[i])*1e-6
                        else:
                            substep_decay_heat[i] = float(heat_list[i])*1e-6

                        if i == 0 or i == range(len(heat_list))[-1]:
                            ave_decay_heat += 1. / \
                                (2.*len(powers))*float(heat_list[i])*1e-6
                        else:
                            ave_decay_heat += 1. / \
                                (len(powers))*float(heat_list[i])*1e-6
                timepower = 0
                for i in range(len(powers)):
                    # transform substep_deacay_heat from substep to EFPD*MW
                    timepower += ((substep_decay_heat[i] + substep_decay_heat[i+1])/2. +
                                  float(powers[i]))*depletion.corrector_time[iteration]/len(powers)

                    substep_decay_heat[timepower] = substep_decay_heat[i]
                    del substep_decay_heat[i]
                    # add to power and flux attributes
                    depletion.mpos[iteration].ORIGEN_fissionpower[
                        cell] += float(powers[i])/len(powers)
                    depletion.mpos[iteration].ORIGEN_fluxes[
                        cell] += float(fluxes[i])/len(fluxes)
                    if i in depletion.global_decayheat.keys():
                        depletion.global_decayheat[
                            i] += 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6
                    else:
                    #	global_fissionpower_density[i] = float(powers[i])
                        depletion.global_decayheat[
                            i] = 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6

                depletion.mpos[iteration].ORIGEN_power[cell] = depletion.mpos[
                    iteration].ORIGEN_fissionpower[cell] + ave_decay_heat
                power_density[cell] = depletion.mpos[iteration].ORIGEN_fluxes[
                    cell]/depletion.mpos[iteration].ORIGEN_fissionpower[cell]

                # move outputs in moi_files

                # move outputs in moi_files CORRECTOR
                move = 'cat TAPE12.OUT TAPE6.OUT > moi_files/cor.%s.%d.out' % (
                    cell, jteration)
                print move
                os.system(move)

                # move origen input deck, TAPE5.INP
                move = '\mv TAPE5.INP moi_files/cor.%s.%d.inp' % (
                    cell, jteration)
                print move
                os.system(move)

                # move CORRECTOR discharge composition, TAPE7.OUT
                move = '\mv TAPE7.OUT moi_files/cor.%s.%d.pch' % (
                    cell, jteration)
                print move
                os.system(move)

            else:
                # This is the first step where the BOC composition is the charge

                if depletion.start_ocf == 1 and jteration == 1:
                    # No initial timestep is expected, thus the concentration
                    # taken from the MCNP outp file is converted to moles and
                    # used as the charge for the next timestep

                    TAPE6 = open('TAPE6.OUT').read()
                    TAPE12 = open('TAPE12.OUT').read()
                    tokens = re.compile('\s\d+\s{13,13}THERMAL').findall(TAPE12)

                    k = TAPE6.find('NEUT. FLUX')
                    fluxes = TAPE6[k:].split()[3:13]

                    k = TAPE6.find('SP POW,MW')
                    powers = TAPE6[k:].split()[3:13]
                    depletion.mpos[iteration].ORIGEN_power[cell] = 0
                    depletion.mpos[iteration].ORIGEN_fissionpower[cell] = 0
                    depletion.mpos[iteration].ORIGEN_fluxes[cell] = 0

                    decay_heat = []
                    substep_decay_heat = {}

                    for token in tokens:
                        skeleton_token = 'PAGE     '
                        token1 = skeleton_token[
                            :-len(token.split()[0])] + token.split()[0]

                        k = TAPE6.find(token1)
                        l = TAPE6[k:].find('0TOTAL')

                        decay_heat.append(TAPE6[k+l:].split()[1:12])

                    ave_decay_heat = 0
                    for heat_list in decay_heat:
                        for i in range(len(heat_list)):
                            if i in substep_decay_heat.keys():
                                substep_decay_heat[
                                    i] += float(heat_list[i])*1e-6
                            else:
                                substep_decay_heat[i] = float(heat_list[i])*1e-6

                            if i == 0 or i == range(len(heat_list))[-1]:
                                ave_decay_heat += 1. / \
                                    (2.*len(powers))*float(heat_list[i])*1e-6
                            else:
                                ave_decay_heat += 1. / \
                                    (len(powers))*float(heat_list[i])*1e-6

                    timepower = 0
                    for i in range(len(powers)):
                        # transform substep_deacay_heat from substep to EFPD*MW
                        timepower += ((substep_decay_heat[i] + substep_decay_heat[i+1])/2. + float(
                            powers[i]))*depletion.corrector_time[iteration]/len(powers)
                        substep_decay_heat[timepower] = substep_decay_heat[i]
                        del substep_decay_heat[i]
                        # add to power and flux attributes
                        depletion.mpos[iteration].ORIGEN_fissionpower[
                            cell] += float(powers[i])/len(powers)
                        depletion.mpos[iteration].ORIGEN_fluxes[
                            cell] += float(fluxes[i])/len(fluxes)
                        if i in depletion.global_decayheat.keys():
                            depletion.global_decayheat[
                                i] += 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6
                        else:
                            depletion.global_decayheat[
                                i] = 1./2.*(float(heat_list[i]) + float(heat_list[i+1]))*1e-6

                    depletion.mpos[iteration].ORIGEN_power[cell] = depletion.mpos[
                        iteration].ORIGEN_fissionpower[cell] + ave_decay_heat
                    power_density[cell] = depletion.mpos[iteration].ORIGEN_fluxes[
                        cell]/depletion.mpos[iteration].ORIGEN_fissionpower[cell]

                    mat = material()
                    mat.comp = depletion.mpos[iteration].concentration[cell]
                    mat = mat * \
                        (depletion.mpos[iteration].volume[cell]/.602214078)
                    mat_loc = 'moi_files/cor.%s.%d.pch' % (cell, jteration)
                    mat.make_ocf(mat_loc, lib=depletion.library)

                elif depletion.start_ocf == 1 and jteration != 1:
                    print 'ERROR'
                else:
                    # There is a zero length corrector time step and mocup
                    # expects a written pch for beginning of this timestep, thus
                    # the previous corrector step is copied
                    copy = '\cp moi_files/cor.%s.%d.pch moi_files/cor.%s.%d.pch' % (
                        cell, iteration, cell, jteration)
                    print copy
                    os.system(copy)

        sum_power = 0
        sum_fissionpower = 0
        for cell in depletion.mpos[iteration].cell:
            sum_power = depletion.mpos[iteration].ORIGEN_power[cell]
            sum_fissionpower = depletion.mpos[
                iteration].ORIGEN_fissionpower[cell]

        # check for convergence of power
        if power_convergence == 'no':

            if math.fabs(1. - sum_power / depletion.power[iteration]) <= depletion.power_err:
                # the calculated power is within error of the target power
                power_convergence = 'yes'

                # clean workspace
                remove = '\\rm TAPE*'
                print remove
                os.system(remove)

                for i in range(10):
                    depletion.global_decayheat[i] = sum_power - sum_fissionpower

            else:
                # the calculated power is not within error of the target
                power_convergence = 'no'

        else:
            # clean workspace
            print remove
            os.system(remove)


def compBRO(mpo, new_loc=''):
    # This method generates a new MCNP5 input deck based on the ORIGEN PCH
    # files and current MCNP5 input deck

    import os
    import re

    # ensure there are the files you need:
    #	inpiteration
    #   mpo files for the current timestep
    #   moi.cell.iteration.pch files

    iteration = mpo.timestep
    jteration = mpo.timestep + 1
    a = Depletion()
    nct = a.nct

    os.system('\\rm tmp')
    os.system('ls >> tmp')
    main_files = open('tmp').read().split()

    os.system('\\rm tmp')
    os.system('ls moi_files/ >> tmp')
    moi_files = open('tmp').read().split()

    inp_loc = 'inp%d' % iteration
    if inp_loc in main_files:
        # There is an input deck in the main directory
        inp = open(inp_loc).read()
    else:
        # There is not an input deck in the main directory
        print '%s is not in the main directory' % inp_loc
        exit()

    re_XS = re.compile('\d{4,5}[.]\d{2,2}c')

    for cell in mpo.cell:
        new = material()
        pch = material()
        ocf_loc = 'moi_files/moi.%s.%d.pch' % (cell, iteration)
        pch.import_ocf(ocf_loc)

        # normalize to atoms/bn-cm
        pch = pch*(.602214078/mpo.volume[cell])

        # get library and temperature

        library = mpo.mcnp_nuclides.values()[0].replace('.', ' ').split()[1]

        if mpo.transport_module[0] == 'M' or mpo.transport_module[0] == 'm':

            token = '\nm%s ' % mpo.mat_ID[cell]
            k = inp.find(token)

            # isolate the material string by finding the begining of the next
            # card type

            cards = ['\nm', '\nM', '\nf', '\nF', '\nk', '\nK', '\ns', '\nS']

            l = len(inp) - k

            for card in cards:
                m = inp[k+2:].find(card)+2
                if m > 1 and m < l:
                    # there is a card of this type in the system and it is cloer
                    # than the previous type
                    l = m
            mat_card = inp[k+1:k+l]

            nuclides = re_XS.findall(mat_card)

            # filter out trace isotopics
            for nucl in nuclides:
                if nct[nucl] in pch.comp.keys() and nucl in mpo.mcnp_nuclides.values():
                    # this isotopes concentration will be taken from the OCF
                    # file
                    new.comp[nct[nucl]] = pch.comp[nct[nucl]]
                else:
                    # this isotopes concentration will be taken from the
                    # previous input deck
                    new.comp[nct[nucl]] = mpo.concentration[cell][nct[nucl]]

            print mpo.temperature.keys()
            print dir(mpo)
            if 'temperature' in dir(mpo):
                new.mcf(mpo.mat_ID[cell], tmp=mpo.temperature[cell])
            else:
                new.mcf(mpo.mat_ID[cell], xs=library)

            inp = inp.replace(mat_card, new.mat_card)

            # update the density

            token = '\n%s ' % cell
            k = inp.find(token)
            # confirm this is the correct location by checking the material
            # number

            mat_num = inp[k:].split()[1]
            if mat_num == mpo.mat_ID[cell]:
                # this is the correct location!
                l = inp[k:].find(inp[k:].split()[2]) + len(inp[k:].split()[2])
                token = inp[k:k+l]
                new_density = '\n%s %s %1.5E' % (
                    cell, mpo.mat_ID[cell], new.moles())
                inp = inp.replace(token, new_density)
            else:
                print 'could not find cell material density combination - proceed with previous estimate of density'

        elif mpo.transport_module[0] == 'S' or mpo.transport_module[0] == 's':

            token = '\nmat %s' % cell
            k = inp.find(token)

            # isolate the material string by finding the begining of the next
            # card type

            cards = [
                'cell',
                'CELL',
                'Cell',
                'Dep',
                'dep',
                'DEP',
                'DET',
                'det',
                'Det',
                'Disp',
                'DISP',
                'disp',
                'ENE',
                'ene',
                'Ene',
                'Inc',
                'INC',
                'inc',
                'Lat',
                'LAT',
                'lat',
                'mat',
                'Mat',
                'MAT',
                'Mesh',
                'MESH',
                'mesh',
                'NEST',
                'nest',
                'Nest',
                'Par',
                'par',
                'PAR',
                'PBED',
                'pbed',
                'Pbed',
                'Pin',
                'pin',
                'PIN',
                'PLOT',
                'plot',
                'Plot',
                'SET',
                'Set',
                'set',
                'src',
                'SRC',
                'Src',
                'Surf',
                'surf',
                'SURF',
                'Therm',
                'THERM',
                'therm',
                'Trans',
                'TRANS',
                'trans']

            l = len(inp) - k

            for card in cards:
                m = inp[k+3:].find(card)+3
                if m > 2 and m < l:
                    # there is a card of this type in the system and it is cloer
                    # than the previous type
                    l = m
            mat_card = inp[k:k+l-1]

            nuclides = re_XS.findall(mat_card)

            # filter out trace isotopics
            for nucl in nuclides:
                if nct[nucl] in pch.comp.keys() and nucl in mpo.mcnp_nuclides.values():
                    # this isotopes concentration will be taken from the OCF
                    # file
                    new.comp[nct[nucl]] = pch.comp[nct[nucl]]
                else:
                    # this isotopes concentration will be taken from the
                    # previous input deck
                    new.comp[nct[nucl]] = mpo.concentration[cell][nct[nucl]]

            k = mat_card.find('tmp')
            if k != -1:
                temperature = float(mat_card[k:].split()[1])
                new.scf(
                    cell,
                    volume=mpo.volume[cell],
                    temperature=temperature,
                    density='atom/bn-cm',
                    burn='1',
                    library=library)
            else:
                new.scf(
                    cell,
                    volume=mpo.volume[cell],
                    density='atom/bn-cm',
                    burn='1',
                    library=library)

            inp = inp.replace(mat_card, new.mat_card)

    # write new input file
    if new_loc == '':
        new_loc = 'inp%d' % (iteration + 1)
        print 'vi %s' % new_loc
    open(new_loc, 'w').write(inp)


def analysis(depletion):

    # This method prints the mocup depletion analysis output files
    #
    # formats are:
    #	text
    #	csv
    #	python
    #	matlab

    art = 'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMI7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMM$.MMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMI MMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMM.MMMMMMMM=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMM~.MMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMM.,MMMM 7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMZOM   ..MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMM .~MMMMMMN  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMM.MMMMMMMMMMMMM..MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMM MMMMMMMMMMMMMMMMM  OMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMM.MMMMMMMMMMMMMMMMMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMN ?MMMMMMMMMMMMMM  8MMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMM.,MMMMM.MMMMMMMMMMN..MMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNMMMMMMMMMMMM\nMMMMMM..NMMMMMN,M MMMMMMMMMM  MMMM MMM...O= : 8. MMM, MM..MMM. MM .MM MMMM.MM.. MN..MMMMMMM .~M..NM .MMMN MMMM\nMMM  NMMMMMMM .NMM MMMMMMMMMM MMMM.MMM..MMM..MMM..M. MMMM .M  MMMM,MM.MMMM.MM. MMMM  MMMMMM MMMM, MN.MMM.=MMMM\nMMN=MMMMM...MMMMMM.~MMMMMMMMM DMMM.MMM..MMM. MMM..M .MMMM  M..MMMMMMM.MMMM MM. MMMM. MMMMMM MMMMM.MM? MM MMMMM\nMMM. ..$MMMMMMMMMMN MMMMMMMMM..MMM MMM..MMM. MMM..M  MMMM..M..MMMM MM.MMMM.MM. MMMM.IMMMMMM MMMM. MMM M.DMMMMM\nMMMMMMMMMMMMMMMMMN. ,MMMMMMMMM.MMM.MMM..MMM. MMM..MM..MM..MMM..MI $MM. 8 . MM.  ND .MMM.DMM .+M  MMMM$. MMMMMM\nMMMMMMMMMMMM$...MMMMMMMMMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM..MMMMMMMMMMMM MMMMMMMMMM MMMMMMM\nMMMMMMMMM  MMMMMMMMMMMMMMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM..MMMMMMMMMMMM MMMMMMMNM ?MMMMMMM\nMMMMMMMMM.MMMMMMMMMMMMMMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMM= MMMMMMMM8 .MMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMM.MMMM7 .MMMM MMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMZ?MMMM MMMMM.NMMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMM MMMM,MMMMMM MMMMMM= MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMM.MMMMM MMMMMN.MMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMM MMMM.MMMMMM MMMMMM=DMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMM NMMMM,MMMMMM.MMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMM.MMMM MMMMMM$~MMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMM M   MMMMMMMM.MMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMMMMMMMMMMMM MMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMMMMMMMMMMMM 7MMMO MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMMMMMMMMMMMMM.M..MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMMMMMMMMMMMMMM.MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\nMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM'

    # initial output text files
    time_out_cum = []
    time_out_cum2 = []
    k_vector = []

    # Determine the evoluotion of keff as a function of time
    time_sum = 0

    time = []
    time2 = []
    if depletion.method[0] == 'b' or depletion.method[0] == 'B':
        time[:] = depletion.time
        time2[:] = depletion.time
    elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
        time[:] = depletion.delta_time
        time2[:] = depletion.delta_time
        del time2[0]
        del time2[-1]

    for iteration in range(len(time)):

        time_sum += time[iteration]

        time_out_cum.append(float(time_sum))

        jteration = iteration + 1

        # grab k-effective

        if depletion.transport_module[0] == 'm' or depletion.transport_module[0] == 'M':

            # import mcnp output file
            outp_loc = 'outp%d' % jteration
            OUTP = open(outp_loc).read()

            # grab keffective from mcnp output file
            k = OUTP.find('keff = ')
            k_vector.append(float(OUTP[k:].split()[2].replace(')', '')))

        elif depletion.transport_module[0] == 's' or depletion.transport_module[0] == 'S':

            # import serpent results file
            outp_loc = 'inp%d_res.m' % jteration
            OUTP = open(outp_loc).read().split()

            # grad keffective from serpent output file
            k = OUTP.index('IMP_KEFF')
            k_vector.append(float(OUTP[k+6]))

    if depletion.method[0] == 'p' or depletion.method[0] == 'P':
        time_out_cum2[:] = time_out_cum[1:-1]

    elif depletion.method[0] == 'b' or depletion.method[0] == 'B':
        time_out_cum2[:] = time_out_cum[:-1]

    flux_matrix = {}
    power_matrix = {}

    # write cell specific output files
    for cell in depletion.mpos[0].cell:

        # initial flux and power output strings
        flux_matrix[cell] = {}
        power_matrix[cell] = {}

        # initialized concentration evolution
        comp_evo = []

        # grab the composition from ORIGEN PCH files
        # grab the flux and power from the ORIGEN output files
        for iteration in range(len(time2)):
            jteration = iteration + 1

            # import discharge fuel vector timestep

            if depletion.method[0] == 'b' or depletion.method[0] == 'B':
                # beginning mode use beginning of timestep
                ocf_loc = 'moi_files/moi.%s.%d.pch' % (cell, iteration)
            elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
                # predictor corrector mode use midpoint composition
                ocf_loc = 'moi_files/moi.%s.%d.pch' % (cell, jteration)
            mat = material()
            mat.import_ocf(ocf_loc)
            comp_evo.append(mat)

            # import ORIGEN output file for timestep
            origen_loc = 'moi_files/moi.%s.%d.out' % (cell, jteration)
            ORIGEN = open(origen_loc).read()

            # grab the neutron flux from the ORIGEN output file
            k = ORIGEN.find('NEUT. FLUX')
            flux_vector = ORIGEN[k:].split()[3:13]
            flux = 0
            for f in flux_vector:
                flux += float(f)
            flux_matrix[cell][jteration] = flux/len(flux_vector)

            # grab the power from the ORIGEN output file
            power_matrix[cell][jteration] = depletion.mpos[
                iteration].ORIGEN_power[cell]

    # generation concentration matrix [time],[cell], [nuclide]

    if depletion.method[0] in ['b', 'B']:
        I = 0
    elif depletion.method[0] in ['p', 'P']:
        I = 1
    concentration = {}
    for i in range(I, len(time_out_cum)-1):
        concentration[time_out_cum[i]] = depletion.mpos[i].concentration

    depletion.output = output(
        depletion.mpos[0].cell,
        concentration,
        flux_matrix,
        k_vector,
        power_matrix,
        time_out_cum[
            I:-1],
        time_out_cum)

    if depletion.format[0] == 'p' or depletion.format[0] == 'P' or depletion.format[0] == 'a' or depletion.format[0] == 'A':
        # This is formated as a python file
        out_loc = 'mocup_out.py'

        out = "#!/usr/bin/env python"

        out += "\n# depletion cells"
        out += "\ncells = ["
        for cell in depletion.mpos[0].cell:
            out += "'%s'," % cell
        out += ']'

        # print time vector for keff evolution
        out += "\n# time (effective full power days)"
        out += "\ntime_keff = ["
        for t in time_out_cum:
            out += "%1.5E," % t
        out += ']'

        # print keff evolution
        out += "\n# keff evolution"
        out += "\nkeff = ["
        for k in k_vector:
            out += "%1.5E," % k
        out += ']'

        out += '\n# cell dependent values'

        if depletion.method[0] == 'b' or depletion.method[0] == 'B':
            _time_ = 'beginning'
            I = 0
        elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
            _time_ = 'middle'
            I = 1
        out += '\n# time (effective full power days) at %s of timestep' % _time_
        out += '\ntime_cell = ['
        for t in time_out_cum2:
            out += '%1.5E,' % t
        out += ']'

        out += '\nconcentration = {}'

        print time_out_cum
        print time_out_cum2

        for cell in depletion.mpos[0].cell:
            out += "\nconcentration['%s'] = {}" % cell
            for nucl in depletion.mpos[0].concentration[cell].keys():
                out += "\nconcentration['%s']['%s'] = {" % (cell, nucl)
                for j in range(I, len(time_out_cum)-1):
                    # if I = 1, predictor corrector
                    i = j-I
                    out += '%1.5E:%1.5E,' % (time_out_cum[j],
                                             depletion.mpos[i].concentration[cell][nucl])
                out += '}'

        # print power
        out += '\n# power vector (MW)\npower = {}'
        for cell in depletion.mpos[0].cell:
            out += "\npower['%s'] = {" % cell
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E:%1.5E,' % (time_out_cum[j-1+I],
                                         power_matrix[cell][j])
            out += '}'

        # print flux
        out += '\n# flux vector (n/cm2s)\nflux = {}'
        for cell in depletion.mpos[0].cell:
            out += "\nflux['%s'] = {" % cell
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E:%1.5E,' % (time_out_cum[j-1+I],
                                         flux_matrix[cell][j])
            out += '}'

        out += "\n\nprint '%s'" % art.replace('\n', '\\n')

        open('output.py', 'w').write(out)

    if depletion.format[0] == 'c' or depletion.format[0] == 'C' or depletion.format[0] == 'a' or depletion.format[0] == 'A':
        # This is formated as a csv

        out_loc = 'output.csv'

        # print depletion cells
        out = art.replace(',', ".") + '\n\ndepletion cells:\n'

        for cell in depletion.mpos[0].cell:
            out += '%s,' % cell

        # print time vector for keff evolution
        out += '\n\ntime (effective full power days):\n'
        for t in time_out_cum:
            out += '%1.5E,' % t

        # print keff evolution
        out += '\n\nkeff:\n'
        for k in k_vector:
            out += '%1.5E,' % k

        if depletion.method[0] == 'b' or depletion.method[0] == 'B':
            _time_ = 'beginning'
            I = 0
        elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
            _time_ = 'middle'
            I = 1

        time_csv = '\n\nnuclide,time (effective full power days) at %s of timestep\n\t,' % _time_
        for t in time_out_cum2:
            time_csv += '%1.5E,' % t

        out += '\n\n'+'-'*80 + '\ncell dependent parameters\n' + '-'*80 + '\n'

        for cell in depletion.mpos[0].cell:
            out += '\ncell: %s' % cell
            out += '\nconcentration evolution at/bn-cm'
            out += time_csv
            for nucl in depletion.mpos[0].concentration[cell].keys():
                out += '\n%s,' % nucl
                for j in range(I, len(depletion.mpos)):
                    out += '%1.5E,' % depletion.mpos[j].concentration[cell][nucl]

            # print power
            out += '\n\npower evolution (MW)\n\t,'
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E,' % power_matrix[cell][j]

            # print flux
            out += '\n\nflux evolution (n/cm2s)\n\t,'
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E,' % flux_matrix[cell][j]

        open(out_loc, 'w').write(out)

    if depletion.format[0] == 'm' or depletion.format[0] == 'M' or depletion.format[0] == 'a' or depletion.format[0] == 'A':
        # This is formated as a m-file

        out_loc = 'output.m'

        # print time vector for keff evolution
        out = '\n\n% time (effective full power days):\ntime_keff = ['
        for t in time_out_cum:
            out += '%1.5E ' % t
        out += '];'

        # print keff evolution
        out += '\n\n% keff evolution:\nkeff = ['
        for k in k_vector:
            out += '%1.5E ' % k
        out += '];'

        # print times for depletion parameters
        if depletion.method[0] == 'b' or depletion.method[0] == 'B':
            _time_ = 'beginning'
            I = 0
        elif depletion.method[0] == 'p' or depletion.method[0] == 'P':
            _time_ = 'middle'
            I = 1

        time_csv = '\n\n%% time (effective full power days) at %s of timestep\ntime_dep = [' % _time_
        for t in time_out_cum2:
            time_csv += '%1.5E ' % t
        time_csv += '];'

        # print nuclide identifiers

        out += '\n% nuclide identifiers' + '\nzai = ['
        id = '\n% nuclide IDs '

        i = 1
        for nucl in sorted(depletion.mpos[0].concentration[cell].keys()):
            out += '\n%s' % nucl
            id += '\ni%s = %d' % (nucl, i)
            i += 1

        out += '\v]'

        out += id

        out += '\n\n% '+'-'*80 + \
            '\n% cell dependent parameters\n%' + '-'*80 + '\n'

        # print concentration matrix

        for cell in depletion.mpos[0].cell:
            out += '\n%% concentration vector for cell: %s (at/bn-cm)' % cell
            out += time_csv
            out += '\nconcentration_%s = [' % cell
            for nucl in sorted(depletion.mpos[0].concentration[cell].keys()):
                out += '\n'
                for j in range(1, len(depletion.mpos)):
                    out += '%1.5E ' % depletion.mpos[j].concentration[cell][nucl]
                out += '; %% %s' % nucl
            out += '\n];'
            # power vectors
            out += '\n%% power vector (MW)\npower_%s = [' % cell
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E ' % power_matrix[cell][j]
            out += '];'

            # flux vectors
            out += '\n%% flux vector (n/cm2s)\n flux_%s = [' % cell
            for j in sorted(power_matrix[cell].keys()):
                out += '%1.5E ' % flux_matrix[cell][j]
            out += '];'

        out += "\ndisp('" + art[1:].replace('\n', "')\ndisp('") + "')"

        open(out_loc, 'w').write(out)
    else:
        # Print a standard text file
        pass


class output:

    def __init__(
            self,
            cells,
            concentration,
            flux,
            keff,
            power,
            time_cell,
            time_keff):

        self.cells = cells
        self.concentration = concentration
        self.flux = flux
        self.keff = keff
        self.power = power
        self.time_cell = time_cell
        self.time_keff = time_keff
