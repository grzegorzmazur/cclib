# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, Grzegorz Mazur
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for niedoida output files"""

from __future__ import print_function
import re

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils


class Niedoida(logfileparser.Logfile):
    """A niedoida log file."""

    def __init__(self, *args, **kwargs):
        super(Niedoida, self).__init__(logname="niedoida", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "niedoida log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'niedoida("%s")' % (self.filename)

    def before_parsing(self):
        self.et_symmetry_prefix = ""

    def after_parsing(self):
        if hasattr(self, "etenergies"):
            self.etenergies, self.etoscs, self.etsyms, self.etsecs = (list(x) for x in zip(
                *sorted(zip(self.etenergies, self.etoscs, self.etsyms, self.etsecs), key=lambda t: t[0])))

    def extract_vector(self, inputfile, size, elems_per_section):
        no_sections = size // elems_per_section
        if size % elems_per_section != 0:
            no_sections = no_sections + 1

        labels = []
        elems = []

        for i in range(0, no_sections):
            labels = labels + inputfile.next().split()
            elems = elems + [float(x) for x in inputfile.next().split()]

        return (labels, elems)

    def extract_matrix_section(self, inputfile, no_rows, no_cols):
        line = inputfile.next().strip()
        if len(line) == 0:
            line = inputfile.next().strip()
        col_labels = line.split()
        row_labels = []
        elems = []
        for i in range(0, no_rows):
            line = inputfile.next().split()
            row_labels.append(line[0])
            elems.append([float(x) for x in line[1:]])

        return (row_labels, col_labels, elems)

    def extract_matrix(self, inputfile, size, cols_per_section):
        no_sections = size[1] // cols_per_section
        if size[1] % cols_per_section != 0:
            no_sections = no_sections + 1

        row_labels = []
        col_labels = []
        elems = size[0] * [[]]

        for i in range(0, no_sections):
            section_row_labels, section_col_labels, section_elems = self.extract_matrix_section(
                inputfile, size[0], min(cols_per_section, size[0] - i * cols_per_section))
            if i == 0:
                row_labels = section_row_labels
            col_labels = col_labels + section_col_labels
            elems = [x[0] + x[1] for x in zip(elems, section_elems)]

        return (row_labels, col_labels, elems)

    def extract_atomic_basis(self, inputfile, no_contractions):
        b = []
        for c in range(0, no_contractions):
            type, no_primitives = inputfile.next().split()
            type = type.upper()
            no_primitives = int(no_primitives)
            primitives = []
            for p in range(0, no_primitives):
                primitives.append(tuple(map(float, inputfile.next().split())))
            b.append((type, primitives))
        return b

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        if line[0:7] == "charge:":
            self.charge = int(line.split(":")[1])

        if line[0:13] == "multiplicity:":
            self.mult = int(line.split(":")[1])

#        if line[0:15] == "begin basis set":
#            line = inputfile.next()
#            self.gbasis = []
#            while line.strip() != "end basis set":
#                no_contractions = int(line.split()[1])
#                self.gbasis.append(self.extract_atomic_basis(inputfile, no_contractions))
#                line = inputfile.next()

        if line[0:12] == "coordinates:":
            size, cols_per_line = line[12:].split(",")
            size = list(map(int, size.split("x")))
            cols_per_line = int(cols_per_line)

            atoms, dummy, atomcoords = self.extract_matrix(
                inputfile, size, cols_per_line)

            atomcoords = [
                [utils.convertor(q, "bohr", "Angstrom") for q in xyz] for xyz in atomcoords]

            if not hasattr(self, "atomcoords"):
                self.atomcoords = []

            self.atomcoords.append(atomcoords)

            if not hasattr(self, "atomnos"):
                pt = utils.PeriodicTable()
                self.atomnos = [pt.number[a.capitalize()] for a in atoms]

            self.natom = len(atoms)

        if line[0:13] == "total_energy:":

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []

            self.scfenergies.append(utils.convertor(
                float(line.split(":")[1]), "hartree", "eV"))

        if line[0:23] == "alpha_orbital_energies:":
            size, elems_per_line = list(map(int, line[23:].split(",")))
            labels, energies = self.extract_vector(
                inputfile, size, elems_per_line)
            self.moenergies = [energies]

        if line[0:22] == "beta_orbital_energies:":
            size, elems_per_line = list(map(int, line[22:].split(",")))
            labels, energies = self.extract_vector(
                inputfile, size, elems_per_line)
            self.moenergies.append(energies)

        if line[0:26] == "alpha_orbital_occupations:":
            size, elems_per_line = list(map(int, line[26:].split(",")))
            self.homos = [sum(self.extract_vector(
                inputfile, size, elems_per_line)[1]) - 1]

        if line[0:25] == "beta_orbital_occupations:":
            size, elems_per_line = list(map(int, line[25:].split(",")))
            self.homos.append(sum(self.extract_vector(
                inputfile, size, elems_per_line)[1]) - 1)

        if line[0:22] == "alpha_mo_coefficients:":
            size, cols_per_line = line[22:].split(",")
            size = list(map(int, size.split("x")))
            cols_per_line = int(cols_per_line)

            self.nbasis = size[0]
            self.nmo = size[1]

            row_labels, col_labels, coeffs = self.extract_matrix(
                inputfile, size, cols_per_line)

            self.mocoeffs = [coeffs]

        if line[0:21] == "beta_mo_coefficients:":
            size, cols_per_line = line[21:].split(",")
            size = list(map(int, size.split("x")))
            cols_per_line = int(cols_per_line)

            row_labels, col_labels, coeffs = self.extract_matrix(
                inputfile, size, cols_per_line)

            self.mocoeffs.append(coeffs)

        if line[0:20] == "begin singlet states":
            self.et_symmetry_prefix = "Singlet"

        if line[0:20] == "begin triplet states":
            self.et_symmetry_prefix = "Triplet"

        if line[0:9] == "energies:":
            size, elems_per_line = list(map(int, line[9:].split(",")))
            labels, energies = self.extract_vector(
                inputfile, size, elems_per_line)

            energies = [utils.convertor(x, "hartree", "cm-1")
                        for x in energies]

            if not hasattr(self, "etenergies"):
                self.etenergies = []
                self.etoscs = []
                self.etsyms = []
                self.etsecs = []

            self.etenergies = self.etenergies + energies
            self.etsyms = self.etsyms + \
                len(energies) * [self.et_symmetry_prefix]
            if self.et_symmetry_prefix == "Triplet":
                self.etoscs = self.etoscs + len(energies) * [0.0]

        if line[0:13] == "coefficients:":
            size, elems_per_line = list(map(int, line[13:].split(",")))
            labels, coeffs = self.extract_vector(
                inputfile, size, elems_per_line)
            confs = [map(int, l.split("->")) for l in labels]
            o, v = zip(*confs)
            self.etsecs = self.etsecs + \
                [[(x[0], 0), (x[1], 0), x[2]] for x in zip(o, v, coeffs)]

        if line[0:20] == "oscillator_strength:":
            size, elems_per_line = list(map(int, line[20:].split(",")))
            labels, f = self.extract_vector(inputfile, size, elems_per_line)
            self.etoscs = self.etoscs + list(map(float, f))

        if line[0:17] == "total_mp2_energy:":

            if not hasattr(self, "mpenergies"):
                self.mpenergies = []

            self.mpenergies.append([])

            mp2energy = self.float(line.split(":")[1])
            self.mpenergies[-1].append(utils.convertor(mp2energy,
                                                       "hartree", "eV"))


if __name__ == "__main__":
    import doctest
    import niedoidaparser
    import sys

    if len(sys.argv) == 1:
        doctest.testmod(niedoidaparser, verbose=False)

    if len(sys.argv) >= 2:
        parser = niedoidaparser.Niedoida(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))
