#    This file is part of 3D TMH Complexity.
#
#    3D TMH Complexity is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    3D TMH Complexity is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with 3D TMH Complexity.  If not, see <http://www.gnu.org/licenses/>.


from pymol import cmd, stored, math
import sys

def complexb(mol, startaa=0, visual="Y"):
    """
    Replaces B-factors with a list of values contained in a plain txt file

    usage: complexb mol, [startaa, [source, [visual]]]

    mol = any object selection (within one single object though)
    startaa = number of first amino acid in 'new B-factors' file (default=1)
    source = name of the file containing new B-factor values (default=newBfactors.txt)
    visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)

    example: complexb 1LVM and chain A
    """
    for subchains in cmd.get_chains(mol):
        source="bfactors_%s_%s.txt" % (mol, subchains)
        obj = cmd.get_object_list(mol)[0]
        cmd.alter(mol, "b=-1.0")
        infile = open(source, 'r')
        counter = int(startaa)
        bfacts = []
        for line in infile.readlines():
            bfact = float(line)
            bfacts.append(bfact)
            cmd.alter("%s and resi %s and n. CA" % (mol, counter), "b=%s" % bfact)
            counter = counter + 1
        if visual == "Y":
            cmd.show_as("cartoon", mol)
            cmd.cartoon("putty", mol)
            cmd.spectrum("b", "rainbow", "%s and n. CA " % mol,0,4)
            cmd.ramp_new("count", obj, [0, 4], "rainbow")
            cmd.recolor()

# This is needed to load the script in pymol
cmd.extend("complexb", complexb)
