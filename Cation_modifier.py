from pymatgen.core import Structure
from pymatgen.core import Element
from glob import glob


struct = Structure.from_file("Test_structures/FPB_bulk_cubic.vasp")

#Need to do:
# Rotation of cations
#   1) Identify molecule
#       How? --> Bond information (maybe from VESTA)

# Read bond data
# Input data from user: what elements composing the molecule
elements_list=['C','H','N']



# When replacing inorganic cation to organic molecule,
# I need to operate with the center of mass
# command: Element.C.atomic_mass