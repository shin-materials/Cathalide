# Cathalide

Cathalide intends to modulate oragnic molecules in perovskite structures. Unlike inorganic materials, organic molecules are often considered as coarse grain. But conventional DFT codes cannot group atoms to form a group, which makes the operation of translation and rotation difficult.

This uses modified version of pymatgen functions. 
Currently, only VASP formats are allows (POSCAR) files, but Quantum Espresso will be soon to be implemented.

Features to be implemented:
1. Assigning reference points more effectively. For example, if I have xyz file of molecule, it is hard to align to certain direction in the cell.
2. Modifying the molecular structures. For C-C chains, an user might want to rotate the part of molecules.
3. Prepare high-symmetry molecules. For example, MA (H3C-NH3) can be prepared with three-fold rotation symmetry.
4. Preparation of perovskites matrix with different distortion patterns.

