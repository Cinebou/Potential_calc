# Instructions
You can calculate the Coulombic and VDW potential of each specific molecule in the framwork separately.
Import your force field ('pseudo_atoms.def', 'force_fieldmixing_rules.def') to the sub derectory.
Then, prepare the CO2 pdb file and framework pdb file as shown below.

Finally, you can run the calucaltion.
'''
python calc_potential.py
'''

# How to get framework pdb file with atom numbering
1. Install Pymol
2. Open cif file which has atom numbering
3. From 'File', open 'export moleculaes'. 
4. In 'Generic Options', click 'original atom order'.
5. In 'PDB options', click 'Retain atom ids'.
6. Save as pdb files.


### Attention
1. In the CO2 pdb files, remove the header. Only atom line can be accepted.
2. Electric potential is not calculated based on Ewald method.
3. If you want to run this calculation with different gas molecule, you need to change the function '__modify_name()'
