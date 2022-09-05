# Instructions
You can calculate the Coulombic and VDW potential of each specific molecule in the framwork separately.
Import your force field ('pseudo_atoms.def', 'force_fieldmixing_rules.def') to the sub derectory.
Then, prepare the CO2 pdb file and framework pdb file as shown below.

Finally, you can run the calucaltion for all molecules in the system.
'''
python calc_pot_total.py
'''

Or, if you want to calculate for each molecule, you can specify them and calculte the potential for each of them.
'''
python calc_pot_each.py
'''
You can find the molecules you are intereted in from VMD and record their position. From the positions, you can find the specific number of the atom in PDB file. Then run this code with the number.


# How to get framework pdb file with atom numbering
1. Install Pymol
2. Open cif file which has atom numbering
3. From 'File', open 'export moleculaes'. 
4. In 'Generic Options', click 'original atom order'.
5. In 'PDB options', click 'Retain atom ids'.
6. Save as pdb files.


## Attention
1. What is calculated is energy, not enthalpy. 
2. In the CO2 pdb files, remove the header. Only atom line can be accepted.
3. Electric potential is not calculated based on Ewald method.
4. If you want to run this calculation with different gas molecule, you need to change the function '__modify_name()'
5. The periodic boundary codition is not considered. (Later)