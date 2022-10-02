import pandas as pd
from math import sqrt, ceil
import matplotlib.pyplot as plt
import numpy as np
from calc_pot_total import Potential
import log_output


def main():
    CO2_file = './force_field/VMD/298.15_100000.pdb'
    frame_file = './force_field/MIL-101_all.pdb'
    Results_file_name = './Results/example.csv'
    # specific atom number in pdb file
    # it should be the whole molecule, which is the set of three atoms

    # 1 bar, 298.15 K, MIL-101
    #c_atoms = [176]  # 1 molecules in ST
    #c_atoms = [182]  # 1 molecules in ST
    #c_atoms = [269, 227]  # 2 molecules in ST
    #c_atoms = [275, 131]  # 2 molecules in ST
    c_atoms = [104, 140, 137]  # 3 molecules in ST

    specified_atoms = C_to_O(c_atoms)
    
    #specified_atoms = [7,8,9]
    pot = Potential_each(CO2_file, frame_file, specified_atoms, Results_file_name)
    pot.potential_each()
    return 0
    

def C_to_O(C_list):
    co2_list = []
    for c in C_list:
        co2_list += [c-1, c, c+1]
    return co2_list


class Potential_each(Potential):
    def __init__(self,CO2_file_, frame_file_, specified_atom_, Results_file_name_):
        super().__init__(CO2_file_, frame_file_,Results_file_name_)
        self.co2_pickup(specified_atom_)

    # for calc_each_pot
    def co2_pickup(self, atom_numbers):
        # devide the CO2 position file into two parts, one is used as center atom, the other is used as a frame
        self.not_terget_CO2 = self.pos_CO2.query('number != @atom_numbers')

        self.target_CO2 = pd.concat([self.pos_CO2[self.pos_CO2['number']==r] for r in atom_numbers], axis=0)
        self.new_pos_frame = pd.concat([self.pos_frame, self.not_terget_CO2], axis=0)
        
        # initialize the number of atoms
        self.init_const_each()
        return 0

    # constant as a function of number of CO2 molecule
    def init_const_each(self):
        self.avogadro = 6.02214e23  # /mol
        self.num_CO2 = len(self.target_CO2) / self.num_atom_in_molecule
        self.R = 8.31446261815324 / 1000  /self.num_CO2 # kJ/K/mol
        self.k_charge = 8.987552e9 * (1.60217733e-19) * (1.60217733e-19) * 1e10 * self.avogadro /1000 /self.num_CO2 # kJ/mol
        self.mol_amount = self.num_CO2 / self.avogadro


    # LJ potential loop for molecule pairs of frame and CO2
    def LJ_CO2_frame_each(self):
        LJ_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2, rowCO2 in self.target_CO2.iterrows():
            CO2_atom = rowCO2

            # for loop (frameworks)
            for indexFrame, rowFrame in self.pos_frame.iterrows():
                frame_atom = rowFrame
                r2 = self.distance_2(CO2_atom, frame_atom)
                if r2 <= self.cutoff_2:
                    LJ_pot_all += self.LJ_each(CO2_atom, frame_atom, r2)

                    log_output.log_any_msg("{} : {}".format(self.LJ_each(rowCO2, rowFrame,r2), r2))
        return LJ_pot_all * self.R

    # electric potential loop for all molecules pairs of frame and CO2
    def charge_CO2_frame_each(self):
        charge_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2, rowCO2 in self.target_CO2.iterrows():
            CO2_atom = rowCO2
            # for loop (frameworks)
            for indexFrame, rowFrame in self.pos_frame.iterrows():
                frame_atom = rowFrame
                r2 = self.distance_2(CO2_atom, frame_atom)
                charge_pot_all += self.charge_each(CO2_atom, frame_atom, r2)
        return charge_pot_all * self.k_charge



    # LJ potential loop for all molecules pairs of CO2 and CO2
    def LJ_CO2_CO2_each(self):
        LJ_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2_a, rowCO2_a in self.target_CO2.iterrows():
            CO2_atom_a = rowCO2_a
            # for loop (CO2 molecules)
            for indexCO2_b, rowCO2_b in self.pos_CO2.iterrows():
                # atoms in the same molecule should be excluded
                if (rowCO2_a['number']-1)//self.num_atom_in_molecule==(rowCO2_b['number']-1)//self.num_atom_in_molecule:
                    continue
                else:
                    CO2_atom_b = rowCO2_b
                    r2 = self.distance_2(CO2_atom_a, CO2_atom_b)
                    if r2 <= self.cutoff_2:
                        LJ_pot_all += self.LJ_each(CO2_atom_a, CO2_atom_b, r2)
        return LJ_pot_all * self.R

    # LJ potential loop for all molecules pairs of CO2 and CO2
    def charge_CO2_CO2_each(self):
        charge_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2_a, rowCO2_a in self.target_CO2.iterrows():
            CO2_atom_a = rowCO2_a
            # for loop (CO2 molecules)
            for indexCO2_b, rowCO2_b in self.pos_CO2.iterrows():
                # atoms in the same molecule should be excluded
                if (rowCO2_a['number']-1)//self.num_atom_in_molecule==(rowCO2_b['number']-1)//self.num_atom_in_molecule:
                    continue
                else:
                    CO2_atom_b = rowCO2_b
                    r2 = self.distance_2(CO2_atom_a, CO2_atom_b)
                    charge_pot_all += self.charge_each(CO2_atom_a, CO2_atom_b, r2)
        return charge_pot_all * self.k_charge

    
    # show the summation of potential
    def potential_each(self):
        LJ_CO2 = self.LJ_CO2_CO2_each()
        LJ_frame = self.LJ_CO2_frame_each()
        print('VDW potential between framework and CO2 molecules  :   ', LJ_frame  , 'kJ/mol')
        print('VDW potential between CO2 and CO2 molecules        :   ', LJ_CO2, 'kJ/mol')
        print('VDW potential overall                              :   ', LJ_CO2 + LJ_frame, 'kJ/mol\n')

        charge_CO2 = self.charge_CO2_CO2_each()
        charge_frame = self.charge_CO2_frame_each()
        print('electric pot between framework and CO2 molecules  :   ', charge_frame  , 'kJ/mol')
        print('electric pot between CO2 and CO2 molecules        :   ', charge_CO2, 'kJ/mol')
        print('electric pot overall                              :   ', charge_CO2 + charge_frame, 'kJ/mol\n')

        print('total heat of adsorption                          :   ',LJ_CO2 + LJ_frame + charge_CO2 + charge_frame, 'kJ/mol\n')
        return 0
    
    
if __name__ == '__main__':
    main()