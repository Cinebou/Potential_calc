import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np

"""
LJ potential calculation for CO2 adsorption
if you want to do this calculation for other molecules like H2O, you should change the '__modify_name(self)'
"""
class Potential():
    def __init__(self):
        self.avogadro = 6.02214e23  # /mol
        self.R = 8.31446261815324 / 1000  # kJ/K/mol
        self.k_charge = 8.987552e9

        self.num_atom_in_molecule = 3

        self.CO2_file = './force_field/CO2_example.pdb'
        self.frame_file = './force_field/MIL-101_all.pdb'
        self.__import_position()
        self.__import_potential()

        self.cutoff = 12 # A
        self.cutoff_2 = self.cutoff**2
        self.mol_amount = self.num_CO2 / self.avogadro



    # read the pdb file
    def __import_position(self):
        f_cols=['A','number','atom','M','x','y','z','b','c','atom spec']
        self.pos_CO2 = pd.read_csv(self.CO2_file, delim_whitespace = True, skiprows=2, header=None,names=f_cols)
        self.num_CO2 = len(self.pos_CO2) / self.num_atom_in_molecule
        self.pos_frame = pd.read_csv(self.frame_file, delim_whitespace = True, skiprows=1, header=None,names=f_cols)[:-1]


    # read the parameters for potential functions
    def __import_potential(self):
        col_names = ['atom_type', 'function type', 'epsilon', 'sigma']
        self.LJ_pot = pd.read_csv('./force_field/force_field_mixing_rules.def', skiprows=7,  delim_whitespace = True, names=col_names)
        col_names = ["atom_type","print","as","chem","oxidation" ,"mass","charge","polarization","B-factor"," radii","connectivity","anisotropic","anisotropic-type","tinker-type"]
        self.charges = pd.read_csv('./force_field/pseudo_atoms.def', skiprows=3,  delim_whitespace = True, names=col_names)
        self.__modify_name()


    # the name of atom in pdb file and force field file are differnt, modify them
    def __modify_name(self):
        self.LJ_pot = self.LJ_pot.replace('O_co2','O')
        self.LJ_pot = self.LJ_pot.replace('C_co2','C')
        self.charges = self.charges.replace('O_co2','O')
        self.charges = self.charges.replace('C_co2','C')



    # take out LJ params for each atom
    def params_LJ(self, atom_name):
        epsilon = float(self.LJ_pot[self.LJ_pot['atom_type']==atom_name]['epsilon']) # ε / kb
        sigma = float(self.LJ_pot[self.LJ_pot['atom_type']==atom_name]['sigma'])
        return epsilon, sigma


    # take out the charges params for each atom
    def params_charge(self, atom_name):
        charge = float(self.charges[self.charges['atom_type']==atom_name]['charge'])
        return charge


    # distance between molecules powers to 2
    def distance_2(self, molA, molB):
        dx = molA['x'] - molB['x']
        dy = molA['y'] - molB['y']
        dz = molA['z'] - molB['z']
        if (dx <= self.cutoff) and (dy <= self.cutoff) and (dz <= self.cutoff):
            return dx**2 + dy**2 + dz**2
        else:
            return 1000000
            # just to calculate fast, 

    
    # Lorentz-Berthelot mixing rules
    def lorentz_berthelot(self, atom_a, atom_b):
        ea, sa = self.params_LJ(atom_a['atom'])
        eb, sb = self.params_LJ(atom_b['atom'])
        e_LB = sqrt(ea*eb)
        s_LB = (sa+sb)/2
        return e_LB, s_LB

    
    # calculate LJ potential function from distance
    def LJ_function(self, r2, epsilon, sigma):
        if (r2 <= self.cutoff_2):
            r_s_6 = (sigma**2/r2)**3
            pot = 4 * epsilon * (r_s_6**2 - r_s_6)
            return pot
        # when the atom pairs are too far away, retrun potential zero
        else:
            return 0


    # calculate LJ potential function from distance
    def charge_function(self, r2, qA, qB):
        if (r2 <= self.cutoff_2):
            r = sqrt(r2)
            pot = -self.k_charge * qA * qB / r
            return pot
        # when the atom pairs are too far away, retrun potential zero
        else:
            return 0
        

    # LJ potential calculation for each pairs
    def LJ_each(self, molA, molB):
        eps, sig = self.lorentz_berthelot(molA, molB)
        r2_ = self.distance_2(molA, molB)
        pot_each = self.LJ_function(r2_, eps, sig)
        return pot_each


    # LJ potential loop for all molecules pairs of frame and CO2
    def LJ_CO2_frame(self):
        LJ_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2, rowCO2 in self.pos_CO2.iterrows():
            CO2_atom = rowCO2
            # for loop (frameworks)
            for indexFrame, rowFrame in self.pos_frame.iterrows():
                frame_atom = rowFrame
                LJ_pot_all += self.LJ_each(CO2_atom, frame_atom)
        return LJ_pot_all / self.R



    # LJ potential loop for all molecules pairs of CO2 and CO2
    def LJ_CO2_CO2(self):
        LJ_pot_all = 0
        # for loop (CO2 molecules)
        for indexCO2_a, rowCO2_a in self.pos_CO2.iterrows():
            CO2_atom_a = rowCO2_a
            # for loop (CO2 molecules)
            for indexCO2_b, rowCO2_b in self.pos_CO2.iterrows():
                # atoms in the same molecule should be excluded
                if (rowCO2_a['number']-1)//self.num_atom_in_molecule==(rowCO2_b['number']-1)//self.num_atom_in_molecule:
                    continue
                else:
                    CO2_atom_b = rowCO2_b
                    LJ_pot_all += self.LJ_each(CO2_atom_a, CO2_atom_b)
        return LJ_pot_all / self.R


    # show the summation of potential
    def potential(self):
        LJ_CO2 = self.LJ_CO2_CO2()
        LJ_frame = self.LJ_CO2_frame()
        print('LJ potential between framework and CO2 molecules  :   ', LJ_CO2  , 'kJ/mol')
        print('LJ potential between CO2 and CO2molecules         :   ', LJ_frame, 'kJ/mol')
        print('LJ potential overall                              :   ', LJ_CO2 + LJ_frame, 'kJ/mol')


    # LJ potential curve
    def LJ_graph(self):
        atom_A = 'Cr1'
        atom_B = 'F1'
        ea, sa = self.params_LJ(atom_A)
        eb, sb = self.params_LJ(atom_B)
        e_LB = sqrt(ea*eb)
        s_LB = (sa+sb)/2

        r_list = np.linspace(s_LB-0.2, s_LB*2, 100)
        pot_list = [self.LJ_function(r**2, e_LB, s_LB)*self.R for r in r_list]

        fig_pot = plt.figure()
        ax_pot = fig_pot.add_subplot(111)
        ax_pot.plot(r_list, pot_list, label=atom_A+'-'+atom_B)
        ax_pot.set_xlabel('distance  A')
        ax_pot.set_ylabel('potential  kJ/mol')
        ax_pot.legend()

        plt.show()

def main():
    mil101 = Potential()
    #mil101.LJ_graph()
    mil101.potential()
    
if __name__=='__main__':
    main()