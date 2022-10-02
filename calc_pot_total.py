import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from math import radians, cos, sin
from os import path, mkdir
import log_output

cell_L_x = 62.840 
cell_L_y = 62.840 
cell_L_z = 62.840 

cell_angle_a = radians(60)	
cell_angle_b = radians(60)
cell_angle_c = radians(60)

global ax, bx, by, cx, cy, cz, file_name
ax = cell_L_x
bx = cell_L_y*cos(cell_angle_c)
by = cell_L_y*sin(cell_angle_c)
cx = cell_L_z*cos(cell_angle_b)
cy = cell_L_z*(cos(cell_angle_a)-cos(cell_angle_b)*cos(cell_angle_c))/sin(cell_angle_c)
cz = sqrt(cell_L_z**2 - cx**2 - cy**2)

if not path.exists('./Results'):
    mkdir('./Results')

def main():
    temp_pres = '298.15_100000'
    CO2_file = './force_field/VMD/'+temp_pres+'.pdb'
    frame_file = './force_field/MIL-101_all.pdb'
    Results_file_name = './Results/'+temp_pres+'.csv'
    mil101 = Potential(CO2_file, frame_file, Results_file_name)
    #mil101.LJ_graph()
    #mil101.charge_graph()
    mil101.pot_for_all()
    
"""
LJ potential calculation for CO2 adsorption
"""
class Potential(object):
    def __init__(self, CO2_file, frame_file, res_file):
        # to be specified for each case
        self.CO2_file = CO2_file
        self.frame_file = frame_file
        self.res_file = res_file
        self.num_atom_in_molecule = 3
        self.cutoff = 12 # A, only for VDW interaction
        self.surface_depth = self.cutoff
        
        # constant and data reader
        self.cutoff_2 = self.cutoff**2
        self.__import_position()
        self.__import_potential()
        self.init_const()
        print('Finish reading all the positions from files')


    # read the pdb file
    def __import_position(self):
        self.f_cols=['A','number','atom','M','x','y','z','b','c','atom spec']
        self.pos_CO2 = pd.read_csv(self.CO2_file, delim_whitespace = True,header=None,names=self.f_cols)
        self.pos_frame_or = pd.read_csv(self.frame_file, delim_whitespace = True, skiprows=1, header=None,names=self.f_cols)[:-1]
        self.pos_CO2_PBD = self.__repro_pbd(self.pos_CO2)
        self.pos_frame = self.__repro_pbd(self.pos_frame_or)

    # read the parameters for potential functions
    def __import_potential(self):
        col_names = ['atom_type', 'function type', 'epsilon', 'sigma']
        self.LJ_pot = pd.read_csv('./force_field/force_field_mixing_rules.def', skiprows=7,  delim_whitespace = True, names=col_names)
        col_names = ["atom_type","print","as","chem","oxidation" ,"mass","charge","polarization","B-factor"," radii","connectivity","anisotropic","anisotropic-type","tinker-type"]
        self.charges = pd.read_csv('./force_field/pseudo_atoms.def', skiprows=3,  delim_whitespace = True, names=col_names)
        self.__modify_name()


    # the name of atom in pdb file and force_field file should be the same
    def __modify_name(self):
        self.LJ_pot = self.LJ_pot.replace('O_co2','O')
        self.LJ_pot = self.LJ_pot.replace('C_co2','C')
        self.charges = self.charges.replace('O_co2','O')
        self.charges = self.charges.replace('C_co2','C')

    
    # reproduce the position of the atoms into PBC
    def __repro_pbd(self, pos):
        pos_num = len(pos)
        for k in range(pos_num):
            for x in range(3):
                for y in range(3):
                    for z in range(3):
                        if x==1 & y==1 & z==1:
                            continue
                        imT = pos.iloc[k].copy(deep=True)
                        imT['x'] += cx*(z-1) + bx*(y-1) + ax*(x-1)
                        imT['y'] += cy*(z-1) + by*(y-1)
                        imT['z'] += cz*(z-1)

                        # only the surface of the surrounding cell can be listed
                        if self.__judge_surface(imT['x'], imT['y'], imT['z']):
                            pos = pos.append(imT)
        return pos


    # judge the surrounding surface within the distance of cutoff radius
    def __judge_surface(self, x, y, z):
        # surface depth is extended in the virtical direction of the wall 
        x_min_surface = (bx*cz*y-(bx*cy-by*cx)*z)/by/cz - self.surface_depth*10.0
        x_max_surface = x_min_surface + ax + self.surface_depth*10.0
     
        y_min_surface = z*cy/cz - self.surface_depth*10.0
        y_max_surface = y_min_surface + by + self.surface_depth*10.0

        # z surface is pararell to the original xyz coordinates 
        z_min_surface = 0 - self.surface_depth
        z_max_surface = cz + self.surface_depth
        
        if (x_min_surface <= x <= x_max_surface) & (y_min_surface <= y <= y_max_surface) & (z_min_surface <= z <= z_max_surface):
            return True
        else:
            return False
        
    
    # constant as a function of number of CO2 molecule
    def init_const(self):
        self.avogadro = 6.02214e23  # /mol
        self.num_CO2 = len(self.pos_CO2) / self.num_atom_in_molecule
        self.R = -8.31446261815324 / 1000  # kJ/K/mol
        self.k_charge = -8.987552e9 * (1.60217733e-19) * (1.60217733e-19) * 1e10 * self.avogadro /1000 # kJ/mol
        self.mol_amount = self.num_CO2 / self.avogadro



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
        return dx**2 + dy**2 + dz**2
    
    # Lorentz-Berthelot mixing rules
    def lorentz_berthelot(self, atom_a, atom_b):
        ea, sa = self.params_LJ(atom_a['atom'])
        eb, sb = self.params_LJ(atom_b['atom'])
        e_LB = sqrt(ea*eb)
        s_LB = (sa+sb)/2
        return e_LB, s_LB

    
    # calculate LJ potential function from distance
    def LJ_function(self, r2, epsilon, sigma):
        r_s_6 = (sigma**2/r2)**3
        pot = 4 * epsilon * (r_s_6**2 - r_s_6)
        return pot
        

    # calculate LJ potential function from distance
    def charge_function(self, r2, qA, qB):
        pot = qA * qB / sqrt(r2)
        return pot
        
    # LJ potential calculation for each pairs
    def LJ_each(self, molA, molB,r2):
        eps, sig = self.lorentz_berthelot(molA, molB)
        pot_each = self.LJ_function(r2, eps, sig)
        return pot_each

    # electric potential calculation for each pairs
    def charge_each(self, molA, molB,r2):
        qA = self.params_charge(molA['atom'])
        qB = self.params_charge(molB['atom'])
        charge = self.charge_function(r2, qA, qB)
        return charge
        
    # pot for each atoms
    def pot_for_all(self):
        # prepare results file
        res_potential = pd.DataFrame(columns=['atom_num','x','y','z', 'LJ_CO2', 'LJ_frame', 'Elec_CO2', 'Elec_frame', 'Total_pot'])
        res_potential.to_csv(self.res_file) 

        # for loop (CO2 molecules)
        for indexCO2, rowCO2 in self.pos_CO2.iterrows():
            LJ_frame = 0
            Elec_frame = 0

            # frame iteration
            for indexFrame, rowFrame in self.pos_frame.iterrows():
                r2 = self.distance_2(rowCO2, rowFrame)
                Elec_frame += self.charge_each(rowCO2, rowFrame,r2)
                if r2 <= self.cutoff_2:
                    LJ_frame += self.LJ_each(rowCO2, rowFrame,r2)
            
            Elec_CO2 = 0
            LJ_CO2 = 0
            # CO2 iteration
            for indexCO2_b, rowCO2_b in self.pos_CO2_PBD.iterrows():
                r2 = self.distance_2(rowCO2, rowCO2_b)
                # atoms in the same molecule should be excluded
                if (rowCO2['number']-1)//self.num_atom_in_molecule==(rowCO2_b['number']-1)//self.num_atom_in_molecule:
                    continue
                Elec_CO2 += self.charge_each(rowCO2, rowCO2_b,r2)
                if r2 <= self.cutoff_2:
                    LJ_CO2 += self.LJ_each(rowCO2, rowCO2_b,r2)

            # const for potential
            LJ_CO2 = LJ_CO2 * self.R
            LJ_frame = LJ_frame * self.R
            Elec_CO2 = Elec_CO2 * self.k_charge
            Elec_frame = Elec_frame * self.k_charge
            Total_pot = LJ_CO2 + LJ_frame + Elec_CO2 + Elec_frame

            # write_all results
            res_potential = pd.DataFrame([rowCO2['number'], rowCO2['x'], rowCO2['y'], rowCO2['z'], LJ_CO2, LJ_frame, Elec_CO2, Elec_frame, Total_pot]).T
            res_potential.to_csv(self.res_file, mode='a', header=False)

            print('calculation going...  : ',int(rowCO2['number']/len(self.pos_CO2)*100), '  %')
        return 0

    #Auxually
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
        ax_pot.set_title('VDW')
        ax_pot.set_xlabel('distance  A')
        ax_pot.set_ylabel('potential  kJ/mol')
        ax_pot.legend()
        plt.show()
        return 0

    
    # electric potential curve
    def charge_graph(self):
        atom_A = 'Cr1'
        atom_B = 'F1'
        qa = self.params_charge(atom_A)
        qb = self.params_charge(atom_B)

        r_list = np.linspace(2.5, self.cutoff, 100)
        pot_list = [self.charge_function(r**2, qa, qb)*self.k_charge for r in r_list]

        fig_pot = plt.figure()
        ax_pot = fig_pot.add_subplot(111)
        ax_pot.plot(r_list, pot_list, label=atom_A+'-'+atom_B)
        ax_pot.set_title('charge')
        ax_pot.set_xlabel('distance  A')
        ax_pot.set_ylabel('potential  kJ/mol')
        ax_pot.legend()
        plt.show()
        return 0

    
if __name__=='__main__':
    main()