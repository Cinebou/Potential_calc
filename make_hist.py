import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
import csv

def main(file_name, legend_name,ax):
    num_atom_in_molecule = 3
    hist_bin = 45
    max_heat = 45

    res_pot = pd.read_csv(file_name)
    res_pot = res_pot[0:(len(res_pot)//3)*3]
    atom_num = len(res_pot)
    mol_num = int(atom_num / num_atom_in_molecule)
    print('number of atoms in the file is  : ',atom_num)
    print('all heat of adsorption  : ',sum(res_pot['Total_pot'])/mol_num)

    # CO2 consists of 3 atoms, you need to take the sum of the three atoms for each molecules
    tot_heat_list = [0 for m in range(mol_num)]
    for res_index, res_row in res_pot.iterrows():
        tot_heat_list[int((res_row['atom_num']-1)//3)] += res_row['Total_pot']

    pot_hist, pot_bins = np.histogram(tot_heat_list, hist_bin,range=(0, max_heat),density = True)
    pot_bins_correct = [(pot_bins[i]+pot_bins[i-1])/2.0 for i in range(1,len(pot_bins))]

    ax.scatter(pot_bins_correct,pot_hist)
    ax.plot(pot_bins_correct,pot_hist, linestyle = "--",label=legend_name)
    ax.legend()
    ax.set_xlabel('Interaction energy: kJ/mol', fontsize = 14)
    ax.set_ylabel('fraction of molecule: a.u.', fontsize = 14)
    return 0

def small_cav_dist(file_name):
    num_atom_in_molecule = 3
    hist_bin = 45
    max_heat = 45

    res_pot = pd.read_csv(file_name)
    res_pot = res_pot[0:(len(res_pot)//3)*3]
    atom_num = len(res_pot)
    mol_num = int(atom_num / num_atom_in_molecule)
    print('number of atoms in the file is  : ',atom_num)
    print('all heat of adsorption  : ',sum(res_pot['Total_pot'])/mol_num)

    # CO2 consists of 3 atoms, you need to take the sum of the three atoms for each molecules
    tot_heat_list = [0 for m in range(mol_num)]
    for res_index, res_row in res_pot.iterrows():
        tot_heat_list[int((res_row['atom_num']-1)//3)] += res_row['Total_pot']
    # small cavity analysis
    Scav_heat_list_all = [0 for m in range(mol_num)]
    Scav = [0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]
    index=0
    for cav in Scav:
        if  cav==1:
            Scav_heat_list_all[index] = tot_heat_list[index]
        index+=1
    Scav_heat_list = [s for s in Scav_heat_list_all if s!=0]
    S_pot_hist, S_pot_bins = np.histogram(Scav_heat_list, hist_bin,range=(0, max_heat),density = True)
    S_pot_bins_correct = [(S_pot_bins[i]+S_pot_bins[i-1])/2.0 for i in range(1,len(S_pot_bins))]

    fig2 = plt.figure(figsize=(12, 8))
    ax2 = fig2.add_subplot(111)
    ax2.scatter(S_pot_bins_correct,S_pot_hist)
    ax2.plot(S_pot_bins_correct,S_pot_hist, linestyle = "--")
    ax2.set_xlabel('Interaction energy: kJ/mol', fontsize=18)
    ax2.set_ylabel('fraction of molecule: a.u.', fontsize=18)
    

def show_results(file_name):
    num_atom_in_molecule = 3
    res_pot = pd.read_csv(file_name)
    res_pot = res_pot[0:(len(res_pot)//3)*3]
    atom_num = len(res_pot)
    mol_num = int(atom_num / num_atom_in_molecule)
    print("overall heat of adsorption    : ",sum(res_pot['Total_pot']/mol_num))
    print("VDW potential of CO2-CO2      : ",sum(res_pot['LJ_CO2']/mol_num))
    print("VDW potential of frame-CO2    : ",sum(res_pot['LJ_frame']/mol_num))
    print("Electric potential CO2-CO2    : ",sum(res_pot['Elec_CO2']/mol_num))
    print("Electric potential CO2-frame  : ",sum(res_pot['Elec_frame']/mol_num))


if __name__ == '__main__':
    file_name_1 = './Results/MIL101/MIL101_298_1bar.csv'
    file_name_2 = './Results/NH2-MIL101_298/MIL-101-NH2_1bar_298.csv'
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    main(file_name_1,'MIL101',ax)
    main(file_name_2,'MIL101-NH2',ax)
    plt.savefig('./Fig/1bar_mil101_nh2mil101.tiff')

    print(file_name_1)
    show_results(file_name_1)
    print(file_name_2)
    show_results(file_name_2)
    plt.show()