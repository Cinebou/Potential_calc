import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main():
    file_name = './Results/298.15_100000.csv'
    num_atom_in_molecule = 3
    hist_digit = 40
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

    plt.hist(tot_heat_list)
    plt.show()

if __name__ == '__main__':
    main()