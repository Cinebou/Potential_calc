import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from calc_pot_total import Potential
import pandas as pd

def main():
    CO2_file = './force_field/VMD/298.15_100000.pdb'
    frame_file = './force_field/MIL-101_all.pdb'
    mil101 = Potential(CO2_file, frame_file)
    plot_3D_frame(mil101)
    #write_PDB(mil101)


def write_PDB(pot_mil101):
    pot_mil101.pos_frame.to_csv('./force_field/VMD/PBC_framework.pdb',sep='\t',index=False)


def plot_3D_frame(pot_mil101):
    x = pot_mil101.pos_frame['x']
    y = pot_mil101.pos_frame['y']
    z = pot_mil101.pos_frame['z']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(x, y, z)
    ax.set_title("MIL-101 structure")
    plt.show()



if __name__ == '__main__':
    main()