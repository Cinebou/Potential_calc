import pandas as pd
import csv
from split import modify_PBDfile, read_line, read_PDBfile
from pprint import pprint

# judge whether it is in the sphere of the tetrahidra structure
def sphere_judge(r, xC,yC,zC, x_ST, y_ST, z_ST):
    dx = abs(xC - x_ST)
    dy = abs(yC - y_ST)
    dz = abs(zC - z_ST)

    # make the calcualtion Faster
    if (dx > r) or (dy > r) or (dz > r):
        return False

    dr = (dx**2 + dy**2 + dz**2)**(1/2)
    if dr <= r:
        return True
    else:
        return False


# count the histogram of the number of CO2 in the cell
def make_hist(counter):
    num_ST = len(counter)
    max_num_in_ST = 7
    hist = [0 for r in range(max_num_in_ST)]
    for c in counter:
        hist[c] += 1/num_ST*100 # normalized value in % representation
    return hist

    

def countST():
    # reading the position of the center of tetrahidra structure
    cavity_file = './Configure/MIL101/tetra_center_PBC.csv'
    CO2_file = './force_field/VMD/298.15_100000.pdb'
    tetra = pd.DataFrame(pd.read_csv(cavity_file))

    # Radius of ST
    R = 4
    DataCo2 = read_PDBfile(CO2_file) 

    InOutSmall = [0 for i in range(int(len(DataCo2)/3))]
    # calcuate for each CO2 molecluaes in the cell 
    for index, posCO2 in DataCo2.iterrows():
        if (index-1)%3 == 0: # calculate only C atom in the trajectory
            x = float(posCO2[1])
            y = float(posCO2[2])
            z = float(posCO2[3])
                
            # iterate for each small cavity in the cell
            for k, posST in tetra.iterrows():
                xST = float(posST['center x'])
                yST = float(posST['center y'])
                zST = float(posST['center z'])

                # judge whether it is in the ST or not
                if sphere_judge(R, x,y,z, xST, yST, zST):
                    InOutSmall[int(index/3)] += 1
    print('number of molecules in the cell   : ', len(InOutSmall))
    with open('./Results/MIL101/small_cav.csv', 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(InOutSmall)

    

if __name__ == '__main__':
    countST()