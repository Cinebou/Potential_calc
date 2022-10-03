import pandas as pd

# cut difinition
def cut_pdb(x):
    x=x.drop(['A','atom','M','b','c','d'],axis=1)
    return x


# cut the pdb file into csv file to make it readable
def modify_PBDfile(skip_step, row_num, fname = 'input.csv', step = 100):
    f_cols=['A','number','atom','M','x','y','z','b','c','d']
    # re-read the file after the second last endmol
    f_reader = pd.read_csv(fname, delim_whitespace = True, chunksize=100, header=None, names = f_cols, skiprows=skip_step, nrows=row_num)
    
    # cut into pieces to read as csv file
    df =pd.concat((cut_pdb(r) for r in f_reader), ignore_index = True)
    return df

# just read all of them not skip
def read_PDBfile(fname):
    f_cols=['A','number','atom','M','x','y','z','b','c','d']
    # re-read the file after the second last endmol
    f_reader = pd.read_csv(fname, delim_whitespace = True, chunksize=100, header=None, names = f_cols)
    
    # cut into pieces to read as csv file
    df =pd.concat((cut_pdb(r) for r in f_reader), ignore_index = True)
    return df

# read the division line of steps frompdb file
def read_line(fname):
    f_cols=['A','number','atom','M','x','y','z','b','c','d']
    f_reader = pd.read_csv(fname, delim_whitespace = True, chunksize=100, header=None, names = f_cols)
    
    # find 'MODEL as a sign of the beginning of the step
    model_col = []
    for r in f_reader:
         model_col += r.index[r['A']=='MODEL'].tolist()
    return model_col


if __name__ == '__main__':
    start_s = 1
    final_s = 100
    modify_PBDfile(start_s, final_s)