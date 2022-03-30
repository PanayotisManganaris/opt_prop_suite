import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
import re

def isadigit(string):
    return bool(re.match(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)', string))

def makePlot():
    filepath = "*.vc-relax.out"
    txt = glob.glob(filepath)
    for file in txt:
        output = open(file, "r")
        lines = output.readlines()
        E = []
        for line in lines:
            if '!' in line:
                data = line.split()
                # if isadigit(data[3]):
                E.append(data[4])

        output.close()

    # print(E)
    E = [float(e) for e in E]
    # pd.series(E)
    # print(E)
    # E = [e for e in E if e.isnumeric()]
    s = np.linspace(1, len(E), len(E))

    plt.figure()
    plt.plot(s,E)
    plt.xlabel("ionic step")
    plt.ylabel("energy (Ry)")
    plt.title("Energy Convergence")
    plt.savefig("E.png")
    plt.close()

def addInputs(filename, df):

    inputs = open("inputs.yaml","r")
    lines = inputs.readlines()
    for i, line in enumerate(lines):
        if 'ecutwfc' in line:
            ecutwfc = line.split()[1]
        elif 'kpoints' in line:
            kpoints = line.split()[1]
        elif 'pp' in line:
            pp1 = lines[i+1].split()[1]
            pp2 = lines[i+2].split()[1]
    print(filename,ecutwfc,kpoints,pp1,pp2)
    df.loc[len(df.index)] = [filename,ecutwfc,kpoints,pp1,pp2]
    print("added inputs")
    inputs.close()

    return df

outputs = ['dir','ecutwfc','kpoints','pp1','pp2']

df = pd.DataFrame(columns=outputs)

for filename in os.listdir('RUNS'):

    f = os.path.join('RUNS',filename)
    if os.path.isdir(f):
        os.chdir(f)
        try:
            makePlot()
            df = addInputs(filename, df)
        except:
            print("failed for "+filename)
        os.chdir('../..')
df.to_csv("inputs.csv")
