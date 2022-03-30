import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd

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
                E.append(data[4])

        output.close()

    s = np.linspace(1, len(E), len(E))

    plt.figure()
    plt.plot(s,E)
    plt.xlabel("ionic step")
    plt.ylabel("energy (eV)")
    plt.title("Energy Convergence")
    plt.savefig("E.png")

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
#
# data = subprocess.check_output("grad2 OUTCAR",shell=True)
# shift = data.split().index(b'1')
#
# E = np.array(data.split()[2+shift:-1:18]).astype(float);
# avgF = np.array(data.split()[8+shift:-1:18]).astype(float);
# maxF = np.array(data.split()[10+shift:-1:18]).astype(float);
# s = np.linspace(1,len(E),len(E))
#
# plt.figure()
# plt.plot(s,maxF,label="Max F")
# plt.plot(s,avgF,label="Avg F")
# plt.xlabel("ionic step")
# plt.ylabel("force (eV/Angstrom)")
# plt.title("Force Convergence")
# plt.legend()
# plt.savefig("F.png")
#
# plt.figure()
# plt.plot(s,E)
# plt.xlabel("ionic step")
# plt.ylabel("energy (eV)")
# plt.title("Energy Convergence")
# plt.savefig("E.png")
#
# plt.figure()
# plt.plot(sP,P)
# plt.xlabel("ionic step")
# plt.ylabel("external pressure (kB)")
# plt.ylim([0,500])
# plt.title("Pressure Convergence")
# plt.savefig("convergence_P.png")
#
# plt.close('all')
