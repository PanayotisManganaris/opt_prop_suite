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
    ecutwfc = 0
    ecutrho = 0
    kpoints = 0
    elements = []
    pps = []
    for line in lines:
        if 'ecutwfc' in line:
            ecutwfc = line.split()[1]
        elif 'ecutrho' in line:
            ecutrho = line.split()[1]
        elif 'kpoints' in line:
            kpoints = line.split()[1]
        elif 'element' in line:
            elements.append(line.split()[2])
        elif '.upf' in line:
            pps.append(line.split()[1])
    inputs.close()

    log = open('run.log','r')
    lines = log.readlines()
    vcr = False
    scf = False
    ph = False
    dm = False
    date = lines[0].split()[1]
    time = lines[0].split()[2]
    for line in lines:
        if 'reached cell relaxation' in line:
            vcr = True
        elif 'reached self consistent field calculation' in line:
            scf = True
        elif 'reached phonon calculation' in line:
            ph = True
        elif 'reached dynamical matrix calculation' in line:
            dm = True
        log.close()

    df.loc[len(df.index)] = [filename, elements, ecutwfc, ecutrho, kpoints, pps, date, time, vcr, scf, ph, dm]
    print(f"added inputs successfully for {filename}")

    return df

outputs = ['filename', 'elements', 'ecutwfc', 'ecutrho', 'kpoints', 'pps', 'date', 'time', 'vcr', 'scf', 'ph', 'dm']

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
