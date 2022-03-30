import plotly.graph_objects as go
import numpy as np
import os
import glob
import pandas as pd
import re

def isadigit(string):
    return bool(re.match(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)', string))

def makePlot():
    filepath = os.path.expanduser("/home/panos/src/opt_prop_suite/*.vc-relax.out")
    txt = glob.glob(filepath)
    print(txt)
    for file_name in txt:
        output = open(file_name, "r")
        lines = output.readlines()
        E = []
        for line in lines:
            if '!' in line:
                data = line.split()
                # if isadigit(data[3]):
                E.append(data[4])

        output.close()

    print(E)
    E = [float(e) for e in E]
    # pd.series(E)
    # print(E)
    # E = [e for e in E if e.isnumeric()]
    s = np.linspace(1, len(E), len(E))

    #pick plot that works... f.write_image doesn't work anyway....
#    plt.figure()
#    plt.plot(s,E)
#    plt.xlabel("ionic step")
#    plt.ylabel("energy (Ry)")
#    plt.title("Energy Convergence")
#    plt.savefig("E.png")
#    plt.close()

    f = go.Figure(go.Scatter(x=s,y=E))
    f.update_layout(
        title = "Energy Convergence",
        xaxis_title = "ionic steps",
        yaxis_title = "energy (Ry)"
    )
    f.write_image("e.png")
    
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

# for filename in os.listdir('RUNS'):
# 
#     f = os.path.join('RUNS',filename)
#     if os.path.isdir(f):
#         os.chdir(f)
#         makePlot()
#         df = addInputs(filename, df)
#         #except:
#         #    print("failed for "+filename)
#         os.chdir('../..')
df.to_csv("inputs.csv")
