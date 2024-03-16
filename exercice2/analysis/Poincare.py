# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:22:14 2024

@author: Administrator
"""

import numpy as np
import subprocess
import matplotlib.pyplot as plt

# Parameters
# TODO adapt to what you need (folder path executable input filename)

repertoire = './'
repertoire = "C:\\Users\\Administrator\\physnum\\2024\\EX2\\STUDENT\\"

 # Path to the compiled code (NB: ./ is not required on Windows)
executable = "lolo.exe"  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


thetadot = np.array([16.0, 10]) # TODO change
nsimul   = len(thetadot)        # Number of simulations to perform

nsteps   = 100 #change according what you need


paramstr  = 'thetadot0'  # Parameter name to scan
paramstr2 = 'nsteps'  # Parameter name to scan
paramstr3 = 'sampling'  # Parameter name to scan
param     = thetadot  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
vy_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    
    cmd = f"{repertoire}{executable} {input_filename} {paramstr3}={nsteps:.15g} {paramstr2}={nsteps:.15g} {paramstr}={param[i]:.15g} output={output_file}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)
plt.figure(figsize=(10, 6))

for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t = data[:, 0]

    xx = data[:, 1] 
    yy = data[:, 2]   
    plt.plot(np.mod(xx+np.pi,2*np.pi)-np.pi, yy, '.')
    
plt.xlabel(r'$\theta$', fontsize=20)  # LaTeX for theta
plt.ylabel(r'$\dot{\theta}$', fontsize=20)  # Lating the y-axis label with fontsize 20

plt.xticks(fontsize=15)  
plt.yticks(fontsize=15)  

plt.grid(True)  
plt.show()  
    

