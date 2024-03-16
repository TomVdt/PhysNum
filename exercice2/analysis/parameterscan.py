import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb


# TODO adapt to what you need (folder path executable input filename)

repertoire = './'
repertoire = "C:\\Users\\Administrator\\physnum\\2024\\EX2\\STUDENT\\" # change as you need

 # Path to the compiled code (NB: ./ is not required on Windows)
executable = "l.exe"  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([300, 1000, 5000, 10000]) # TODO change as you need
nsimul = len(nsteps)  # Number of simulations to perform



dt = 1/nsteps # define it as you need

# Analysis
# TODO: Insert here the expressions for the exact final solution


paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  
vy_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)
plt.figure(figsize=(10, 6))
xx_list = []
yy_list = []
for i in range(nsimul):  # Iterate through the results of all simulations
    
    data = np.loadtxt(outputs[i])  
    t  = data[:, 0]
    xx = data[-1, 1] 
    yy = data[-1, 2]  
    xx_list.append(xx)
    yy_list.append(yy)
    # compute the error if you can ..
    error = 0
    

plt.plot(dt, xx_list) # this is an example, change it if you need
plt.xlabel("$dt$", fontsize=20)  # LaTeX for theta
plt.ylabel(r'${\theta}$', fontsize=20)  # Lating the y-axis label with fontsize 20
plt.xticks(fontsize=15)  
plt.yticks(fontsize=15)  
plt.grid(True)

# TODO add also the other plots that you need


  
plt.show()  
    

