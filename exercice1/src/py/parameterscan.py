import numpy as np
import subprocess
import matplotlib.pyplot as plt

# Parameters
# TODO adapt to what you need (folder path executable input filename)
repertoire = './'  # Path to the compiled code (NB: ./ is not required on Windows)
executable = 'exe'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([400, 800]) # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 600  # TODO: Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps

# Analysis
# TODO insert the values
m = 0  
v = 0
omega = 0  
# add the other variables
# TODO: Insert here the expressions for the exact final solution
x_th  = 0 
y_th  = 0
vx_th = 0
vy_th = 0
"""
... and other parameters
"""

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
vy_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


error = np.zeros(nsimul)

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    xx = data[-1, 1]  # final position, velocity, energy
    yy = data[-1, 2]
    vx = data[-1, 3]
    vy = data[-1, 4]
    En = data[-1, 5]
    vy_list.append(vy)
    # TODO compute the error for each simulation
    error[i] =  100 

lw = 1.5 # line width. TODO: adjust if needed
fs = 16  # font size. TODO: adjust if needed

fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[:, 1], data[:, 2])
ax.set_xlabel('x [m]', fontsize=fs)
ax.set_ylabel('y [m]', fontsize=fs)



# uncomment the following 2 lines if you want debug
#import pdb
#pbd.set_trace()
plt.figure()
plt.loglog(dt, error, 'r+-', linewidth=lw)
plt.xlabel('\Delta t [s]', fontsize=fs)
plt.ylabel('final position error [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, TODO: modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1  # TODO: Modify if needed

plt.figure()
plt.plot(dt**norder, vy_list, 'k+-', linewidth=lw)
plt.xlabel('\Delta t [s]', fontsize=fs)
plt.ylabel('v_y [m/s]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

plt.show()
