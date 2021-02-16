import numpy as np
import MAS_library as MASL
import matplotlib.pyplot as plt
import Pk_library as PKL

#the number of snapshots that are compared
num_snaps = int(input("Please enter the number of snapshots you wish to analyse \n"))

#creates a list containing the paths to the snapshots
snapshots = [None]*(num_snaps)

#basic path where all files are placed
path = "/madfs/nucome/"

#the path to the reference run
ref = input("Please enter the path to the reference run that the snapshots are compared to \n" + path)


grid     = 256                     #grid size
ptypes   = [1]                     #we only are interested in the CDM + baryon spectrum
MAS      = 'CIC'                   #Cloud-in-Cell
do_RSD   = False                   #dont do redshif-space distortions
axis     = 0                       #axis along which place RSD; not used here
verbose  = True   #whether print information on the progress
BoxSize = 512
threads = 4

#calculate Pk for the reference run
ref_delta = MASL.density_field_gadget(path+ref, ptypes, grid, MAS, do_RSD, axis, verbose)
ref_delta /= np.mean(ref_delta, dtype = np.float64); ref_delta -= 1.0
ref_Pk = PKL.Pk(ref_delta, BoxSize, axis, MAS, threads, verbose)
ref_k = ref_Pk.k3D
ref_Pk0 = ref_Pk.Pk[:,0]

#set up the figure
fig, ax = plt.subplots(2,2, figsize = (7,7))
shaded = np.linspace(0, 100,3)
ax[0,0].loglog(ref_k, ref_Pk0, label = "reference run")

#loop over all snapshots
for i in range(num_snaps):
    
    snapshots[i] = input("Enter the path to snapshot " + str(i+1) + "\n" + path)
    name = input("Give a name for the snapshot simulation \n")
    print("Analyzing snapshot " + str(i+1) + "... \n")
    
    delta = MASL.density_field_gadget(path + snapshots[i], ptypes, grid, MAS, do_RSD, axis, verbose)
    delta /= np.mean(delta, dtype = np.float64); delta -= 1.0
    Pk = PKL.Pk(delta, BoxSize, axis, MAS, threads, verbose)
    k = Pk.k3D
    Pk0 = Pk.Pk[:,0]

    ax[0,0].loglog(k, Pk0, label = name)
    ax[0,0].set_xlabel(r"$k$ [h/Mpc]")
    ax[0,0].set_ylabel(r"$P(k)$ [Mpc$^3$/h$^3$]")
    ax[0,0].legend()

    ax[1,0].plot(k, Pk0/ref_Pk0, label = name)
    ax[1,0].fill_between(shaded, 0.99, 1.01, color = "grey", alpha = 0.3)
    ax[1,0].set_xlabel(r"$k$ [h/Mpc]")
    ax[1,0].set_xscale("log")
    ax[1,0].set_ylabel(r"$P(k)/P(k)_{reference}$")
    ax[1,0].set_ylim([0.9,1.1])
    ax[1,0].set_xlim([0.01,2])
    ax[1,0].legend()

    ax[1,1].plot(k, (Pk0 - ref_Pk0)/ref_Pk0, label = name)
    ax[1,1].fill_between(shaded, -0.01, 0.01, color = "grey", alpha = 0.3)
    ax[1,1].set_xlabel(r"$k$ [h/Mpc]")
    ax[1,1].set_xscale("log")
    ax[1,1].set_ylabel(r"$(P(k) - P(k)_{reference}) / P(k)_{reference}$")
    ax[1,1].set_ylim([-0.1,0.1])
    ax[1,1].set_xlim([0.01,2])
    ax[1,1].legend()

plt.subplots_adjust(wspace = 0.4, hspace = 0.2)
fig.savefig("analysis.pdf")
