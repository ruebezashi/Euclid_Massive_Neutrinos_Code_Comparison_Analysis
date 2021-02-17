import numpy as np
import MAS_library as MASL
import matplotlib.pyplot as plt
import Pk_library as PKL

#the number of snapshots that are compared
num_snaps = int(input("Please enter the number of snapshots you wish to analyse \n"))

#creates a list containing the paths to the snapshots
snapshots = [None]*(num_snaps)
snapshots_0eV = [None]*num_snaps
mass = input("Enter the neutrino mass in the simulations you wish to analyse (e.g. 0.0eV, 0.15eV...)\n")+"/"
#basic path where all files are placed
path = "/madfs/nucome/"

#the path to the reference run

#ref = input("Please enter the path to the 0.0eV reference run that the snapshots are compared to \n" + path)
ref_mass = path + mass + "Gadget3/output/snapdir_005/snap_005"
ref = "/madfs/nucome/0.0eV/Gadget3/output/snapdir_005/snap_005"
#ref_150 = /madfs/nucome/0.15eV/Gadget3/output/snapdir_005/snap_005

grid     = 256                     #grid size, chose a small one for now because of RAM issues
ptypes   = [1]                     #we investigate the CDM + baryon power spectrum 
MAS      = 'CIC'                   #Cloud-in-Cell
do_RSD   = False                   #dont do redshif-space distortions
axis     = 0                       #axis along which place RSD; not used here
verbose  = True   #whether print information on the progress
BoxSize = 512
threads = 4

#calculate Pk for the 0.0eV reference run
ref_delta = MASL.density_field_gadget(ref, ptypes, grid, MAS, do_RSD, axis, verbose)
ref_delta /= np.mean(ref_delta, dtype = np.float64); ref_delta -= 1.0
ref_Pk = PKL.Pk(ref_delta, BoxSize, axis, MAS, threads, verbose)
ref_k = ref_Pk.k3D
ref_Pk0 = ref_Pk.Pk[:,0]

#calculate Pk for the reference run with given neutrino mass
ref_mass_delta = MASL.density_field_gadget(ref_mass, ptypes, grid, MAS, do_RSD, axis, verbose)
ref_mass_delta /= np.mean(ref_mass_delta, dtype = np.float64); ref_mass_delta -= 1.0
ref_mass_Pk = PKL.Pk(ref_mass_delta, BoxSize, axis, MAS, threads, verbose)
ref_mass_k = ref_mass_Pk.k3D
ref_mass_Pk0 = ref_mass_Pk.Pk[:,0]

#set up the figure
fig, ax = plt.subplots(2,2, figsize = (7,7))
fig.suptitle(mass[:-1])
shaded = np.linspace(0, 100,3)

ax[0,0].loglog(ref_mass_k, ref_mass_Pk0, label = "reference run")
ax[0,1].plot(ref_mass_k, ref_mass_Pk0/ref_Pk0, label = "reference run")

#loop over all snapshots
for i in range(num_snaps):
    
    snapshots[i] = input("Enter the path to snapshot " + str(i+1) + "\n" + path + mass)
    snapshots_0eV[i] = input("Enter the path to the corresponding 0.0eV simulation \n" + path + "0.0eV/")
    name = input("Give a name for the snapshot simulation \n")
    print("Analyzing snapshot " + str(i+1) + "... \n")

    #calculate Pk for 0eV
    z_delta = MASL.density_field_gadget(path+"0.0eV/"+snapshots_0eV[i], ptypes, grid, MAS, do_RSD, axis, verbose)
    z_delta /= np.mean(z_delta, dtype = np.float64); z_delta -= 1.0
    z_Pk = PKL.Pk(z_delta, BoxSize, axis, MAS, threads, verbose)
    z_k = z_Pk.k3D
    z_Pk0 = z_Pk.Pk[:,0]

    #calculate Pk for given neutrino mass
    delta = MASL.density_field_gadget(path + mass + snapshots[i], ptypes, grid, MAS, do_RSD, axis, verbose)
    delta /= np.mean(delta, dtype = np.float64); delta -= 1.0
    Pk = PKL.Pk(delta, BoxSize, axis, MAS, threads, verbose)
    k = Pk.k3D
    Pk0 = Pk.Pk[:,0]



    ax[0,0].loglog(k, Pk0, label = name)
    ax[0,0].set_xlabel(r"$k$ [h/Mpc]")
    ax[0,0].set_ylabel(r"$P(k)$ [Mpc$^3$/h$^3$]")
    ax[0,0].legend()

    ax[0,1].plot(k, Pk0/z_Pk0, label = name)
    ax[0,1].set_xlabel(r"$k$ [h/Mpc]")
    ax[0,1].set_ylabel(r"$P(k)_{0.15eV}/P(k)_{0eV}$")
    ax[0,1].legend()
    ax[0,1].set_xscale("log")
    ax[0,1].set_ylim([0.9,1])
    ax[0,1].set_xlim([0.01,2])

    ax[1,0].plot(k, Pk0/ref_mass_Pk0, label = name)
    ax[1,0].fill_between(shaded, 0.99, 1.01, color = "grey", alpha = 0.3)
    ax[1,0].set_xlabel(r"$k$ [h/Mpc]")
    ax[1,0].set_xscale("log")
    ax[1,0].set_ylabel(r"$P(k)/P(k)_{reference}$")
    ax[1,0].set_ylim([0.9,1.1])
    ax[1,0].set_xlim([0.01,2])
    ax[1,0].legend()

    ax[1,1].plot(k, (Pk0 - ref_mass_Pk0)/ref_mass_Pk0*100, label = name)
    ax[1,1].fill_between(shaded, -1, 1, color = "grey", alpha = 0.3)
    ax[1,1].set_xlabel(r"$k$ [h/Mpc]")
    ax[1,1].set_xscale("log")
    ax[1,1].set_ylabel(r"$(P(k) - P(k)_{reference}) / P(k)_{reference}$ [%]")
    ax[1,1].set_ylim([-3,3])
    ax[1,1].set_xlim([0.01,2])
    ax[1,1].legend()

plt.subplots_adjust(wspace = 0.4, hspace = 0.2)
fig.savefig("analysis.pdf")
# Compute the effective number of particles/mass in each voxel
#delta1 = MASL.density_field_gadget(snapshot1, ptypes, grid, MAS, do_RSD, axis, verbose)
#delta2 = MASL.density_field_gadget(snapshot2, ptypes, grid, MAS, do_RSD, axis, verbose)

# compute density contrast: delta = rho/<rho> - 1
#delta1 /= np.mean(delta1, dtype=np.float64);  delta1 -= 1.0
#delta2 /= np.mean(delta2, dtype=np.float64);  delta2 -= 1.0


# compute power spectrum
#Pk1 = PKL.Pk(delta1, BoxSize, axis, MAS, threads, verbose)
#Pk2 = PKL.Pk(delta2, BoxSize, axis, MAS, threads, verbose)
# Pk is a python class containing the 1D, 2D and 3D power spectra, that can be retrieved as

# 1D P(k)
#k1D1      = Pk1.k1D
#Pk1D1     = Pk1.Pk1D
#Nmodes1D1 = Pk1.Nmodes1D

#k1D2 = Pk2.k1D
#Pk1D2 = Pk2.Pk1D
#Nmodes1D2 = Pk2.Nmodes1D

# 2D P(k)
#kpar1     = Pk1.kpar
#kper1     = Pk1.kper
#Pk2D1     = Pk1.Pk2D
#Nmodes2D1 = Pk1.Nmodes2D

#kpar2 = Pk2.kpar
#kper2 = Pk2.kper
#Pk2D2 = Pk2.Pk2D
#Nmodes2D2 = Pk2.Nmodes2D

# 3D P(k)
#k1       = Pk1.k3D
#Pk01     = Pk1.Pk[:,0] #monopole
#Pk21     = Pk1.Pk[:,1] #quadrupole
#Pk41     = Pk1.Pk[:,2] #hexadecapole
#Pkphase1 = Pk1.Pkphase #power spectrum of the phases
#Nmodes1  = Pk1.Nmodes3D

#k2 = Pk2.k3D
#Pk02 = Pk2.Pk[:,0]
#Pk22 = Pk2.Pk[:,1]
#Pk42 = Pk2.Pk[:,2]
#Pkphase2 = Pk2.Pkphase
#Nmodes2 = Pk2.Nmodes3D


#file = np.loadtxt("/data/jadame/nucome/0.0eV/output/test1024hybrid_pk007_delta.dat")

#fig = plt.figure()
#plt.plot(k1, Pk01/Pk02)
#plt.xscale("log")
#plt.xlabel(r"k [$h/Mpc$]")
#plt.ylabel("P(k) gevolution / P(k) reference run") 
#plt.ylim([0.96, 1.01])
#fig.savefig("Pk0_"+snapshot_name+".pdf")
#plt.close(fig)

#fig = plt.figure()
#plt.loglog(k1, Pk01)
#plt.loglog(k2, Pk02)
#plt.loglog(file[:,0], file[:,1]/file[:,0]**3*2*np.pi**2)
#plt.xlabel(r"k [$h/Mpc$]")
#plt.ylabel(r"P(k) [$Mpc^3/h^3$]")
#fig.savefig("Pk0_"+snapshot_name+"_comparison.pdf")
