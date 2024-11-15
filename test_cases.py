import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simps
import time
from matplotlib.patches import Patch
import nglview as nv
from twodanalysis import twod_analysis




top = "dopcchol_Charmm.pdb"
traj = "dopcchol_Charmm.pdb"
tpr = "veamos.tpr"


top = "membrane.gro"
traj = "membrane.xtc"

# Creating the class
membrane = twod_analysis(top,
                         traj,
                        tpr=tpr,
                        v_min = -10,
                        v_max = 180,
                        verbose = True,
                        add_radii = True)






lipid_list = list(membrane.lipid_list)
first_lipids = membrane.first_lipids

######### Lipid order 2d code related ########

layers = ["top", "bot", "both"]
lipid_list.remove("CHL1")
nbins = 50
lipids = membrane.chain_info
"""
for layer in layers:
    for key in lipid_list:
        H, edges = membrane.order_histogram(key, layer, nbins, lipids[key])
        print(key, layer, nbins, lipids[key], 0, 180)
        plt.imshow(H,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
        plt.colorbar(cmap = "Spectral")
        plt.savefig(f"{key}_test1_{layer}.png")
        plt.close()
plt.show()
"""
layer = "top"
#mat_top, edges = membrane.all_lip_order("top", nbins,
#                        start = 0,
#                        final = 100,
#                        step = 1)#


#mat_bot, edges = membrane.all_lip_order("bot", nbins,
#                        start = 0,
#                        final = 100,
#                        step = 1)

#mat_both, edges = membrane.all_lip_order("both", nbins,
#                        start = 0,
#                        final = 100,
#                        step = 1)

#plt.imshow(mat_top, cmap = "Spectral")
#plt.show()
#plt.imshow(mat_top, cmap = "Spectral")
#plt.show()




##### Packing deffects related ######


# Adding POPE lipids that are not taken into account
membrane.non_polar_dict["POPE"].append("H101")
membrane.non_polar_dict["POPE"].append("H91")

""" Print selections to test packing deffects with VMD
lipid_polar = {}
for key in first_lipids.keys():
    polar_dict = {"polar":[], "nonpolar": []}
    #print(membrane.non_polar_dict[key])
    atoms = membrane.u.select_atoms(f"resid {first_lipids[key]}")
    names = atoms.names
    for name in names:
        #print(membrane.non_polar_dict[key])
        if name in membrane.non_polar_dict[key]:
            polar_dict["nonpolar"].append(name)
        else:
            polar_dict["polar"].append(name)

    lipid_polar[key] = polar_dict


for key in lipid_polar.keys():
    print(key)
    print(f"###### resname {key} and {membrane.build_name(  lipid_polar[key]['polar'])}    ######\n\n")
    print(f"###### resname {key} and {membrane.build_name(lipid_polar[key]['nonpolar'])}    ######\n\n")

"""


#membrane.visualize_polarity()

#print(membrane.non_polar_dict["POPE"])
#membrane.visualize_polarity()
for ts in membrane.u.trajectory[:20]:
    matrix, matrix_height = membrane.packing_defects(start = 0, final = 10, step =1,nbins = 180, layer = "top", height = True)
    fig,ax = plt.subplots(1,2)
    ax[0].imshow(np.rot90(matrix))
    ax[1].imshow(np.rot90(matrix_height))
    plt.show()




#### Membrane thickness related code ##########

#mat_both, edges = membrane.all_lip_order("both", nbins,
#                        start = 0,
#                        final = 100,
#                        step = 1)

#mat_thi, edges = membrane.thickness(50, start = 0, final = 100, step = 1)
#plt.scatter(mat_thi.flatten(), mat_both.flatten())
#plt.colorbar(cmap = "Spectral")
#plt.show()
#plt.close()
#plt.scatter(mat_both.flatten(), mat_thi.flatten(), alpha = 0.5)
#plt.savefig("corr_thilip.png")
#plt.close()
#print(membrane.lipid_list)
