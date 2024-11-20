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
import imageio
import os



top = "dopcchol_Charmm.pdb"
traj = "dopcchol_Charmm.pdb"
tpr = "veamos.tpr"

#top = "../../centered_prot.gro"
#traj = "../../centered_prot.xtc"
#tpr = "../../veamos.tpr"

top = "membrane.gro"
traj = "membrane.xtc"

# Creating the class

membrane = twod_analysis(top,
                         traj,
                        tpr=tpr,
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
#membrane.non_polar_dict["POPE"].append("H101")
#membrane.non_polar_dict["POPE"].append("H91")

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
"""
membrane.visualize_polarity()
membrane.non_polar_dict["POPI24"].append("H91")
membrane.non_polar_dict["POPI24"].append("H101")
membrane.non_polar_dict["POPI15"].append("H91")
membrane.non_polar_dict["POPI15"].append("H101")
membrane.visualize_polarity()
count = 0
filenames = []
for ts in membrane.u.trajectory[51::3]:
    matrix, matrix_height = membrane.packing_defects(
                                    nbins = 200,
                                    layer = "bot",
                                    height = True,

                                    )
    fig,ax = plt.subplots(1,2)
    ax[0].imshow(np.rot90(matrix))
    ax[1].imshow(np.rot90(matrix_height))
    ax[0].set_title(f"Frame {count}")
    plt.savefig(f"frame_{count}.png")

    filenames.append(f"frame_{count}.png")
    count += 1

with imageio.get_writer("gif_packing.gif", mode = "I", duration = 0.3) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

for filename in filenames:
    os.remove(filename)

"""




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




###### Code to test apl ######


voronoi_dict = membrane.voronoi_apl(layer = "top")


colors = {"DODMA": "blue",
        "POPE": "red",
        "DSPC": "green",
        "CHL1": "yellow"
}
"""
plt.scatter(*voronoi_dict["points"].T, c =[colors[lipid] for lipid in voronoi_dict["resnames"]])
count = 0
for vertices in voronoi_dict["vertices"]:

    temp = vertices
    temp = np.concatenate([temp, [temp[0]]])
    plt.plot(*temp.T, c = colors[voronoi_dict["resnames"][count]])
    count+=1
print(membrane.print_dict(voronoi_dict["apl"]))
plt.show()
"""

#membrane.map_voronoitest(voronoi_dict["vertices"], voronoi_dict["areas"], 100, [membrane.v_min, membrane.v_max, membrane.v_min, membrane.v_max])
#membrane.map_voronoi(voronoi_dict["points"], voronoi_dict["areas"], 300, [membrane.v_min, membrane.v_max, membrane.v_min, membrane.v_max])

membrane.grid_apl(layer = "top", start = 0, final = 100, step = 1, lipid_list = None)
