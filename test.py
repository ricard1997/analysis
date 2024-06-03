import MDAnalysis as mda
import sys
sys.path.append("/projects/academic/vmonje/ricardox/github/analysis")
from analysis_class import analysis
import matplotlib.pyplot as plt
import numpy as np

top = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.gro"
traj = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.xtc"

#top = "/projects/academic/vmonje/ricardox/200_rna_project/0model/rep2/centered_prot.gro"
#traj = "/projects/academic/vmonje/ricardox/200_rna_project/0model/rep2/centered_prot.xtc"
start = 1900
system = analysis(top, traj, start = start, final = start + 100, lipid_list = ["DOPC","DOPE", "POPI1", "POPI2", "CHL1"])

#data = system.surface_list(lipids = ["DOPC", "DOPE", "POPI1", "POPI2"])

lipids = ["DOPC", "DOPE", "POPI1", "POPI2"]
lipids = {
        "DOPC" : [16,16],
        "DOPE": [16,16],
        "POPI1": [14,16],
        "POPI2": [14,16],
    }
Hs = []
for key in lipids.keys():
    print(key)
    H, edges = system.order_histogram(key, "bot", 50, lipids[key], v_min  = 0, v_max = 180)
    plt.imshow(H,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
    plt.colorbar(cmap = "Spectral")
    plt.savefig(f"{key}test{start}.png")
    plt.close()
    print(H)
    #H[H == np.nan] = 0
    Hs.append(H)

Hs = np.array(Hs)
print(Hs)
print(Hs.shape)
Hs = np.nanmean(Hs, axis = 0)
plt.imshow(Hs,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
plt.colorbar(cmap = "Spectral")
plt.savefig(f"test{start}.png")
plt.close()



#print(data)
