import MDAnalysis as mda
import sys
sys.path.append("/projects/academic/vmonje/ricardox/github/analysis")
from analysis_class import analysis
import matplotlib.pyplot as plt
import numpy as np

bot = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.gro"
traj = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.xtc"
replica = sys.argv[1]
layer = sys.argv[2]
bot = f"/projects/academic/vmonje/ricardox/200_rna_project/0model/{replica}/centered_prot.gro"
traj = f"/projects/academic/vmonje/ricardox/200_rna_project/0model/{replica}/centered_prot.xtc"
start = 900
nbins = 50
lipid_list = ["DSPC", "POPE", "DODMA", "CHL1"]

#bot = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.gro"
#traj = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.xtc"
#lipid_list = ["DOPC", "DOPE", "POPI1", "POPI2", "CHL1"]
system = analysis(bot, traj, start = start, final = start + 100, lipid_list = lipid_list)

data = system.surface_list(lipids = lipid_list[:-1], layer = layer)

lipids = lipid_list
lipids = {
        "DOPC" : [16,16],
        "DOPE": [16,16],
        "POPI1": [14,16],
        "POPI2": [14,16],
        "DSPC" : [17,17],
        "POPE" : [15,17],
        "DOPS" : [17,17],
        "DODMA" : [17,17],
        "POPS" : [15,17],
        "DSPE" : [17,17],
    }

Hs = []
for key in lipid_list[:-1]:
    print(key)
    H, edges = system.order_histogram(key, layer, nbins, lipids[key], v_min  = 0, v_max = 180)
    #plt.imshow(H,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
    #plt.colorbar(cmap = "Spectral")
    #plt.savefig(f"{key}test{start}.png")
    #plt.close()
    print(H)
    #H[H == np.nan] = 0
    Hs.append(H)

Hs = np.array(Hs)
print(Hs)
print(Hs.shape)
Hs = np.nanmean(Hs, axis = 0)
plt.imshow(Hs,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
plt.colorbar(cmap = "Spectral")
plt.savefig(f"test{start}{layer}{replica}.png")
plt.close()

Hs1 = []
for key in data.keys():
    print(data[key])
    H, X, Y = np.histogram2d(x = data[key]["x"], y = data[key]["y"],bins = nbins, weights = data[key]["z"], range = [[0,180], [0,180]])
    H1, X, Y = np.histogram2d(x = data[key]["x"], y = data[key]["y"],bins = nbins, range = [[0,180], [0,180]])

    H1[H1 == 0] = np.nan
    H_avg = H/H1
    H_avg = np.rot90(H_avg) 
    #plt.imshow(H_avg,cmap = "Spectral", extent = [X[0], X[-1], X[0], X[-1]])
    #plt.colorbar(cmap = "Spectral")
    #plt.savefig(f"{key}height{start}.png")
    #plt.close()
    Hs1.append(H_avg)
    print(H_avg, X, Y)

Hs1 = np.array(Hs1)
print(Hs1)
print(Hs1.shape)
Hs1 = np.nanmean(Hs1, axis = 0)
plt.imshow(Hs1,cmap = "Spectral", extent = [X[0], X[-1], X[0], X[-1]])
plt.colorbar(cmap = "Spectral")
plt.savefig(f"height{start}{replica}{layer}.png")
plt.close()

height = Hs1.flatten()
order = Hs.flatten()

def normal(vector):
    min_v = np.nanmin(vector)
    max_v = np.nanmax(vector)
    normalized = (vector-min_v)/(max_v-min_v)
    return normalized


height = normal(height)
order = normal(order)
print(height)
print(order)

plt.scatter(height, order, alpha =0.2)
plt.xlabel("height")
plt.ylabel("order")
plt.savefig(f"normaltrend{replica}{layer}.png")

#print(data)
