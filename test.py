import MDAnalysis as mda
import sys
sys.path.append("/projects/academic/vmonje/ricardox/github/analysis")
from analysis_class import analysis
import matplolib.pyplot as plt


top = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.gro"
traj = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.xtc"

#top = "/projects/academic/vmonje/ricardox/200_rna_project/0model/rep2/centered_prot.gro"
#traj = "/projects/academic/vmonje/ricardox/200_rna_project/0model/rep2/centered_prot.xtc"

system = analysis(top, traj, start = 500, final = 600, lipid_list = ["DOPC","DOPE", "POPI1", "POPI2", "CHL1"])

#data = system.surface_list(lipids = ["DOPC", "DOPE", "POPI1", "POPI2"])

H, edges = system.order_histogram("DOPC", "top", 50, [16,16], v_min  = 0, v_max = 180)
print(H, edges)

plt.imshow(H, extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
plt.savefig("test.png")
#print(data)
