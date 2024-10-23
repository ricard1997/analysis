import MDAnalysis as mda
import sys
sys.path.append("/projects/academic/vmonje/ricardox/github/analysis")
from analysis_class import analysis



top = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.gro"
traj = "/projects/academic/vmonje/ricardox/trajetories_mlkl/rep3/waterno/centered_wowatertrep1.xtc"


#system = analysis(top, traj, start = 500, final = 600, lipid_list = ["DOPC", "DOPE", "POIP1", "POIP2", "CHL1"])

data = system.surface_list(lipids = ["DOPC", "DOPE", "POPI1", "POPI2"])
print(data)
