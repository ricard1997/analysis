import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simps
import time
from matplotlib.patches import Patch





class twod_analysis:
    def __init__(
                self,
                top,
                traj,
                lipid_list = None,
                tpr = None,
                info = False,
                guess_chain_l = True,
                chain_info = None,
            ):



        if tpr:
            self.u = mda.Universe(tpr, traj)
        else:
            self.u = mda.Universe(top, traj)
        #print(self.u.resids)

        if not lipid_list: # Select only elements of the membrane
            self.memb = self.u.select_atoms("all and not protein and not (resname URA or resname GUA or resname ADE or resname CYT)")
            self.lipid_list = set(self.memb.residues.resnames)

        self.working_lip = {
                                "CHL1" : {"head" :"O3", "charge" : 0},
                                "DODMA" : {"head" :"N1", "charge" : -0.21},
                                "DSPC" : {"head" :"P", "charge" : 1.1},
                                "POPE" : {"head" :"P", "charge" : 1.1},
                                "DOPS" : {"head" :"P", "charge" : 0.1},
                                "POPS" : {"head" :"P", "charge" : 0.1},
                                "DSPE" : {"head" :"P", "charge" : 1.3},
                                "DOPC" : {"head" :"P", "charge" : 1.3},
                                "DOPE" : {"head" :"P", "charge" : 1.3},
                                "POPI1" : {"head" :"P", "charge" : 1.3},
                                "POPI2" : {"head" :"P", "charge" : 1.3},
                            } #List of known lipids and lipids head people usually use to work

        self.chain_info = chain_info
        if guess_chain_l: # Guess the chain lenght of lipids. Chain sn2 start with C2 and chain sn1 start with C3
            self.chain_info = {}
            for lipid in self.lipid_list:
                first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
                actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
                actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
                actual_sn1 = actual_sn1.names
                actual_sn2 = actual_sn2.names
                self.chain_info[lipid] = [len(actual_sn1) - 2, len(actual_sn2) - 2]

        self.all_head = self.u.select_atoms(self.build_resname(self.lipid_list) + " and name P")


    @staticmethod
    def get_individual(lista
                    ):
        r"""This function gets a list with a specific carbon (e.g. C34 or C22)
        and its respective hidrogens (e.g. H4X, H4Y). It computes the vectors
        that connect the carbons and the hydrogens and computes the :math:`cos(\theta)^2`, where :math:`\theta` is
        the angle between each vector and the z-axis. Finally, this function returns a vector with the individual (per lipid)
        :math:`\braket{cos(\theta)^2}`, where the mean is computed over the hydrogens of each carbon.

        Parameters
        ----------
        lista : list
            Vector of the shape :math:`[C*i, HiX,HiY, HiZ]`, the minimun len is 2 (when the carbon
            only have one hydrogen) and the maximun is 4 (when there is three hydrogens)
            Note: If there is N lipids, there will be N carbons :math:`C*i`, and the i represents
            the position of the carbon in the lipid tail.

        Returns
        -------
        order : array(n_lipids)
            Float with the mean of :math:`\braket{cos(\theta)^2}`

        Notes
        -----
        The average of the angle of the i-th carbon for all the lipids in the selection is computed
        as follows:

        .. math:: \braket{cos(\theta_i)^2}

        where :math:`\theta_i` is the angle between the z- axis and the vector that connects the i-th carbon and the hydrogen.



        """
        angles = [] # Store the angles for the working carbon
        for i in (range(len(lista)-1)): # Accounts for variable number of list (Change if the carbon has or not double bonds)
            vectores = lista[i+1].positions - lista[0].positions # Hidrogen - Carbons; output of shape (n_lipids, 3)
            costheta = vectores[:,2]**2/np.linalg.norm(vectores, axis = 1)**2 # Compute the costheta^2
            angles.append(costheta) # dim (n_lipids,)
        angles = np.array(angles) # dim ((1, 2 or 3),n_lipids)
        #print("angles", angles.shape)
        angles = np.mean(angles, axis = 0) # output is dim n_lipids, it means the cos^2(theta) or the Carbon passed for each lipid
        return angles


    # Get the cos^2(theta) for each carbon in the selection, for sn1
    def individual_order_sn1(self, sel, lipid, n_chain):
        r"""

        Code to loop over the number of carbons in the lipid tail and get a list with the carbon and its
        hydrogens for each carbon in the lipid tail: :math:`[C3i, HiX, HiY, ...]`. This list is passed to get_vectors which
        return the averages of each i-th carbon. This code returns an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

        Parameters
        ----------
        lipid : str
            Name of the lipid to compute the order parameters
        n_chain : int
            Number of carbons in the lipid tail sn1

        Returns
        -------
        chains : ndarray
            Vector dim n_chains with the mean :math:`\braket{cos(\theta_i)^2}`

        Notes
        -----
        The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

        .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}]

        The index starts at 2 because that is the carbon the lipid tail starts with.

        """

        # Define list to store the chain cos^2(theta)
        chains = []

        # Loop over carbons
        for i in range(n_chain):
            # Define selections for H and C in the chain
            print(f"Value of the chain {i} sn1")
            selections = [
                            f"name C3{i+2}",
                            f"name H{i+2}X and not name HX",
                            f"name H{i+2}Y and not name HY",
                            f"name H{i+2}Z and not name HZ"
                        ]
            #print(selections)


            # Define a list to store atoms
            lista = []

            for selection in selections:
                atoms = sel.select_atoms(selection)


                if atoms.n_atoms != 0:
                    lista.append(atoms)
            # Call get_individual that computes the cos^2(theta) for each carbon.
            chains.append(self.get_individual(lista))
            #print(i, self.get_individual(lista).shape, self.get_individual(lista))
        chains = np.array(chains) # Expect array of dim (n_chain, n_lipids)
        return chains


    # Get the cos^2(theta) for each carbon in the selection, for sn2
    def individual_order_sn2(self, sel, lipid, n_chain):

        r"""

        Code to loop over the number of carbons in the lipid tail and get a list with the carbon and its
        hydrogens for each carbon in the lipid tail: :math:`[C2i, HiX, HiY, ...]`. This list is passed to get_vectors which
        return the averages of each i-th carbon. This code returns an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

        Parameters
        ----------
        lipid : str
            Name of the lipid to compute the order parameters
        n_chain : int
            Number of carbons in the lipid tail sn1

        Returns
        -------
        chains : ndarray
            Vector dim n_chains with the mean :math:`\braket{cos(\theta_i)^2}`

        Notes
        -----
        The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

        .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}]

        The index starts at 2 because that is the carbon the lipid tail starts with.

        """
        # Define list to store the chain cos^2(theta)
        chains = []
        # Loop over carbons
        max_v = 0
        for i in range(n_chain):
            # Define selections for H and C in the chain
            #print(f"Value of the chain {i} sn2")
            selections = [
                            f"name C2{i+2}",
                            f"name H{i+2}R and not name HR",
                            f"name H{i+2}S and not name HS",
                            f"name H{i+2}T and not name HT"
                        ]
            if lipid == "POPE" or lipid == "POPS" or lipid == "POPI1" or lipid == "POPI2":
                if selections[0] == "name C29":
                    selections[1] = "name H91"
                if selections[0] == "name C210":
                    selections[1] = "name H101"
            # Define a list to store atoms
            lista = []

            for selection in selections:
                atoms = sel.select_atoms(selection)


                if atoms.n_atoms != 0:
                    lista.append(atoms)
            # Call get_individual that computes the cos^2(theta) for each carbon.
            angles = self.get_individual(lista)
            if len(angles) > max_v:
                max_v = len(angles)
            chains.append(angles)


        chains = np.array(chains) # Expect array of dim (n_chain, n_lipids)
        return chains


    @staticmethod
    def count_order(data, min_lenght, n_chain):
        columns = ["index"]
        carbons_sn2 = False
        try:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain[0]) ]
            carbons_sn2 = [f"sn2-{i+2}" for i in range(n_chain[1])]
            columns = columns + carbons_sn1 + carbons_sn2
        except:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain) ]
            columns = columns + carbons_sn1


        df = pd.DataFrame(data, columns = columns)
        result = []

        for i in range(min_lenght):
            temp = df[df["index"] == i]
            if len(temp) > 0 :
                sn1 = temp[carbons_sn1]
                sn1 = sn1.mean()
                sn1 = 1.5 * sn1 - 0.5
                sn1 = sn1.abs()
                sn1 = sn1.mean()
                final = sn1
                if carbons_sn2:
                    sn2 = temp[carbons_sn2]
                    sn2 = sn2.mean()
                    sn2 = 1.5 * sn2 - 0.5
                    sn2 = sn2.abs()
                    sn2 = sn2.mean()

                    final = (final + sn2)*0.5
                result.append([i,final])
            else:
                result.append([i,0])
        result = np.array(result)
        result = result[:,1]
        return result


    # Method to average vector to pseudovector program
    @staticmethod
    def average_vector(data, min_lenght):
        columns = ["index", "x", "y", "z"] # Data expected is an np array with columns ["index", "x", "y", "z"]

        df = pd.DataFrame(data, columns = columns)
        result = []

        for i in range(min_lenght):
            temp = df[df["index"] == i]
            if len(temp) > 0 :
                bin_vect = temp[columns[1:]]
                bin_vect = bin_vect.mean()
                result.append(bin_vect.to_list())
            else:
                result.append([np.nan, np.nan, np.nan])
        result = np.array(result)

        return result

    # Computes the average vector for each bin, sample are the raw x,y positions and weights are the vectors related to the head
    def pseudohistogram2D(self,sample1, weights, bins = 10, v_min = None, v_max = None):
        if v_min == None:
            v_min = np.min(sample1)
        if v_max == None:
            v_max = np.max(sample1)

        #print(v_min, v_max)
        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(v_min, v_max, bins +1)
            nbin[i] = len(edges[i]) + 1

        Ncount = (tuple(np.searchsorted(edges[i], sample1[:,i], side = "right") for i in range(2)))

        for i in range(2):
            on_edge = (sample1[:,i] == edges[i][-1])
            Ncount[i][on_edge] -= 1


        xy = np.ravel_multi_index(Ncount, nbin)
        xy_test = xy.reshape(-1,1)

        xy_test = np.concatenate((xy_test, weights), axis = 1)
        hist = self.average_vector(xy_test, nbin.prod())
        nbin = (nbin[0], nbin[1], 3)
        hist = hist.reshape(nbin)
        hist = hist.astype(float, casting = "safe")
        core = 2*(slice(1,-1),)
        hist = hist[core]

        return hist, edges

    # Computes teh histogram of the average order parameters in each bin
    def histogram2D(self,sample1, weights, n_chain, bins = 10, v_min = None, v_max = None):
        if v_min == None:
            v_min = np.min(sample1)
        if v_max == None:
            v_max = np.max(sample1)

        #print(v_min, v_max)
        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(v_min, v_max, bins +1)
            nbin[i] = len(edges[i]) + 1

        Ncount = (tuple(np.searchsorted(edges[i], sample1[:,i], side = "right") for i in range(2)))

        for i in range(2):
            on_edge = (sample1[:,i] == edges[i][-1])
            Ncount[i][on_edge] -= 1


        xy = np.ravel_multi_index(Ncount, nbin)
        xy_test = xy.reshape(-1,1)

        xy_test = np.concatenate((xy_test, weights), axis = 1)

        hist = self.count_order(xy_test, nbin.prod(), n_chain)
        hist = hist.reshape(nbin)
        hist = hist.astype(float, casting = "safe")
        core = 2*(slice(1,-1),)
        hist = hist[core]

        return hist, edges



    def order_histogram(self, lipid, layer, n_grid,
                        n_chain,
                        v_min = None,
                        v_max = None,
                        all_head = None,
                        start = None,
                        final = None,
                        step = 1):

        if all_head == None:
            all_head = self.all_head
        if start == None:
            start = self.start
        if final == None:
            final = self.final

        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "

        try:
            n_chain1 = n_chain[0]
            n_chain2 = n_chain[1]
        except:
            n_chain1 = n_chain
            n_chain2 = 0

        matrix = [] # this will store a matrix of the shape (2+n_chain,
        for ts in self.u.trajectory[start:final:step]:
            z = all_head.positions[:,2]
            z_mean = z.mean() # get middel of the membrane
            #Pick atoms in the layer
            if layer == "both":
                layer = self.u.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")
            else:
                layer = self.u.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {z_mean})")
            #print("Info:", all_head.n_atoms, z_mean, layer.n_atoms)

            only_p = layer.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            angles_sn1 = self.individual_order_sn1(layer, lipid, n_chain1)
            angles_sn1 = angles_sn1.T

            #print(angles_sn1.T.shape, positions.shape)
            #print(angles_sn1.shape, positions.shape)
            to_write = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                angles_sn2 = self.individual_order_sn2(layer, lipid, n_chain2)
                angles_sn2 = angles_sn2.T
                to_write = np.concatenate([to_write, angles_sn2], axis = 1)

            matrix.append(to_write) # Expect dim (n_lipids, 2+n_chain1+n_chain2)
            #print("Frame:",to_write.shape)

        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)


        H, edges = self.histogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = n_grid, v_min = v_min, v_max = v_max)
        H = np.rot90(H)
        H[H==0] = np.nan

        return H, edges

    @staticmethod
    def build_resname(resnames_list):
        resnames_list = list(resnames_list)
        string = " (resname " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or resname " + resname

        string = string + ") "
        return string

    @staticmethod
    def build_name(resnames_list):
        string = " (name " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or name " + resname

        string = string + ") "
        return string




top = "membrane.gro"
traj = "membrane.xtc"
membrane = twod_analysis(top, traj)
u = membrane.memb

print(membrane.lipid_list)







