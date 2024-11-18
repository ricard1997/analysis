import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simps
import time
from matplotlib.patches import Patch
import nglview as nv




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
                v_min = 0,
                v_max = 180,
                add_radii = False,
                verbose = False,
            ):



        # Read trajectory depending if tpr is provided or not
        if tpr:
            self.u = mda.Universe(tpr, traj)
        else:
            self.u = mda.Universe(top, traj)

        # Select elements in the membrane (in principle only lipids)
        if not lipid_list: # Select only elements of the membrane
            self.memb = self.u.select_atoms("all and not protein and not (resname URA or resname GUA or resname ADE or resname CYT)")
            self.lipid_list = set(self.memb.residues.resnames)
        else:
            self.memb = self.u.select_atoms(f"{self.build_resname(list(lipid_list))}")


        # Set radius sizes of different elements
        self.radii_dict = 0
        if add_radii:
            self.radii_dict = {"H": 0.7,
                            "N": 1.85,
                            "C": 2.06,
                            "P": 2.15,
                            "O": 1.65,
                            }


            # Add radii as a topology attribute for Mdanalysis
            string_array = self.memb.elements
            radii_array = np.array([self.radii_dict[element] for element in string_array])
            self.u.add_TopologyAttr("radii")
            self.memb.radii = radii_array

            # May use to build a different way to select polar atoms
            #polar_motif = "N HN1 HN2 HN3 C12 H12A C13 O13A O13B C11 H11A H11B"
            #polar_PS = "N HN1 HN2 HN3 C12 H12A C13 O13A O13B C11 H11A H11B"
            #polar_PI = "C12 H2 O2 HO2 C13 H3 O3 HO3 C14 H4 O4 HO4 C15 H5 O5 HO5 C16 H6 O6 HO6 C11 H1"
            #polar_PA = "H12 "
            #polar_PC = "N C12 C13 C14 C15 H12A H12B H13A H13B H13C H14A H14B H14C H15A H15B H15C C11 H11A H11B"
            #polar_PE = "N HN1 HN2 HN3 C12 H12A H12B C11 H11A H11B"
            #polar_CHL = "O3 H3'"


            #polar_chains = [polar_motif, polar_PS, polar_PI, polar_PA, polar_PC, polar_PE]
            #polar_atoms = [chain.split() for chain in polar_atoms]
            #dspc = self.memb.select_atoms("(resname DSPC and not (name C3* or name H*X or name H*Y or name C2* or name H*R or name H*S)) or (resname DSPC and(name C3 or name HX or name HY or name C2 or name HR or name HS))")
            #print(set(dspc.atoms.names))







        self.v_min = v_min
        self.v_max = v_max
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
            self.non_polar_dict = {}
            self.first_lipids = {}
            self.non_polar_visualize = {}
            for lipid in self.lipid_list:
                first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
                actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
                actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
                actual_sn1 = actual_sn1.names
                actual_sn2 = actual_sn2.names
                self.chain_info[lipid] = [len(actual_sn1) - 2, len(actual_sn2) - 2]
                self.first_lipids[lipid] = first_lipid

                if lipid == "CHL1":
                    non_polar = self.memb.select_atoms(f"resid {first_lipid} and not (name O3 or name H3')")
                    all_lip = self.memb.select_atoms(f"resid {first_lipid}")

                else:
                    non_polar = self.memb.select_atoms(f"resid {first_lipid} and (name *C3* or name H*Y or name H*X or name H*Z  or name *C2* or name H*R or name H*S or name H*T) and not (name C3 or name C31 or name HY or name HX or name HZ  or name C2 or name C21 or name HR or name HS or name HT)")
                    all_lip = self.memb.select_atoms(f"resid {first_lipid}")
                self.non_polar_dict[lipid] = list(non_polar.names)

                self.non_polar_visualize[lipid] = [all_lip, non_polar]



        self.all_head = self.u.select_atoms(self.build_resname(self.lipid_list) + " and name P")
        self.start = 0
        self.final = 100
        self.step = 1

        if verbose:
            print(f"This system contains the following lipids : {self.lipid_list}\n\n")
            print(f"The chain lenght is : \n{self.print_dict(self.chain_info)}\n")
            print(f"We will use the following heads and charges for the following lipids. If the lipid is not here we will use P as head as default \n{self.print_dict(self.working_lip)}\n")
            print("Note: To compute the middle of the membrane we use only P heads\n\n")
            print(f"The default start frame is {self.start}, final {self.final}, step {self.step}\n\n")

    # Method to print dictionaries
    @staticmethod
    def print_dict(dict):
        string = ""
        for key in dict.keys():
            string += f"{key} : {dict[key]}\n"
        return string


    def visualize_polarity(self, lipids = "all"):
        """This function is used to visualize what atoms are considered in polarity

        Args:
            lipids (str or list of str, optional): Lipids to show polarity. Defaults to "all".
        """
        aspect_ratio = [1, 1, 1]

        # Get lipids to work
        if lipids == "all":
            lipids = self.lipid_list
        else:
            if isinstance(lipids, list):
                lipids = [lipids]


        #Guess bonds if needed
        try:
            self.u.bonds
        except:
            bonds = mda.topology.guessers.guess_bonds(self.u.atoms, self.u.atoms.positions)
            self.u.add_TopologyAttr("bonds", bonds)
            #print(mda.topology.guessers.guess_bonds(first_lipid, first_lipid.positions))

        # Make a 3D plot of the individual lipids and shows it
        if not self.non_polar_dict:
            print("Non polar atoms are not yet established, please set self.non_polar_dict first")
            return
        else:
            nplots = len(lipids)

            fig = plt.figure(figsize=(5 * nplots, 5))


            count = 0
            for lipid in lipids:
                first_lipid = self.memb.select_atoms(f"resid {self.first_lipids[lipid]} and resname {lipid}")
                lipid_pos = first_lipid.positions
                lipid_ats = first_lipid.names
                polarity = ["blue" if name in self.non_polar_dict[lipid] else "red" for name in lipid_ats]
                #print(polarity, lipid_ats)
                ax = fig.add_subplot(1,nplots, count+1, projection = "3d")
                ax.scatter(*lipid_pos.T, s = 10, c = polarity)
                dists = []
                mids = []
                for i in range(3):
                    value = np.max(lipid_pos[:,i]) - np.min(lipid_pos[:,i])
                    dists.append(value)
                    mids.append(np.min(lipid_pos[:,i] + value/2))
                dists = max(dists)
                dists = dists + 0.01*dists
                for atom in first_lipid:
                    bonds = atom.bonded_atoms
                    for bond in bonds:
                        vector = np.array([atom.position, bond.position])
                        plt.plot(*vector.T, color = "black")
                for i in range(len(lipid_ats)):
                    if lipid_ats[i] in self.non_polar_dict[lipid]:
                        ax.text(lipid_pos[i,0], lipid_pos[i,1],  lipid_pos[i,2], lipid_ats[i], color = "black", fontsize = 6)
                    else:
                        ax.text(lipid_pos[i,0], lipid_pos[i,1], lipid_pos[i,2], lipid_ats[i], color = "black", fontsize = 6)
                ax.set_box_aspect(aspect_ratio)
                ax.set_title(lipid)
                ax.set_xlim(mids[0]-dists/2, mids[0]+dists/2)
                ax.set_ylim(mids[1]-dists/2, mids[1]+dists/2)
                ax.set_zlim(mids[2]-dists/2, mids[2]+dists/2)
                count += 1
        plt.show()







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

        Code to loop over the number of carbons_summary_ in the lipid tail and get a list with the carbon and its
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
        print(lipid)

        # Loop over carbons
        for i in range(n_chain):
            # Define selections for H and C in the chain
            #print(f"Value of the chain {i} sn1")
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


            angles = self.get_individual(lista)
            #print(angles, "############################")
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

    @staticmethod
    def get_highest(data, min_lenght):
        """Code to get the highest value given two columns that are to be ordered in a 2D grid

        Args:
            data (ndarray(:,2)): Array with two columns (column1: map to a 2D grid, column2: values)
            min_lenght (int): Size of squares in the 2D grid

        Returns:
            ndarray(:,2): With the maximun of each grid square
        """
        columns = ["index", "weight"] # Data expected is an np array with columns ["index", "x", "y", "z"]
        df = pd.DataFrame(data, columns = columns)
        result = df.groupby('index', as_index=False)['weight'].max()
        result_dict = dict(zip(result['index'], result['weight']))
        hist = []
        for i in range(min_lenght):
            try:
                hist.append(result_dict[i])
            except:
                hist.append(np.nan)
        return np.array(hist)

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

    # Computes  and return the indexes of the data if where arranegd in a 2d histogram
    def get_indexes(self,
                    data,
                    bins = 10,
                    v_min = None,
                    v_max = None,
                    matrix_height = False):

        if v_min == None:
            v_min = self.v_min
        if v_max == None:
            v_max = self.v_max


        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(v_min, v_max, bins +1)
            nbin[i] = len(edges[i]) + 1


        if not matrix_height:

            indexes = (tuple(np.searchsorted(edges[i], data[:,i], side = "right") for i in range(2)))
        else:
            indexes = (tuple(np.searchsorted(edges[i], data[:,i], side = "right") for i in range(2)))

            xy = np.ravel_multi_index(indexes, nbin)
            xy_test = xy.reshape(-1,1)

            #print("last shape", data[:,2].reshape(-1,1).shape, xy_test.shape)
            xy_test = np.concatenate((xy_test, data[:,2].reshape(-1,1)), axis = 1)
            hist = self.get_highest(xy_test, nbin.prod())


            hist = hist.reshape(nbin)
            hist = hist.astype(float, casting = "safe")
            hist[np.isnan(hist)] = 0
            #core = 2*(slice(1,-1),)
            #hist = hist[core]
            #print("here", hist[20,:])
            return indexes, hist

        #for i in range(2):
        #    on_edge = (data[:,i] == edges[i][-1])
        #    Ncount[i][on_edge] -= 1
        #print(np.min(Ncount[0]), np.max(Ncount[0]), np.min(Ncount[1]), np.max(Ncount[1]))
        #print(len(edges[0]), "edges len")




        return indexes



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
            n_chain1 = n_chain[1]
            n_chain2 = n_chain[0]
        except:
            n_chain1 = n_chain
            n_chain2 = 0

        matrix = [] # this will store a matrix of the shape (2+n_chain,
        for ts in self.u.trajectory[start:final:step]:
            z = all_head.positions[:,2]
            z_mean = z.mean() # get middel of the membrane

            #Pick atoms in the layer
            if layer == "both":
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")
            else:
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {z_mean})")
            #print("Info:", all_head.n_atoms, z_mean, layer.n_atoms)

            only_p = layer_at.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            angles_sn1 = self.individual_order_sn1(layer_at, lipid, n_chain1)
            angles_sn1 = angles_sn1.T

            print(angles_sn1.T.shape, positions.shape)
            print(angles_sn1.shape, positions.shape)
            to_write = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                angles_sn2 = self.individual_order_sn2(layer_at, lipid, n_chain2)
                angles_sn2 = angles_sn2.T
                to_write = np.concatenate([to_write, angles_sn2], axis = 1)

            matrix.append(to_write) # Expect dim (n_lipids, 2+n_chain1+n_chain2)
            #print("Frame:",to_write.shape)

        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)
        v_min = self.v_min
        v_max = self.v_max

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

    def build_resname_head(self,resnames_list):
        resnames_list = list(resnames_list)
        string = f"( (resname {resnames_list[0]}  and name {self.working_lip[resnames_list[0]]['head']}) "

        for resname in resnames_list[1:]:
            string += f" or (resname {resnames_list[0]}  and name {self.working_lip[resnames_list[0]]['head']}) "

        string +=  " ) "
        return string

    @staticmethod
    def build_name(resnames_list):
        string = " (name " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or name " + resname

        string = string + ") "
        return string


    def all_lip_order(self, layer, nbins,
                        v_min = None,
                        v_max = None,
                        all_head = None,
                        start = None,
                        final = None,
                        step = 1,
                        plot = False):
        """Find the 2D order parameters for all lipids

        Args:
            layer (str): Layer, can be top, bot, both
            nbins (int): number of bins
            v_min (float, optional): min value for the grid. Defaults to None.
            v_max (float, optional): max value for the grid. Defaults to None.
            all_head (_type_, optional): _description_. Defaults to None.
            start (int, optional): start frame. Defaults to None.
            final (int, optional): final frame. Defaults to None.
            step (int, optional): step. Defaults to 1.
            plot (Bool, optional): plot the resulting matrix. Defaults to False

        Returns:
            ndarray(n,n), ndarray(n+1): matrix containing the 2d order, edges of the matrix
        """

        lipid_list = list(self.lipid_list)
        lipid_list.remove("CHL1")
        lipids = self.chain_info

        v_min = self.v_min
        v_max = self.v_max

        matrices = []
        for key in lipid_list:
            print(key)
            H, edges = self.order_histogram(key, layer, nbins, lipids[key],v_min = v_min,
                        v_max = v_max,
                        all_head = all_head,
                        start = start,
                        final = final,
                        step = step)

            matrices.append(H)
        matrices = np.array(matrices)
        matrices = np.nanmean(matrices, axis = 0)

        if plot:
            plt.close()
            plt.imshow(matrices[1:-1,1:-1] ,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
            plt.colorbar(cmap = "Spectral")
            plt.savefig(f"all_lip1_{layer}.png")
            plt.close()

        return matrices, edges

    def surface(self,
                start = None,
                final = None,
                step = None,
                lipid = "DSPC",
                layer = 'top',
                filename = None, include_charge = False):
        """Code to loop over the trajectory and print [x,y,z(referenced to zmean), charge] in a file.

        Args:
            start (int, optional): Start Frame. Defaults to None.
            final (int, optional): Final frame. Defaults to None.
            step (int, optional): Frames to skip. Defaults to None.
            lipid (str, optional): Lipid to work. Defaults to "DSPC".
            layer (str, optional): Layer to work. Defaults to 'top'.
            filename (str, optional): filename to write data. Defaults to None.
            include_charge (bool, optional): Include or not charge. Defaults to False.

        Returns:
            _type_: _description_
        """
        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step
        if filename == None:
            filename = f"{lipid}_{layer}_{start}_{final}.dat"
        lipid_list = self.lipid_list


        print("######### Running surface function ########## ")
        print(f"We will compute the surface files for a membrane with there lipids {lipid_list}")
        print(f"Currently working on: {lipid}")
        print(f"Layer: {layer}")
        print(f"Writing under the name of {filename}")


        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()

            # Selects the lipid head and of the working lipid
            if layer == "both":
                selection_string = f"(resname {lipid} and name {self.working_lip[lipid]['head']})"
            else:
                selection_string = f"(resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {str(mean_z)}"

            # Find the positions of the P atoms
            atoms = self.u.select_atoms(selection_string)

            ### Get positions
            atom_pos = atoms.positions


            ### Get resids
            atom_resid = atoms.resids
            atom_resid = atom_resid[:,np.newaxis]

            atom_pos = np.concatenate((atom_pos, atom_resid), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2]-mean_z)

            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        df_data = pd.DataFrame(pos_data, columns = ["x", "y", "z", "id"])
        df_data["id"] = df_data["id"].astype(int)
        #if include_charge:
        #    df_data["charge"] = self.charge_li[lipid]
        #    df_data.to_csv(f"pd_{filename}", index = False)
        df_data.to_csv(f"pd_{filename}", index = False)

        return df_data   # Maybe have to change, it does not make sense to return thi


    def height_matrix(self, lipids, layer,start = None, final = None, step = None, nbins = 50, clean = True):
        """Code to divide the space in a 2D grid and compute the height referenced to zmean

        Args:
            lipids (list(str)): Lipids to include in the height analysis
            layer (str): Working layer for thickness
            start (int, optional): Frame to start analysis. Defaults to None.
            final (int, optional): Final frame for the analysis. Defaults to None.
            step (int, optional): Steps to skip. Defaults to None.
            nbins (int, optional): Number of bins to divide the grid space. Defaults to 50.
            clean (bool, optional): Decide if rerun and overwrite surface generated files. Defaults to True.

        Returns:
            ndarray(nbins,nbins): Retun a matrix with the height information
        """

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step


        print(f"Computing matrix for {layer} in frames {start}-{final}")
        data = []
        for lipid in lipids:
            filename = f"{lipid}_{layer}_{start}-{final}.dat"
            if not clean:
                try:
                    df_data = pd.read_csv(f"pd_{filename}")
                except:
                    self.surface(lipid = lipid, layer = layer, filename = filename, include_charge = True, start = start, final = final)
                    df_data = pd.read_csv(f"pd_{filename}")
            else:
                self.surface(lipid = lipid, layer = layer, filename = filename, include_charge = True, start = start, final = final)
                df_data = pd.read_csv(f"pd_{filename}")

            data.append(df_data)

        data = pd.concat(data, axis = 0)
        print(data)
        xmin = data["x"].min()
        xmax = data["x"].max()
        ymin = data["y"].min()
        ymax = data["y"].max()

        xmin = self.v_min
        xmax = self.v_max
        ymin = self.v_min
        ymax = self.v_max

        H_height, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["z"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_height/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        np.savetxt(f'Height_{layer}_{start}_{final}.dat', H_avg, fmt = '%.2f')
        np.savetxt(f"edges_{layer}_{start}_{final}.dat", x_edges, fmt = "%.2f")
        return H_avg, x_edges



    def thickness(self, nbins, start = 0, final=-1, step = 1):
        """Find the thichness mapped in a 2d grid
        Args:
            nbins (int): number of bins for thickness
            start (int, optional): Start frame. Defaults to 0.
            final (int, optional): Final frame. Defaults to -1.
            step (int, optional): Step frame. Defaults to 1.

        Returns:
            np.array, np.array: Matrix with the thickness, edeges for the matrix
        """
        lipids = list(self.lipid_list)
        lipids.remove("CHL1")
        matrix_up, edges = self.height_matrix(lipids,
                        "top",
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)
        matrix_bot, edges = self.height_matrix(lipids,
                        "bot",
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)
        mat_thickness = np.sum(np.array([matrix_up, matrix_bot]),axis = 0)
        mat_thickness[mat_thickness == 0] = np.nan
        print(mat_thickness,mat_thickness.shape,matrix_bot.shape,[edges[0], edges[-1], edges[0], edges[-1]])
        plt.close()
        plt.imshow(mat_thickness[1:-1,1:-1] ,cmap = "Spectral", extent = [edges[0], edges[-1], edges[0], edges[-1]])
        plt.colorbar(cmap = "Spectral")
        plt.savefig(f"all_lip_thick_.png")
        plt.close()
        return mat_thickness, edges



    def guess_minmax_space(self):
        """Check the minimun and max position in x,y

        Returns:
            float, float: minimun and max position in x,y
        """
        positions = self.memb.positions[:,2]
        vmin = np.min(positions)
        vmax = np.max(positions)
        return vmin,vmax

    @staticmethod
    def create_circle_array(grid_size, radius_A, center=None):
        """Create a small matrix with a inner circle of the size of radius_A

        Args:
            grid_size (float): define the grid size to create the optimun grid (Amstrongs)
            radius_A (float): define the radius for the matrix (amstrongs)
            center (bool, optional): Bool to set the center of the circle. Defaults to None.

        Returns:
            np.array : array with a circle of ones
        """

        n = int(2*radius_A /grid_size) + 1
        if n % 2 == 0:
            n += 1

        # Create an n x n array initialized to 1
        array = np.zeros((n, n))

        # Default the center to the middle of the array if not provided
        if center is None:
            center = (n // 2, n // 2)

        # Generate a grid of indices
        y, x = np.ogrid[:n, :n]

        #   Calculate the distance from each grid point to the center
        distance_from_center = (x - center[1])**2 + (y - center[0])**2


        # Set values to 2 within the circle of the given radius
        array[distance_from_center <= (radius_A/grid_size)**2] = 1

        return array

    @staticmethod
    def add_small_matrix(big_matrix, small_matrix, center_i, center_j):
        """Add smmall matrix to a big matrix

        Args:
            big_matrix (ndarray(n,n)): big matrix where a small  matrix will be added
            small_matrix (ndarray(m,m)): small matrix to be added
            center_i (int): i coordinate
            center_j (int): j coordunate

        Returns:
            big_matrix (ndarray(n,n)): big matrix modified
        """

    # Calculate the top-left corner of the submatrix in big_matrix
        start_i = center_i - small_matrix.shape[0] // 2
        start_j = center_j - small_matrix.shape[1] // 2
        end_i = start_i + small_matrix.shape[0]
        end_j = start_j + small_matrix.shape[1]

    # Handle boundaries to ensure indices stay within big_matrix
        big_start_i = max(0, start_i)
        big_start_j = max(0, start_j)
        big_end_i = min(big_matrix.shape[0], end_i)
        big_end_j = min(big_matrix.shape[1], end_j)

    # Calculate the overlapping region for small_matrix
        small_start_i = big_start_i - start_i
        small_start_j = big_start_j - start_j
        small_end_i = small_start_i + (big_end_i - big_start_i)
        small_end_j = small_start_j + (big_end_j - big_start_j)

    # Add the overlapping region of small_matrix to the big_matrix
        big_matrix[big_start_i:big_end_i, big_start_j:big_end_j] += small_matrix[small_start_i:small_end_i, small_start_j:small_end_j]

        return big_matrix

    def add_deffects(self,
                    matrix,
                    indexes,
                    elements,
                    names,
                    lipid,
                    mat_radii_dict):
        """ Code to easily add deffects in the 2d matrix
        Args:
            matrix (ndarray(n,n)): Matrix where the deffects are going to be added
            indexes (ndarray(i,j)): List of indexes i,j in the matrix where the deffects should be added
            elements (list): type of element (needed to put the right radious)
            names (list): names of the atoms (needed to map hydrophobic and not hydrophobic atoms)
            lipid (str): lipid name
            mat_radii_dict (dict): dictionary with the radii

        Returns:
            ndarray: matrix matrix filled with the deffects
        """
        matrix = matrix

        for i in range(len(indexes[0])):

            small_matrix = mat_radii_dict[elements[i]]

            if names[i] in self.non_polar_dict[lipid]:
                small_matrix = small_matrix * 0.0001
            self.add_small_matrix(matrix, small_matrix, indexes[0][i], indexes[1][i])
        return matrix


    def packing_defects(self,
                        start = None,
                        final = None,
                        step = None,
                        layer = 'top',
                        nbins = 180,
                        height = False,
                        ):
        """_summary_

        Args:
            start (int, optional): Frame to start. Defaults to None.
            final (int, optional): Frame to finish. Defaults to None.
            step (int, optional): Frames to skip. Defaults to None.
            layer (str, optional): working layer (top/bot). Defaults to 'top'.
            nbins (int, optional): Number of bins of the xy grid. Defaults to 180.
            height (bool, optional): Store height matrix (To study deepness of th epacking defects). Defaults to False.

        Returns:
            ndarray : If height == Flase: matrix with pcking deffects
            ndarray, ndarray : If height === True: matrix with pcking deffects, amtrix with height information
        """




        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step

        vmin = self.v_min
        vmax = self.v_max
        grid_size = abs(vmin - vmax)/nbins



        lipid_list = list(self.lipid_list)
        non_polarity = self.non_polar_dict


        print("######### Running packing defects function ########## ")
        print(f"We will compute packing defects for a membrane with lipids {lipid_list}")





        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []

        mat_radii_dict = {}
        for atom in self.radii_dict.keys():
            mat_radii_dict[atom] = self.create_circle_array(grid_size, self.radii_dict[atom])


        matrix = np.zeros((nbins+2, nbins+2))
        if height:
            matrix_height = np.zeros((nbins+2, nbins+2))
        positions = all_p.positions[:,2]
        mean_z = positions.mean()

        for lipid in self.lipid_list:
            selection_string = f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {mean_z})"

            layer_at = self.memb.select_atoms(selection_string)


            pos_ats = layer_at.positions
            if not height:
                indexes = self.get_indexes(pos_ats[:,:2], nbins)
            else:
                indexes, matrix_temp = self.get_indexes(pos_ats, nbins, matrix_height = True)

                matrix_height = np.maximum(matrix_height.copy(), matrix_temp.copy())

            elements = layer_at.elements
            names = layer_at.names

            matrix = self.add_deffects(matrix, indexes,elements, names, lipid, mat_radii_dict)


        deffects = np.where(matrix < 1, matrix, np.nan)


        if height:
            matrix_height[matrix_height == 0 ] = np.nan
            return deffects, matrix_height
        return deffects












