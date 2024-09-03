import numpy as np
import qiskit as qk
from typing import List,Tuple,Union
from collections import Counter
import itertools
import networkx as nx
import copy
from enum import Enum
import random
import re
import gen_qsim_input_from_qasm


#grouping for QAOA (grouping_cascades) and for VQE state preparation (grouping_without_ordering)

class AllowedValues(Enum):
    ZERO = 0
    ONE = 1
    TWO = 2

class grouping_cascades():
    """attempts to group commuting gates into cascades which can then be cut"""

    def __init__(self, edges:List[Tuple[int]], q:int, cut_loc:int):
        self.edges = edges
        self.q = q
        self.cut_loc = cut_loc
        self.dct_weights=None #only necessary if weights are added to the maxcut

    def count_cuts(self):
        """counts the separate cuts and cuts with blocks from the graphs dict (only QAOA)"""
        count_blocks = 0 #number of blocks 
        count_remaining_sep = 0 #if blocks are done, how many separate cuts are left?
        count_sep = 0 #how many separate cuts are present if no blocks were done

        for key in self.edges_cascade_lim.keys():
            if key=="remainder":
                for tup in self.edges_cascade_lim[key]:
                    i=tup[0]
                    j=tup[1]
                    if (i<=self.cut_loc and j>self.cut_loc) or (j<=self.cut_loc and i>self.cut_loc):
                        count_sep += 1
                        count_remaining_sep += 1
            else:
                num_block = len(self.edges_cascade_lim[key])
                if  num_block == 1:
                    count_sep += 1
                    count_remaining_sep += 1
                else:
                    count_sep += num_block
                    count_blocks += 1

        return count_blocks, count_remaining_sep, count_sep


    def perform_grouping(self, perm_thres:int = 0):
        """
        performs the grouping procedure, depending on the location of the cut. cut_loc determines the qubit after which the cut is performed
        default value of perm_thres is 0 as in practice the ordering does not make a difference
        """

        lst_dicts = []

        #remove the edges which are not on the cut
        edges_reduced = []
        edges_remainder = []
        for edge in self.edges:
            if (edge[0] <= self.cut_loc and self.cut_loc < edge[1]) or (edge[1] <= self.cut_loc and self.cut_loc < edge[0]):
                edges_reduced.append(edge)
            else:
                edges_remainder.append(edge)
        edges_reduced_copy = edges_reduced.copy()
        cut_edges = len(edges_reduced_copy)

        #flatten the edges and count which indices appear multiple times
        edges_flattened = [item for sublist in edges_reduced for item in sublist]
        counter = Counter(edges_flattened)
        ordered_keys = [key for key, _ in counter.most_common()] #order the keys in decreasing order accorrding their appearance

        count_perm = 0
        for keys_list in itertools.permutations(ordered_keys[::-1]):
            keys_list = keys_list[::-1]
            edges_reduced = edges_reduced_copy.copy()
            #go through the indices which appear multiple times and check whether they can form a cascade
            dct = {}
            for vertex in keys_list:
                #gather together as many edges which are allowed to form a cascade into the dct.
                temp_cascade_candidate = [edge for edge in edges_reduced if vertex in edge] #!pay attention to  double elements, the elements in temp_Cascade_candidate must be removed from edges_reduced
                
                #split this temporary list into a list were all uncommon qubits are smaller than `vertex` and in one where all are larger
                cascade_candidate_above = []
                cascade_candidate_below = []
                for el in temp_cascade_candidate:
                    idx_vertex = list(el).index(vertex)
                    idx_other = (idx_vertex+1)%2
                    if  el[idx_other] > vertex:
                        cascade_candidate_above.append(el)
                    else:
                        cascade_candidate_below.append(el)

                if len(temp_cascade_candidate)==0 and len(edges_reduced)!=0: #if none is retrieved but edges_reduced still nonzero, we have to include those
                    temp_cascade_candidate = edges_reduced
                    for i,el in enumerate(temp_cascade_candidate):
                        dct.update({str(vertex)+"s"+str(i) : [el]})
                        edges_reduced.remove(tuple(el))
                        
                if len(cascade_candidate_above)!=0 and len(cascade_candidate_below)==0:
                    dct.update({str(vertex):cascade_candidate_above})
                    #remove the elements from edges_reduced
                    for tup in cascade_candidate_above:
                        edges_reduced.remove(tup)

                elif len(cascade_candidate_above)==0 and len(cascade_candidate_below)!=0:
                    dct.update({str(vertex):cascade_candidate_below})
                    #remove the elements from edges_reduced
                    for tup in cascade_candidate_below:
                        edges_reduced.remove(tup)

                elif len(cascade_candidate_above)!=0 and len(cascade_candidate_below)!=0:
                    #if cascades are possible in both directions add both to the dict
                    dct.update({str(vertex)+"a":cascade_candidate_above})
                    dct.update({str(vertex)+"b":cascade_candidate_below})
                    #remove the elements from edges_reduced
                    for tup in cascade_candidate_below:
                        edges_reduced.remove(tup)
                    for tup in cascade_candidate_above:
                        edges_reduced.remove(tup)

                elif len(cascade_candidate_above)==0 and len(cascade_candidate_below)==0:
                    pass
                
                else:
                    raise ValueError("Something is wrong with the internal case distinction")

                if len(edges_reduced)==0:
                    break  

            lst_dicts.append(dct)
            count_perm += 1

            if count_perm >= perm_thres:
                break

        #choose elemnt with fewest entries in lst_dicts as dct
        dct = min(lst_dicts, key=lambda d: len(d))

        flattened_final_dct = [tuple for sublist in dct.values() for tuple in sublist]

        edge_counter = Counter(flattened_final_dct)

        assert len(set(flattened_final_dct)) == len(flattened_final_dct), "There are duplicates in the final dct"
        
        assert len(flattened_final_dct) == cut_edges, f"The number of tuples in the final result differ from the original number of cut edges: final={len(flattened_final_dct)}, initial={cut_edges}"
        #after iterating through all possibilities, check which dict in the list has the least keys. this is the desired partitioning.
        
        #add edges which are not in the cut
        dct.update({"remainder" : edges_remainder})
        
        self.edges_cascades = dct

    def limit_active_qubits_per_part(self, lim_active=5):
        """
        qsimh is limited regarding the number of qubits in a partition of a block. the default value for lim_active=5 guarantees that the amplitudes are correct. larger values may cause issues with the results 
        hence the cascades must be chosen such that this limit is not exceeded. the block/cascades are chosen w.r.t the cut location, the rest is pushed to the remainder
        """

        edges_cascade = copy.deepcopy(self.edges_cascades)
        edges_cascade_lim = {"remainder" : copy.deepcopy(edges_cascade["remainder"])}

        for key, lst in edges_cascade.items():
            if key != "remainder":
                temp_lst_within = []
                temp_lst_out = []

                indices_top = set()
                indices_bottom = set()
                for el in lst:
                    qubit0 = el[0]
                    qubit1 = el[1]
                    if qubit0 <= self.cut_loc:
                        indices_top.add(qubit0)
                    else:
                        indices_bottom.add(qubit0)
                    if qubit1 <= self.cut_loc:
                        indices_top.add(qubit1)
                    else:
                        indices_bottom.add(qubit1)

                indices_to_remove = [] #choose indices randomly if elements have to be moved out of the cascade
                if len(list(indices_top)) <= lim_active and len(list(indices_bottom)) <= lim_active:
                    temp_lst_within = lst
                elif len(list(indices_top)) > lim_active and len(list(indices_bottom)) <= lim_active:
                    to_remove = len(list(indices_top))-lim_active
                    while len(indices_to_remove) < to_remove:
                        idx = random.randint(0, len(lst) - 1)
                        if idx not in indices_to_remove:
                            indices_to_remove.append(idx)
                elif len(list(indices_top)) <= lim_active and len(list(indices_bottom)) > lim_active:
                    to_remove = len(list(indices_bottom))-lim_active
                    while len(indices_to_remove) < to_remove:
                        idx = random.randint(0, len(lst) - 1)
                        if idx not in indices_to_remove:
                            indices_to_remove.append(idx)
                elif len(list(indices_top)) > lim_active and len(list(indices_bottom)) > lim_active:
                    raise NotImplementedError("This case is very unlikely as cascades have a single active qubit in onepart and only multiple in the other")
                if len(indices_to_remove) != 0:
                    for i in range(len(lst)):
                        if i in indices_to_remove:
                            temp_lst_out.append(lst[i])
                        else:
                            temp_lst_within.append(lst[i])
                
                if len(temp_lst_within)!=0:
                    edges_cascade_lim.update({key : temp_lst_within})
                for el in temp_lst_out:
                    edges_cascade_lim["remainder"].append(el)
        self.edges_cascade_lim = edges_cascade_lim       


    def create_qaoa_qsim(self, angles:List[float], cascade:AllowedValues, path:str, weights:List[float]=[]):
        """
        creates qaoa circuit for qsim. It writes a file and stores it into the `path`
        angles contains the theta and gamma angles for the problem and mixer layer. for each layer another value, 
        hence the number of layers can be deduced from this list's size
        angles = [theta1, theta2, gamma1, gamma2], theta is problem layer, gamma mixer

        the bool `cascade` decides whether the ordered gates are utilized (cascades=True) or the original random order.
        if `cascade`=True also the brackets are included such that the blocks are used in qsimh
        update: true = 0, false = 1, reduced blocksize = 2

        The optional input `weights` is only relevant if you want to do weighted Max Cut. This adds the weights according to the same order as the edges of the given graph

        .. danger::

            if weights is nontrivial, you need to run cascades=1 first, before you run cascades=2 in order to receive the correct assignments of the weights
        
        """
        
        if len(weights)!=0:
            assert len(self.edges) == len(weights), "Your weights for weighted MaxCut do not fit your Graph."
            bool_weights = True
        else:
            bool_weights = False

        if (not self.dct_weights) and bool_weights==True and cascade!=1:
            raise ValueError("You need to run `cascade=1`, i.e. no block first before you run any version with blocks. This is because the weight assignment must be determined in the first run.")
        elif (not self.dct_weights) and bool_weights==True and cascade==1:
            self.dct_weights={} #assign dict to this var


        assert len(angles)%2==0, "The number of elements in `angles` must be even"
        layers = len(angles)//2

        with open(path, 'w') as file:
            file.write(f"{self.q}\n")
            
            #hadamards in the beginning
            for i in range(self.q):
                file.write(f"0 h {i}\n")

            gate_counter = 1
            for layer in range(layers):
                theta = angles[layer]
                gamma = angles[layer+layers]
                #problem layer
                if cascade == 1:
                    if bool_weights == False:
                        for edge in self.edges:
                            file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta}\n")
                            gate_counter += 1
                    else: #with weights
                        for edge, weight in zip(self.edges, weights):
                            file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta*weight}\n")
                            self.dct_weights.update({f"{edge[0]} {edge[1]}":weight}) #need these entries to correctly assign the weights for cascade=2,0
                            gate_counter += 1

                elif cascade == 0 or cascade == 2:
                    if cascade == 0:
                        casc_dict = self.edges_cascades
                    elif cascade == 2:
                        casc_dict = self.edges_cascade_lim

                    weights_counter = 0 #only necessary if weights presnet
                    for key, casc_list in casc_dict.items():
                        for l, edge in enumerate(casc_list):
                            if len(casc_list)>1 and key!="remainder": # the remainder gates should not be in brackets because they are not grouped
                                if l==0:
                                    file.write("[")
                                if bool_weights == False:
                                    file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta}")
                                else:#with weights
                                    #search correct weight in dct_weights
                                    try:
                                        weight_current = self.dct_weights[f"{edge[0]} {edge[1]}"]
                                    except KeyError:
                                        raise KeyError(f"Error: Weight not found for edge ({edge[0]}, {edge[1]})")
                                    file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta*weight_current}")
                                    weights_counter +=1
                                if l == len(casc_list)-1:
                                    file.write("]\n")
                                else:
                                    file.write("\n")
                            else: #no bracket if only a single element in casc_list
                                if bool_weights==False:
                                    file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta}\n")
                                else:#with weights
                                    #search correct weight in dct_weights
                                    try:
                                        weight_current = self.dct_weights[f"{edge[0]} {edge[1]}"]
                                    except KeyError:
                                        raise KeyError(f"Error: Weight not found for edge ({edge[0]}, {edge[1]})")
                                    file.write(f"{gate_counter} rzz {edge[0]} {edge[1]} {theta*weight_current}\n")
                                    weights_counter +=1
                            gate_counter += 1


                #mixer layer
                for i in range(self.q):
                    file.write(f"{gate_counter} rx {i} {gamma}\n") 
                gate_counter += 1

            print(f"Circuit written into {path}")

    def create_qaoa_qiskit(self, angles:List[float], cascade:AllowedValues):
        """
        creates qaoa circuit for qiskit. 
        angles contains the theta and gamma angles for the problem and mixer layer. for each layer another value, 
        hence the number of layers can be deduced from this list's size
        angles = [theta1, theta2, gamma1, gamma2], theta is problem layer, gamma mixer

        the bool `cascade` decides whether the ordered gates are utilized (cascades=True) or the original random order.
        update: true = 0, false = 1, reduced blocksize = 2
        """
        assert len(angles)%2==0, "The number of elements in `angles` must be even"
        layers = len(angles)//2

        #initialize circuit object
        circ = qk.QuantumCircuit(self.q)

        #hadamards in the beginning
        for i in range(self.q):
            circ.h(i)

        for layer in range(layers):
            theta = angles[layer]
            gamma = angles[layer+layers]
            #problem layer
            if cascade == 1:
                for edge in self.edges:
                    circ.rzz(theta, edge[0], edge[1])
            elif cascade == 0:
                for casc_list in self.edges_cascades.values():
                    for edge in casc_list:
                        circ.rzz(theta, edge[0], edge[1])
            elif cascade == 2:
                for casc_list in self.edges_cascade_lim.values():
                    for edge in casc_list:
                        circ.rzz(theta, edge[0], edge[1])
        
            #mixer layer
            for i in range(self.q):
                circ.rx(gamma, i)
        
        return circ
    
