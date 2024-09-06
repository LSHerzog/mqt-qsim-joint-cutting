import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
from typing import List, Tuple
import sys
sys.path.append("..")
import grouping_cascades

class WeightedMaxCut:
    """
    generates and stores circuits for both joint and separate cutting as well as no cutting for weightedMaxCut
    
    Parameters:
    -----------
    q : int
        number of qubits in the system
        
    seed : int
        choose seed in order to be able to regenerate the graphs
        
    cut_loc : int
        position of the cut. note that the qubit indices start with 0 and the cut is performed AFTER the given qubit
        
    sizes : List[int]
        Contains the partitioning of qubits/vertices in the graph. For instance sizes = [10,10] creates two partitions with 10 vertices/qubits each.
        Within the partitions vertices are connected with probability p_intra and between them with p_inter.
        Currently, one can use 2 or 3 partitions.
         
    p_intra : float
        probability to connect the vertices in the graph (within a partition)
        
    p_inter : float
        probability to connect the vertices in the graph (between the partitions)
        
    angles : List[float]
        Angles of the problem and mixer layers. From the length of this list, the number of QAOA layers is deduced. 
        For instance angles = [theta1, theta2, gamma1, gamma2], theta is problem layer, gamma the mixer layer.
        
    save_plot : bool, optional
        Saves the plot of the problem graph or not (default is False)
        
    verbose : bool, optional
        Decides about verbosity (default is True)
        
    weights : List[float], optional
        The weights for each edge in the problem graph. If initially empty, random weights are generated. Otherwise one can pass weights as input. (default is an empty list)
    """

    def __init__(self, 
                 q: int, 
                 seed: int, 
                 cut_loc: int, 
                 sizes: List[int], 
                 p_intra: float, 
                 p_inter: float,
                 angles: List[float],
                 save_plot: bool = False,
                 verbose : bool = True,
                 weights : List[float] = []):
        
        self.q: int = q
        self.seed: int = seed
        self.cut_loc: int = cut_loc
        self.sizes: List[int] = sizes
        self.p_intra: float = p_intra
        self.p_inter: float = p_inter
        self.angles = angles
        self.save_plot = save_plot
        self.verbose = verbose
        self.weights = weights

        if len(sizes)==2:
            self.p_matrix: List[List[float]] = [
                [self.p_intra, self.p_inter], 
                [self.p_inter, self.p_intra], 
            ]
        elif len(sizes)==3:
            self.p_matrix: List[List[float]] = [
                [self.p_intra, self.p_inter, self.p_inter], 
                [self.p_inter, self.p_intra, self.p_inter], 
                [self.p_inter, self.p_inter, self.p_intra]
            ]
        else:
            raise NotImplementedError("More than three partitions in the graph are not implemented.")
        
        assert np.sum(sizes) == q, "The partitioning of qubits in `sizes` does not fit the total number of qubits `q`"
        
        self.G: nx.Graph = nx.stochastic_block_model(sizes, self.p_matrix, seed=self.seed)
        self.edges: List[Tuple[int, int]] = list(self.G.edges())
        self.grouping = grouping_cascades.grouping_cascades(self.edges, self.q, self.cut_loc)
        
    def perform_grouping_operations(self) -> None:
        self.grouping.perform_grouping()
        self.grouping.limit_active_qubits_per_part()
        
    def print_grouping_results(self) -> None:
        print("Cascades prior to limiting their size")
        for item in self.grouping.edges_cascades.items():
            print(item)
        print("Cascades after limiting their size")
        for item in self.grouping.edges_cascade_lim.items():
            print(item)
        
    def count_and_print_cuts(self) -> None:
        count_blocks: int
        count_remaining_blocks: int
        count_sep: int
        count_blocks, count_remaining_blocks, count_sep = self.grouping.count_cuts()
        print("Number of blocks + remaining separate cut gates: ", count_blocks, count_remaining_blocks)
        print("Total number of separate cut gates: ", count_sep)
        print("Total number of 2-qubit gates (i.e. num of edges if 1 layer QAOA): ", len(self.grouping.edges))
        
    def draw_graph(self, save_plot: bool = False) -> None:
        pos: dict = nx.spring_layout(self.G)
        plt.figure(figsize=(3, 3))
        nx.draw(
            self.G, pos, with_labels=True, 
            node_color=(0.58, 0.0, 0.83, 0.8), edge_color='gray',
            node_size=300, 
            font_family='serif',  
            font_size=10, 
            font_color='white'
        )
        if save_plot:
            plt.savefig(f"graph_q{self.q}_cutloc{self.cut_loc}_{self.sizes}_pintra{self.p_intra}_pinter{self.p_inter}.pdf")
        plt.show()
        
    def generate_weights(self) -> None:
        self.weights = [random.uniform(0, 2 * np.pi) for _ in range(len(self.edges))]
        
    def generate_qsim_file_names(self) -> Tuple[str, str]:
        angles_str: str = "_".join([f"{a:.5f}".replace(".", "") for a in self.angles])
        sizes_str: str = "_".join(map(str, self.sizes))
        pinter_str: str = f"{self.p_inter:.2f}".replace(".", "")
        pintra_str: str = f"{self.p_intra:.2f}".replace(".", "")
        
        block_name: str = f"grouping_block_q{self.q}_cutloc{self.cut_loc}_blockgraph_seed{self.seed}_sizes_{sizes_str}_angles_{angles_str}_pinter{pinter_str}_pintra{pintra_str}_weights"
        noblock_name: str = f"grouping_noblock_q{self.q}_cutloc{self.cut_loc}_blockgraph_seed{self.seed}_sizes_{sizes_str}_angles_{angles_str}_pinter{pinter_str}_pintra{pintra_str}_weights"
        
        return block_name, noblock_name
    
    def create_qsim_files(self) -> None:
        block_name, noblock_name = self.generate_qsim_file_names()
        self.grouping.create_qaoa_qsim(self.angles, 1, noblock_name, self.weights)
        self.grouping.create_qaoa_qsim(self.angles, 2, block_name, self.weights)
    
    def run_all(self) -> None:
        self.perform_grouping_operations()
        if self.verbose:
            self.print_grouping_results()
            self.count_and_print_cuts()
            self.draw_graph(save_plot=self.save_plot)
        if len(self.weights)==0:
            self.generate_weights()
        self.create_qsim_files()

