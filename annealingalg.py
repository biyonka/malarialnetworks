from __future__ import print_function
import math
import random
import simanneal
from simanneal import Annealer

import argparse
from math import log
import networkx as nx
import sys
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances

# def distance(a, b):
#     """Calculates distance between two latitude-longitude coordinates."""
#     R = 3963  # radius of Earth (miles)
#     lat1, lon1 = math.radians(a[0]), math.radians(a[1])
#     lat2, lon2 = math.radians(b[0]), math.radians(b[1])
#     return math.acos(math.sin(lat1) * math.sin(lat2) +
#                      math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2)) * R

TAU = 0.15
PAGE_RANK = 'page_rank'
MODULE_ID = 'module_id'

def log2(prob):
    "Returns the log of prob in base 2"
    return log(prob, 2)

def entropy1(prob):
    """Half of the entropy function, as used in the InfoMap paper.
    entropy1(p) = p * log2(p)
    """
    if prob == 0:
        return 0
    return prob * log2(prob)

def load_and_process_graph(filename):
    """Load the graph, normalize edge weights, compute pagerank, and store all
    this back in node data."""
    # Load the graph
    graph = nx.DiGraph(nx.read_pajek(filename))
    #graph = nx.read_pajek(filename)
    print ("Loaded a graph (%d nodes, %d edges)" % (len(graph),
            len(graph.edges())))
    # Compute the normalized edge weights
    for node in graph:
        edges = graph.edges(node, data=True)
        total_weight = sum([data['weight'] for (_, _, data) in edges])
        for (_, _, data) in edges:
            data['weight'] = data['weight'] / total_weight
    # Get its PageRank, alpha is 1-tau where [RAB2009 says \tau=0.15]
    page_ranks = nx.pagerank(graph, alpha=1-TAU)
    for (node, page_rank) in page_ranks.items():
        graph.node[node][PAGE_RANK] = page_rank
    return graph

def load_coordinates(filename):
    field_names = ['X', 'Y', "w"]
    coords = pd.read_csv(filename, header=None, names=field_names)
    coords = coords.loc[:,["X","Y"]]
    #coords = coords.as_matrix()
    return coords

class Module:
    """Stores the information about a single module"""
    def __init__(self, module_id, nodes, graph):
        self.module_id = module_id
        self.nodes = frozenset(nodes)
        self.graph = graph
        self.prop_nodes = 1 - float(len(self.nodes)) / len(graph)
        # Set the module_id for every node
        for node in nodes:
            graph.node[node][MODULE_ID] = module_id
        # Compute the total PageRank
        self.total_pr = sum([graph.node[node][PAGE_RANK] for node in nodes])
        # Compute q_out, the exit probability of this module
        # .. Left half: tau * (n - n_i) / n * sum{alpha in i}(p_alpha)
        self.q_out = self.total_pr * TAU * self.prop_nodes
        # .. Right half: (1-tau) * sum{alpha in i}(sum{beta not in i}
        #                  p_alpha weight_alpha,beta)
        # This is what's in [RAB2009 eq. 6]. But it's apparently wrong if
        # node alpha has no out-edges, which is not in the paper.
        # ..
        # Implementing it with Seung-Hee's correction about dangling nodes
        for node in self.nodes:
            edges = graph.edges(node, data=True)
            page_rank = graph.node[node][PAGE_RANK]
            if len(edges) == 0:
                self.q_out += page_rank * self.prop_nodes * (1 - TAU)
                continue
            for (_, dest, data) in edges:
                if dest not in self.nodes:
                    self.q_out += page_rank * data['weight'] * (1 - TAU)
        self.q_plus_p = self.q_out + self.total_pr

    def get_codebook_length(self):
        "Computes module codebook length according to [RAB2009, eq. 3]"
        first = -entropy1(self.q_out / self.q_plus_p)
        second = -sum( \
                [entropy1(self.graph.node[node][PAGE_RANK]/self.q_plus_p) \
                    for node in self.nodes])
        return (self.q_plus_p) * (first + second)

class TravellingSalesmanProblem(Annealer):

    """Test annealer with a travelling salesman problem.
    """

    # pass extra data (the distance matrix) into the constructor
    def __init__(self, state, module, graph, coordinates):
       # self.distance_matrix = distance_matrix
        self.graph = graph
        self.total_pr_entropy = sum([entropy1(graph.node[node][PAGE_RANK]) \
                for node in graph])
        self.module = [Module(module_id, mod, graph) \
                for (module_id, mod) in enumerate(state)]

        d = 0
        for mod in module:
            for elem in range(len(mod)):
                mod[elem] = int(mod[elem])    
        for mod in module:
            m = coordinates.loc[mod,]
            d += np.mean(pairwise_distances(m, metric='euclidean'))
        self.d = d 

        super(TravellingSalesmanProblem, self).__init__(state)  # important!
    
    def move(self):
        """Swaps first element of two randomly chosen clusters."""
        a = random.randint(0, len(self.state) - 1)
        b = random.randint(0, len(self.state) - 1)
        self.state[a][0], self.state[b][0] = self.state[b][0], self.state[a][0]

    def energy(self):
        "Compute the MDL of this clustering according to [RAB2009, eq. 4]"
        total_qout = 0
        total_qout_entropy = 0
        total_both_entropy = 0
        for mod in self.module:
            q_out = mod.q_out
            total_qout += q_out
            total_qout_entropy += entropy1(q_out)
            total_both_entropy += entropy1(mod.q_plus_p)
        term1 = entropy1(total_qout)
        term2 = -2 * total_qout_entropy
        term3 = -self.total_pr_entropy
        term4 = total_both_entropy
        term5 = self.d
        
        return term1 + term2 + term3 + term4 + term5

if __name__ == '__main__':
    #networkX Digraph
    graph = load_and_process_graph("houses.net")#(options.graph_filename)

    #coords is pandas dataframe of the coordinates
    coords = load_coordinates("coordinates.csv")
    coords.index = np.arange(0, len(coords))

    # single_nodes is the "trivial" module mapping
    # initial module, as a list of lists
    single_nodes = [[nodes] for nodes in graph]
    single_nodes[0] = ['0','1']
    single_nodes.remove(single_nodes[1]) 
    init_state = single_nodes

    tsp = TravellingSalesmanProblem(state = init_state, module = init_state, graph = graph, coordinates = coords)
    tsp.steps = 100000
    # since our state is just a list, slice is the fastest way to copy
    tsp.copy_strategy = "slice"
    state, e = tsp.anneal()
        
    print()
    print("%i mile route:" % e)
    for city in state:
        print("\t", city) 