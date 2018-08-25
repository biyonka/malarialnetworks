from __future__ import print_function
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.offline as offline
from plotly.graph_objs import Scatter, Layout
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
offline.init_notebook_mode(connected=True)
from discreteMarkovChain import markovChain
import networkx as nx
import csv
import math
import random
import simanneal
from simanneal import Annealer
import argparse
from math import log
import networkx as nx
import sys
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import hamming

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
    return coords

class GeoInfomap(Annealer):

    """Test annealer with a travelling salesman problem.
    """

    # pass extra data (the distance matrix) into the constructor
    def __init__(self, state, module, graph, coordinates):
        self.graph = graph
        self.total_pr_entropy = sum([entropy1(graph.node[node][PAGE_RANK]) \
                for node in graph])
#         self.module = [Module(module_id, mod, graph) \
#                 for (module_id, mod) in enumerate(state)]
        d = 0
        for mod in module:
            for elem in range(len(mod)):
                mod[elem] = int(mod[elem])    
        for mod in module:
            m = coordinates.loc[mod,]
            d += np.mean(pairwise_distances(m, metric='euclidean'))
        self.d = d 

        super(GeoInfomap, self).__init__(state)  # important!
    def move(self):
        #converts list of node lists into a 1D array of community labels
        cluster_labels = []
        label = 0
        for cluster in self.state:
            for elem in cluster:
                cluster_labels.append(label)
            label += 1

    #flatten self.state
        flat_list = [item for sublist in self.state for item in sublist]
        #sort cluster labels list by node label (not community label)
        cluster_labels = [cluster_labels[flat_list.index(i)] for i in flat_list]
        a = random.randint(0, len(cluster_labels)-1)
        change_node = cluster_labels[a] 
        #if current label is 0, change to 1 to increase hamming distance by 1
        if change_node == 0:
            cluster_labels[a] = 1
        else:
            updown = random.randint(0, 1)
            if updown == 0:
                cluster_labels[a] -= 1
            cluster_labels[a] += 1

        #convert back to list of lists
        new_state = []
        labels = set(cluster_labels)
        for j in labels:
            cluster = []
            indices = [i for i, x in enumerate(cluster_labels) if x == j]
            cluster.extend(list(np.array(flat_list)[indices]))
            new_state.append(cluster)

        self.state = new_state

    def energy(self):
        "Compute the MDL of this clustering according to [RAB2009, eq. 4]"
        graph = self.graph
        
        total_qout = 0
        total_qout_entropy = 0
        total_both_entropy = 0
        for mod in self.state:
            nodes = frozenset(mod)
            prop_nodes = 1 - float(len(nodes)) / len(graph)
            total_pr = sum([graph.node[str(node)][PAGE_RANK] for node in nodes])
            q_out = total_pr * TAU * prop_nodes

            for node in mod:
                edges = graph.edges(str(node), data=True)
                page_rank = graph.node[str(node)][PAGE_RANK]
                if len(edges) == 0:
                    q_out += page_rank * prop_nodes * (1 - TAU)
                    continue
                for (_, dest, data) in edges:
                    if dest not in self.state:
                        q_out += page_rank * data['weight'] * (1 - TAU)
                q_plus_p = q_out + total_pr
        
            q_out = q_out
            total_qout += q_out
            total_qout_entropy += entropy1(q_out)
            total_both_entropy += entropy1(q_plus_p)
        term1 = entropy1(total_qout)
        term2 = -2 * total_qout_entropy
        term3 = -self.total_pr_entropy
        term4 = total_both_entropy
        term5 = self.d
        total =  term1 + term2 + term3 + term4 + 0.1*term5
        #print(total)
        return total

def main(argv):
    "Read the supplied graph and modules and output MDL"
    # Read the arguments
    parser = argparse.ArgumentParser(description="Calculate the infomap")
    parser.add_argument('-g', '--graph-filename', type=argparse.FileType('r'),
            help="the .net file to use as the graph", required=True)
    parser.add_argument('-c', '--coordinates', type=argparse.FileType('r'),
        help="the .coordinate file to use as the graph", required=True)
    parser.add_argument('-s', '--steps', default="10000",
        help="number of steps to use for the simulated annealing optimization")
    parser.add_argument('-m', '--module-filename', default="2009_figure3a.mod",
            help="the .mod file to use as the clustering")
    options = parser.parse_args(argv[1:])

    graph = load_and_process_graph(options.graph_filename.name)
    coords = load_coordinates(options.coordinates.name)
    steps = int(options.steps)
    coords.index = np.arange(0, len(coords))
    print(options.graph_filename.name)
    # single_nodes is the "trivial" module mapping
    single_nodes = [[nodes] for nodes in graph]

    #If clustering provided, use it.
    try:
        modules = [line.strip().split() for line in options.module_filename]
    except IOError:
        print (">>", sys.exc_info()[0])
        print (">> No .mod file provided, or error reading it")
        print (">> Using default clustering of every node in its own module")
        modules = single_nodes
    init_state = single_nodes

    # clustering = Clustering(graph, modules)
    # print ("This clustering has MDL %.2f (Index %.2f, Module %.2f)" % \
    #     (clustering.get_mdl(), clustering.get_index_codelength(),
    #             clustering.get_module_codelength()))

    gi = GeoInfomap(state = init_state, module = init_state, graph = graph, coordinates = coords)
    gi.steps = steps
    # since our state is just a list, slice is the fastest way to copy
    gi.copy_strategy = "slice"
    state, e = gi.anneal()
     
    print()
    for city in state:
        print("\t", city)   

    #create data frame assigning each coordinate to a cluster
    cluster = np.array(1)
    i = 0
    for j in state:
        cluster = np.append(cluster, np.array([i for elem in j]))
        i += 1
    cluster = cluster.flatten()[1:]
    coords['cluster'] = cluster
    output = options.graph_filename.name.split('/')[1]
    coords.to_csv('ArtificialLandscapes/ClusterOutput/' + output[:-4] + '_clustering.csv')

if __name__ == "__main__":
    main(sys.argv)
