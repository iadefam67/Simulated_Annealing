# import networkx as nx
import random

# modified from Massey, CS-350
class GraphAL(object):
  def __init__(self, nvertices, edges):
    "Create an adjacency-list sactually adjacency-set) graph."
    self.nvertices = nvertices
    self.neighbors = [set() for _ in range(nvertices)]
    self.subgraph_K = None
    self.nvertices_K = None
    # how to keep track of edges in subgraph...
    self.cost = float('inf')
    for v1, v2 in edges:
      self.neighbors[v1].add(v2)
      self.neighbors[v2].add(v1)
    self.random_subset()
  def random_subset(self):
    "Result should be a random subset of digits from 0 to |V| in a list"
    nodes = [x for x in range(self.nvertices)]
    # FIXME placeholder shuffle
    random.shuffle(nodes)
    N = random.randint(3, self.nvertices)
    self.subgraph_K = set(nodes[:N])
    self.nvertices_K = len(self.subgraph_K)
  def __repr__(self):
    "Print graph representation."
    return f"GraphAL:\n Vert: {self.nvertices}, Neighbors:{self.neighbors}, Sub: {self.subgraph_K}, Sub num:{self.nvertices_K}"

G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
print(G)


# check for induced edge:
def edge_check(Graph, v1, v2):
  if v1 in Graph.neighbors[v2] and v2 in Graph.neighbors[v1]:
    return True
  else: return False

