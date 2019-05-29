# import networkx as nx
import random

# modified from Massey, CS-350

#ok rep K as a list, K_prime is a list, and use edge count, since all we need is len K, K_prime and an edge count. Cost also doesn't need to be part of the object. 


class GraphAL(object):
  def __init__(self, nvertices, edges):
    "Create an adjacency-list style graph, using sets."
    self.nvertices = nvertices
    self.neighbors = [set() for _ in range(nvertices)]
    #SA attributes
    # build adjacency sets 
    for v1, v2 in edges:
      self.neighbors[v1].add(v2)
      self.neighbors[v2].add(v1)
  def __repr__(self):
    "Print graph representation."
    return f"GraphAL:\n Vertices: {self.nvertices}, Neighbors:{self.neighbors}"

class Subgraph(object):
  def __init__(self, Graph):
    self.node_set = None
    self.nvertices= None
    self.random_subset(Graph)
  def random_subset(self, Graph):
    "Result should be a random subset of digits from 0 to |V| in a list"
    nodes = [x for x in range(Graph.nvertices)]
    # FIXME placeholder shuffle
    random.shuffle(nodes)
    N = random.randint(3, Graph.nvertices)
    self.node_set = set(nodes[:N])
    self.nvertices= len(self.node_set)
  def __repr__(self):
    return f"Subgraph node list: {self.node_set}, Cardinality: {self.nvertices}"
  
G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
K = Subgraph(G)
print(G)
print(K)

def count_edges(GraphAL, Subgraph):
  count = 0
  for i in Subgraph.node_set:
    for v in GraphAL.neighbors[i]:
      if v in Subgraph.node_set:
        count += 1 
  assert count % 2 == 0
  return count/2

print(count_edges(G, K))

# check for induced edge:
def edge_check(Graph, v1, v2):
  if v1 in Graph.neighbors[v2] and v2 in Graph.neighbors[v1]:
    return True
  else: return False

