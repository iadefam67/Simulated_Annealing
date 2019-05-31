# import networkx as nx
import random

# modified from Massey, CS-350
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
  """Stores a list of subset of nodes from a GraphAL object, can be randomly generated when created"""
  def __init__(self, Graph, random_subset=False):
    self.node_set = None
    self.nvertices = None
  
    if random_subset == True:
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
    return f"\nSubgraph node list: {self.node_set}, Cardinality: {self.nvertices}"

def count_edges(GraphAL, node_set):
  if node_set == None:
    print("No subgraph created. Run 'random_subset()'. ")
    return False
  count = 0
  for i in node_set:
    for v in GraphAL.neighbors[i]:
      if v in node_set:
        count += 1 
  assert count % 2 == 0
  return count/2

# check for induced edge:
def edge_check(Graph, v1, v2):
  if v1 in Graph.neighbors[v2] and v2 in Graph.neighbors[v1]:
    return True
  else: return False

G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
K = Subgraph(G, random_subset=True)
# print(G)
# print(K)
# print("K edges:", count_edges(G, K.node_set))


c = count_edges(G, {0,1,3})