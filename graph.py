__author__ = 'Lynnae Bryan'

import random


# ~~~~~~~~~~~~~~~~~~ GRAPH CLASSES ~~~~~~~~~~~~~~~~~~
class GraphAL(object):
  """An adjacency-set based graph representation, based on code written for Prof. Barton Massey's CS-350 Algorithms course. """
  def __init__(self, nvertices, edges):
    self.nvertices = nvertices
    self.nedges = len(edges)
    self.neighbors = [set() for _ in range(nvertices)]
    #SA attributes
    # build adjacency sets 
    for v1, v2 in edges:
      self.neighbors[v1].add(v2)
      self.neighbors[v2].add(v1)
  def __repr__(self):
    return f"GraphAL:\n Vertices: {self.nvertices}, Neighbors:{self.neighbors}"

class Subgraph(object):
  """ Subgraph class. Subset of nodes from parent graph. Parent graph is GraphAL."""
  def __init__(self, Graph, random_subset=False):
    self.node_set = None
    self.nvertices = None
    #FIXME check if this is still being used
    self.parent_graph = Graph
    self.nedges = None
    if random_subset == True:
      self.random_subset(Graph)
  def random_subset(self, Graph):
    """Returns a shuffled list of vertices, with a random length."""
    nodes = [x for x in range(Graph.nvertices)]
    # FIXME placeholder shuffle
    random.shuffle(nodes)
    N = random.randint(3, Graph.nvertices)
    self.node_set = set(nodes[:N])
    self.nvertices= len(self.node_set)
    self.nedges = count_edges(self.parent_graph, self.node_set)
  def __repr__(self):
    return (f"\nSubgraph node list: {self.node_set}\nCardinality: {self.nvertices}, Density: {density_ratio(self.nvertices, self.nedges)}")


# ~~~~~~~~~~~~~~~~~~ GRAPH OPERATIONS ~~~~~~~~~~~~~~~~~~
def count_edges(GraphAL, node_set):
  """Count the number of edges induced by a set of nodes (node_set) for a parent graph (GraphAL)."""
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

def density_ratio(nvertices, nedges):
  """ Calculate ratio of edges to total number of possible edges. If density ratio is 0, the node set is independent."""
  try:
    assert nedges <= (nvertices * nvertices - 1)/2
  except AssertionError:
    print('Impossible number of edges.')
  return (2 * nedges)/(nvertices * (nvertices - 1))

def neighbor_union_subtract(G, K):
  """Select a random node v from V. If node already exists in subgraph K, subtract it from K nodes, otherwise union with K node set. """
  v = random.randint(0, G.nvertices - 1)
  if v in K.node_set:
    return ((K.node_set).difference({v}), K.nvertices - 1)
  else:
    return ((K.node_set).union({v}), K.nvertices + 1)

# FIXME remove?
# ~~~~~~~~~~~~~~~~~~ TEST SCRIPT ~~~~~~~~~~~~~~~~~~
if (__name__ == '__main__'):
  density = density_ratio(10, 45)
  print(density)
  G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
  K = Subgraph(G, random_subset=True)
  # print(G)
  # print(K)
  # print("K edges:", K.nedges) 

