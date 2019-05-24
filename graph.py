# import networkx as nx
import random

# modified from Prof. Bart Massey's code from CS350
# FIXME will need to cite later

class GraphAL(object):
    def __init__(self, nvertices, edges):
        "Create an adjacency-list (actually adjacency-set) graph."
        self.nvertices = nvertices
        self.neighbors = [set() for _ in range(nvertices)]
        for v1, v2 in edges:
            self.neighbors[v1].add(v2)
            self.neighbors[v2].add(v1)
    def __repr__(self):
      return f"GraphAL({self.nvertices}, {self.neighbors})"
# end Massey code adaptation

# could adapt to be Graph attribute to keep everything bundled as one Graph object. May or may not be helpful
def random_subset(Graph):
  "Result should be a random subset of digits from 0 to |V| in a list"
  nodes = [x for x in range(Graph.nvertices)]
  print(nodes)
  # may need to reimplement shuffle to analyze
  random.shuffle(nodes)
  # achtung random number; may want to change, jsut seem dumb to ahve a starting K with less than 3 vert
  N = random.randint(3, Graph.nvertices)
  print(N)
  nodes = nodes[:N]
  print(nodes)

G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
print(G)

# check for induced edge:
def edge_check(Graph, v1, v2):
  if v1 in Graph.neighbors[v2] and v2 in Graph.neighbors[v1]:
    return True
  else: return False

print(edge_check(G,4,2))
random_subset(G)
