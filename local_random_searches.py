__author__ = 'Lynnae Bryan'

import networkx as nx

from graph import GraphAL, Subgraph
from graph import count_edges, neighbor_union_subtract 

from simulated_annealing import cost_dense


# ~~~~~~~~~~~~~~~~~~ LOCAL SEARCH ~~~~~~~~~~~~~~~~~~
def naive_local_search(G, K, max_itr):
  """Starts with subgraph K, and searches neighbor states. Only explores local optima. """
  # generate starting solution
  max_cost = float("-inf")
  max_node_set = None
  max_nedges = None
  for i in range(max_itr):
    # get neighbor solution
    K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
    K_prime_nedges = count_edges(G, K_prime)
    K_cost = cost_dense(K.nvertices, K.nedges)
    # calculate costs
    K_prime_cost = cost_dense(K_prime_nvertices, K_prime_nedges)
    cur_max_set, cur_max_cost, cur_nedges = ((K.node_set, K_cost, K.nedges) if K_cost > K_prime_cost else (K_prime, K_prime_cost, K_prime_nedges)) 
    # track best solution found so far
    if cur_max_cost > max_cost:
      max_node_set = cur_max_set
      max_nedges = cur_nedges
      max_cost = cur_max_cost
    # update K 
    if K_prime_cost >= K_cost:
      K.node_set = K_prime
      K.nvertices = K_prime_nvertices
      K.nedges = K_prime_nedges
  return (max_node_set, max_cost, max_nedges, i)

# ~~~~~~~~~~~~~~~~~~ RANDOM SEARCH ~~~~~~~~~~~~~~~~~~
def random_search(G, K, max_itr):
  """Explores possible randomly-generated subgraphs of G, returns best found solution for subgraph K. Not guaranteed to be MIS. """
  best_K = K  
  best_K_cost = cost_dense(K.nvertices, K.nedges)
  for _ in range(max_itr):
    K = Subgraph(G, random_subset=True)
    K_cost = cost_dense(K.nvertices, K.nedges)
    if K_cost > best_K_cost:
      best_K = K
      best_K_cost = K_cost
  return K


if (__name__ == '__main__'):
  #FIXME for testing
  #script to run local search and such....
  nodes = 10
  edges = 23
  e = nx.dense_gnm_random_graph(nodes, edges)
  G = GraphAL(nodes, e.edges)
  K = Subgraph(G, random_subset=True)
  for _ in range(10):
    K = random_search(G, K, 100000)
    print(K.nedges, K.nvertices, K)






