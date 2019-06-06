__author__ = 'Lynnae Bryan'
""" A Simulated Annealing Implementation for CS-584, Spring 2019 """
"""Implementation based on Steven Skiena's psuedocode and algorithm description from
   'The Algorithm Design Manual'. See paper for full citation. """

import math
import random
import time

import networkx as nx

from graph import GraphAL, Subgraph
from graph import count_edges, density_ratio, neighbor_union_subtract 

#Global Parameters
# ~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~
# FIXME need to make sure that temp doesn't get too small
# throws exp error math error
T = 1
LAMBDA = 1          # Positive integer constant; penalty for violating independent set constraint 
ALPHA = 1         # Temperature-reducing coefficient, 0.8 <= ALPHA <= 0.99 
ITR_PER_T = 10      # Number of iterations before reducing temperature
MAX_ITR = 10        # Max number of temperature changes
#FIXME do I want to include this parameter? Need to verify the math...
K_CONS = 1  #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256 (tried range with k = (1...5), higher K values increase likelyhood of accepting worse solution. tested for T = 1 only)
FREEZE = 600        # Return current best solution after FREEZE iterations with no change to best solution 

# ~~~~~~~~~~~~~~~~~~ COST AND NEIGHBORHOOD FUNCTIONS ~~~~~~~~~~~~~~~~~~
def cost_dense(num_nodes_K, num_edges_K, t=None):
  """Objective function to be maximized by simulated annealing. Calculate cost of solution, given |K| and |E_K| where K is a subset of V and E_K is set of edges induced by K. Includes optional temperature parameter. """
  if t:
    return num_nodes_K - (LAMBDA * num_edges_K)/t #Skiena's cost function suggestion
  else:
    return num_nodes_K - (LAMBDA * num_edges_K) #Feo paper version

def accept_neighbor_solution(cost_K, cost_K_prime, T):
  """Determines whether to move to a solution with a lower cost.""" 
  delta = cost_K - cost_K_prime
  exp_term = -delta / (K_CONS * T)
  if math.exp(exp_term) >= random.random(): 
    return True
  else: return False
  
# ~~~~~~~~~~~~~~~~~~ SIMULATED ANNEALING ~~~~~~~~~~~~~~~~~~ 
def simulated_annealing(G, K, alpha, itr_per_t, max_itr, freeze):
  """ Simulated Annealing implementation for the Maximum Independent Set Problem. """
  # FIXME rework globals vs parameters.
  t = 1
  itr = 0
  max_cost = float('-inf')       
  max_node_set = None
  max_nedges = None
  max_last_update = 0
  while(itr < max_itr * itr_per_t): #FIXME
    for _ in range(itr_per_t):
      # construct neighbor via union/subtract
      K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
      # calculate cost of K and K_prime,
      K_prime_nedges = count_edges(G, K_prime)
      K_cost = cost_dense(K.nvertices, K.nedges, t)
      K_prime_cost = cost_dense(K_prime_nvertices, K_prime_nedges, t)
      # See if K or K' cost is better than current maximum
      cur_max_set, cur_max_cost, cur_nedges = ((K.node_set, K_cost, K.nedges) if K_cost > K_prime_cost else (K_prime, K_prime_cost, K_prime_nedges))
      if cur_max_cost > max_cost:
        max_node_set = cur_max_set
        max_nedges = cur_nedges
        max_cost = cur_max_cost
        max_last_update = itr
      # relies on MAXIMIZING cost function. 
      if K_prime_cost >= K_cost:
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
        K.nedges = K_prime_nedges
      # accept worse solution by probability, a_n_s returns t/f
      elif accept_neighbor_solution(K_cost, K_prime_cost, t):
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
        K.nedges = K_prime_nedges
      if (itr - max_last_update) > freeze:
        return (max_node_set, max_cost, max_nedges, itr) 
      itr += 1
    # reduce temp
    #FIXME approx
    if t > 0.01:
      t = ALPHA * t
  return (max_node_set, max_cost, max_nedges, itr)
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


# ~~~~~~~~~~~~~~~~~~ EXPERIMENT SCRIPT ~~~~~~~~~~~~~~~~~~

if (__name__ == '__main__'):
  # say 4 graphs of different sizes (test on 10, 15, 20, 25)
  # at 80% density
  density_list = [.30, .5, .75]
  graph_size_list = [15, 20]
  for dp in density_list:
    runtime_fp = open(f'./Data/avg_time_{dp}.txt', 'w')
    data_fp = open(f'./Data/SA_run_data_{dp}.txt', 'w')
    num_runs = 500
    for G_nodes in graph_size_list:
      numerator = 0
      edges = math.floor(((G_nodes * (G_nodes - 1))/2) * dp)
      e = nx.dense_gnm_random_graph(G_nodes, edges, seed = 0)
      G = GraphAL(G_nodes, e.edges)
      K = Subgraph(G, random_subset=True)
      
      for j in range(num_runs):
        start = time.time() 
        max_node_set, max_cost, max_nedges, itr = simulated_annealing(G, K, ALPHA, ITR_PER_T, MAX_ITR, FREEZE)
        total_sec = time.time() - start
        numerator += total_sec
        data_fp.write(f'{G_nodes},{max_cost},{density_ratio(len(max_node_set), max_nedges)},{itr}\n')
      # print(f'avg {numerator/num_runs}')
      runtime_fp.write(f'{G_nodes},{numerator/num_runs}\n')
      data_fp.flush()
      runtime_fp.flush()
  runtime_fp.close()
  data_fp.close()



  # G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
  # nodes = 5
  # edges = 9 
  # G = GraphAL(nodes, e.edges)
  # K = Subgraph(G, random_subset=True)
  # maxns, maxcost, maxned, i = simulated_annealing(G, K, 1, 10, MAX_ITR, FREEZE)
  # print(i)


  # for i in range(11):
    # best, cost, nedges, itr_stop = simulated_annealing(G,K)
    # print(f'cost {cost}, itr stop: {itr_stop}, cardinality {len(best)}')
    # print(f'Density of final solution: {density_ratio(len(best), nedges)}')
    
  # for _ in range(10):
  #   best, cost, itr_stop = simulated_annealing(G)
  #   print(f'best {best} cost {cost}, itr stop: {itr_stop}')