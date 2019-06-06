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
ITR_PER_T = 100      # Number of iterations before reducing temperature
MAX_ITR = 1000        # Max number of temperature changes
#FREEZE = 600        # Return current best solution after FREEZE iterations with no change to best solution 

# ~~~~~~~~~~~~~~~~~~ COST AND NEIGHBORHOOD FUNCTIONS ~~~~~~~~~~~~~~~~~~
def cost_dense(num_nodes_K, num_edges_K, t=None):
  """Objective function to be maximized by simulated annealing. Calculate cost of solution, given |K| and |E_K| where K is a subset of V and E_K is set of edges induced by K. Includes optional temperature parameter. """
  if t:
    return num_nodes_K - (LAMBDA * num_edges_K)/t #Skiena's cost function suggestion
  else:
    return num_nodes_K - (LAMBDA * num_edges_K) #Feo paper version

def accept_neighbor_solution(cost_K, cost_K_prime, t):
  """Determines whether to move to a solution with a lower cost.""" 
  delta = cost_K - cost_K_prime
  exp_term = -delta / t 
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
  for _ in range(max_itr):
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
  return (max_node_set, max_cost, max_nedges)

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
  return (K, best_K_cost)


# ~~~~~~~~~~~~~~~~~~ EXPERIMENT SCRIPT ~~~~~~~~~~~~~~~~~~

if (__name__ == '__main__'):
  density_list = [.30, .5, .75]
  graph_size_list = [15, 20]
  for dp in density_list:
    # File pointers for data collection
    SA_runtime_fp = open(f'./Data/SA/SA_avg_time_{dp}.txt', 'w')
    SA_data_fp = open(f'./Data/SA/SA_data_{dp}.txt', 'w')
    LS_runtime_fp = open(f'./Data/LS/LS_avg_time_{dp}.txt', 'w')
    LS_data_fp = open(f'./Data/LS/LS_data_{dp}.txt', 'w')
    RS_runtime_fp = open(f'./Data/RS/LS_avg_time_{dp}.txt', 'w')
    RS_data_fp = open(f'./Data/RS/RS_data_{dp}.txt', 'w')

    num_runs = 25         # number of runs to average
    for G_nodes in graph_size_list:
      numerator = 0       # for calculating average
      # set number of edges by max number of edges possible times density percentage
      edges = math.floor(((G_nodes * (G_nodes - 1))/2) * dp)
      # create random graph with NetworkX, |V| = G_nodes
      # seed all graphs with 0 for consistency
      # starting subgraph K is random for each new graph
      e = nx.dense_gnm_random_graph(G_nodes, edges, seed = 0)
      G = GraphAL(G_nodes, e.edges)
      K = Subgraph(G, random_subset=True)
      # runs for average time and data collection
      for _ in range(num_runs):
        # Run SA
        start = time.time() 
        max_node_set, max_cost, max_nedges, itr = simulated_annealing(G, K, ALPHA, ITR_PER_T, MAX_ITR, math.floor(ITR_PER_T * MAX_ITR * 0.15))
        numerator += (time.time() - start)
        SA_data_fp.write(f'{G_nodes},{max_cost},{density_ratio(len(max_node_set), max_nedges)},{itr}\n')
      SA_runtime_fp.write(f'{G_nodes},{numerator/num_runs}\n')
      SA_data_fp.flush()
      SA_runtime_fp.flush()
      # Run Local Search
      numerator = 0
      for _ in range(num_runs):
        start = time.time()  
        max_node_set, max_cost, max_nedges, i = naive_local_search(G, K, (ITR_PER_T * MAX_ITR))
        numerator += (time.time() - start)
        LS_data_fp.write(f'{G_nodes},{max_cost},{density_ratio(len(max_node_set), max_nedges)}\n')
      LS_runtime_fp.write(f'{G_nodes},{numerator/num_runs}')
      LS_data_fp.flush()
      LS_runtime_fp.flush()
      # Run Random Search
      numerator = 0
      for _ in range(num_runs):
        start = time.time()  
        best_K, best_K_cost = random_search(G, K, (ITR_PER_T * MAX_ITR))
        numerator += (time.time() - start)
        RS_data_fp.write(f'{G_nodes},{best_K_cost},{density_ratio(best_K.nvertices, best_K.nedges)}\n')
      RS_runtime_fp.write(f'{G_nodes},{numerator/num_runs}')
      RS_data_fp.flush()
      RS_runtime_fp.flush() 
        
      # print(f'avg {numerator/num_runs}')
  SA_runtime_fp.close()
  SA_data_fp.close()



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