__author__ = 'Lynnae Bryan'
""" A Simulated Annealing Implementation for CS-584, Spring 2019 """
"""Implementation based on Steven Skiena's psuedocode and algorithm description from
   'The Algorithm Design Manual'. See paper for full citation. """

import math
import random
import time
import datetime as dt

import networkx as nx

from graph import GraphAL, Subgraph
from graph import count_edges, density_ratio, neighbor_union_subtract 

#Global Parameters
# ~~~~~~~~~~~~~~~~~~ GLOBALS ~~~~~~~~~~~~~~~~~~
LAMBDA = 30                # Positive integer constant; penalty for violating independent set constraint 
ITR_PER_T = 800          # Number of iterations before reducing temperature
MAX_ITR = 500             # Max number of temperature changes
FREEZE = 10000           # Return current best solution after FREEZE iterations with no change to best solution 
K_CONS = 1 


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
  exp_term = -delta / t * K_CONS 
  if math.exp(exp_term) >= random.random(): 
    return True
  else: return False
  
# ~~~~~~~~~~~~~~~~~~ SIMULATED ANNEALING ~~~~~~~~~~~~~~~~~~ 
def simulated_annealing(G, K, alpha, itr_per_t, max_itr, freeze):
  """ Simulated Annealing implementation for the Maximum Independent Set Problem. """
  t = 1
  itr = 0
  max_cost = float('-inf')       
  max_node_set = None
  max_nedges = None
  max_last_update = 0
  while(itr < max_itr * itr_per_t): 
    for _ in range(itr_per_t):
      # construct neighbor via union/subtract
      K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
      # calculate cost of K and K_prime,
      K_prime_nedges = count_edges(G, K_prime)
      K_cost = cost_dense(K.nvertices, K.nedges)
      K_prime_cost = cost_dense(K_prime_nvertices, K_prime_nedges)
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
      itr += 1
      if (itr - max_last_update) > freeze:
        return (max_node_set, max_cost, max_nedges, itr) 
    # reduce temp
    if t > 0.001:
      t = alpha * t
  return (max_node_set, max_cost, max_nedges, itr)
# ~~~~~~~~~~~~~~~~~~ LOCAL SEARCH ~~~~~~~~~~~~~~~~~~
def naive_local_search(G, K, max_itr):
  """Starts with subgraph K, and searches neighbor states. Only explores local optima. """
  # generate starting solution
  max_cost = float("-inf")
  max_node_set = None
  max_nedges = None
  last_update = 0
  for itr in range(max_itr):
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
      last_update = itr
    if (last_update - itr) > FREEZE:
      return (max_node_set, max_cost, max_nedges, itr)
    # update K 
    if K_prime_cost >= K_cost:
      K.node_set = K_prime
      K.nvertices = K_prime_nvertices
      K.nedges = K_prime_nedges
  return (max_node_set, max_cost, max_nedges, itr)

# ~~~~~~~~~~~~~~~~~~ RANDOM SEARCH ~~~~~~~~~~~~~~~~~~
def random_search(G, K, max_itr):
  """Explores possible randomly-generated subgraphs of G, returns best found solution for subgraph K. Not guaranteed to be MIS. """
  best_K = K  
  best_K_cost = cost_dense(K.nvertices, K.nedges)
  last_update = 0
  for itr in range(max_itr):
    K = Subgraph(G, random_subset=True)
    K_cost = cost_dense(K.nvertices, K.nedges)
    if K_cost > best_K_cost:
      best_K = K
      best_K_cost = K_cost
      last_update = itr
  return (K, best_K_cost, itr)


# ~~~~~~~~~~~~~~~~~~ EXPERIMENT SCRIPT ~~~~~~~~~~~~~~~~~~
if (__name__ == '__main__'):
  print(dt.datetime.now())
  # alpha_list = [.8, .9, .99]
  # density_list = [.25, .75]
  # graph_size_list = [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
  # num_runs = 10         # number of runs to average
  alpha_list = [.99]
  density_list = [.2]
  graph_size_list = [1000] 
  num_runs = 5 
  for a in alpha_list:
    for dp in density_list:
      SA_runtime_fp = open(f'./Data/SA_avg_time_{dp}_alpha{a}_no_temp.txt', 'w')
      SA_data_fp = open(f'./Data/SA_data_{dp}alpha{a}_no_temp.txt', 'w')
      for G_nodes in graph_size_list:
        numerator = 0       # for calculating average
        # set number of edges by max number of edges possible times density percentage
        edges = math.floor((G_nodes * (G_nodes - 1)/2) * dp) #FIXME double check this....
        e = nx.dense_gnm_random_graph(G_nodes, edges, seed=12) #FIXME Removing Seed to see if that's problem
        G = GraphAL(G_nodes, e.edges)
        print('g size edges', G.nedges)
        for _ in range(num_runs):
          K = Subgraph(G, random_subset=True)
          print('k size edges', K.nedges)
          start = time.time() 
          max_node_set, max_cost, max_nedges, itr = simulated_annealing(G, K, a, ITR_PER_T, MAX_ITR, FREEZE)
          numerator += (time.time() - start)
          SA_data_fp.write(f'{G_nodes},{max_cost},{density_ratio(len(max_node_set), max_nedges)},{itr}\n')
          SA_data_fp.flush()
        SA_runtime_fp.write(f'{G_nodes},{numerator/num_runs}\n')
        SA_runtime_fp.flush()
        print(f'done with run {dp}, {G_nodes}, {dt.datetime.now()}')
  SA_runtime_fp.close()
  SA_data_fp.close()
