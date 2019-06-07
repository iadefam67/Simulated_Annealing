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
LAMBDA = 4                # Positive integer constant; penalty for violating independent set constraint 
ITR_PER_T = 1000          # Number of iterations before reducing temperature
FREEZE = 1500          # Return current best solution after FREEZE iterations with no change to best solution 
MAX_ITR = 1000000             # Max number of temperature changes
K_CONS = .25 


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
def simulated_annealing(G, alpha, itr_per_t, max_itr, freeze):
  """ Simulated Annealing implementation for the Maximum Independent Set Problem. """
  K = Subgraph(G, random_subset=True)
  t = 1
  best_K = None
  best_K_cost = float("-inf")
  no_change = 0
  total_itr = 0
  while(True):
    for i in range(max_itr):
      K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
      K_prime_nedges = count_edges(G, K_prime)
      cost_K = cost_dense(K.nvertices, K.nedges)
      cost_K_prime = cost_dense(K_prime_nvertices, K_prime_nedges)
      # print(cost_K, cost_K_prime)
      if cost_K_prime >= cost_K:
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices  
        K.nedges = K_prime_nedges
        no_change = 0
        if cost_K_prime > best_K_cost:
          best_K = K_prime
          best_K_cost = cost_K_prime
          total_itr += 1
          continue
      elif accept_neighbor_solution(cost_K, cost_K_prime, t):
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices  
        K.nedges = K_prime_nedges 
      else:
        no_change += 1
      if no_change >= freeze:
        return (best_K, best_K_cost)
      # print(cost_K)
      total_itr += 1
      if total_itr > max_itr:
        return (best_K, best_K_cost)
    if t > 0.0001:
      t = t * alpha   #reduce temp
    # print(t)

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
  for itr in range(max_itr):
    K = Subgraph(G, random_subset=True)
    K_cost = cost_dense(K.nvertices, K.nedges)
    if K_cost > best_K_cost:
      best_K = K
      best_K_cost = K_cost
  return (K, best_K_cost, itr)

# ~~~~~~~~~~~~~~~~~~ EXPERIMENT SCRIPT ~~~~~~~~~~~~~~~~~~
if (__name__ == '__main__'):
  print(dt.datetime.now())
  # alpha_list = [.8, .9, .99]
  # density_list = [.25, .5, .75]
  # graph_size_list = [300, 400, 500, 600, 700, 800, 90]
  alpha_list = [.99]
  density_list = [.1]
  graph_size_list = [200]
  num_runs = 2 
  X_avg = 0
  for a in alpha_list:
    for dp in density_list:
      # SA_runtime_fp = open(f'./Data/SA_avg_time_{dp}_alpha{a}_with_temp.txt', 'w')
      # SA_data_fp = open(f'./Data/SA_data_{dp}alpha{a}_with_temp.txt', 'w')
      for G_nodes in range(10,200):
        numerator = 0       # for calculating average
        # set number of edges by max number of edges possible times density percentage
        edges = math.floor((G_nodes * (G_nodes - 1)/2) * dp) #FIXME double check this....
        e = nx.dense_gnm_random_graph(G_nodes, edges) #FIXME Removing Seed to see if that's problem
        G = GraphAL(G_nodes, e.edges)
        # print('g size edges', G.nedges)
        X_avg = 0
        for _ in range(num_runs):
          X = nx.maximal_independent_set(e)
          X_avg += len(X)
          K = Subgraph(G, random_subset=True)
          # print('k size edges', K.nedges)
          start = time.time() 
          K_sol, K_sol_cost = simulated_annealing(G, a, ITR_PER_T, MAX_ITR, FREEZE)
          print(f'SA nodes {len(K_sol)}, cost {K_sol_cost}')
          numerator += (time.time() - start)
          # SA_data_fp.write(f'{G_nodes},{max_cost},{density_ratio(len(max_node_set), max_nedges)},{itr}\n')
          # SA_data_fp.flush()
        # SA_runtime_fp.write(f'{G_nodes},{numerator/num_runs}\n')
        # SA_runtime_fp.flush()
        print(f'done with run {dp}, {G_nodes}, {dt.datetime.now()}')
        print(f'X avg len {X_avg/num_runs}')
  # SA_runtime_fp.close()
  # SA_data_fp.close()
