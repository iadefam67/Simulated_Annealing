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
      total_itr += 1
      if total_itr > max_itr:
        return (best_K, best_K_cost)
    if t > 0.0001:
      t = t * alpha   #reduce temp

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
    # calculate costs
    K_cost = cost_dense(K.nvertices, K.nedges)
    K_prime_cost = cost_dense(K_prime_nvertices, K_prime_nedges)
    cur_max_set, cur_max_cost, cur_nedges = ((K.node_set, K_cost, K.nedges) if K_cost > K_prime_cost else (K_prime, K_prime_cost, K_prime_nedges)) 
    # track best solution found so far
    if cur_max_cost > max_cost:
      max_node_set = cur_max_set
      max_nedges = cur_nedges
      max_cost = cur_max_cost
      last_update = itr
    if (itr - last_update) > FREEZE:
      return (max_node_set, max_cost, max_nedges)
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
  for itr in range(max_itr):
    K = Subgraph(G, random_subset=True)
    K_cost = cost_dense(K.nvertices, K.nedges)
    if K_cost > best_K_cost:
      best_K = K
      best_K_cost = K_cost
  return (best_K.nvertices, best_K_cost)

# ~~~~~~~~~~~~~~~~~~ EXPERIMENT SCRIPTS ~~~~~~~~~~~~~~~~~~
if (__name__ == '__main__'):
  # LOCAL SEARCH
  print(dt.datetime.now())
  header = 'Num_G_nodes, Avg_time, Max_cost, Avg_cost, Percent_MIS,Density'
  num_runs = 5 
  density_list = [.1, .7]
  fp = open(f'./Data/RS_data.txt', 'w')
  fp.write(f'{header}\n')
  for dp in density_list:
    for G_nodes in range(15, 21):
      time_avg = 0       # for calculating average
      edges = math.floor((G_nodes * (G_nodes - 1)/2) * dp) 
      e = nx.dense_gnm_random_graph(G_nodes, edges) 
      G = GraphAL(G_nodes, e.edges)
      max_K_set = None
      max_sol_cost = float("-inf")
      avg_cost = 0
      MIS_avg = 0
      for _ in range(num_runs):
        K = Subgraph(G, random_subset=True)
        start = time.time() 
        K_sol_len, K_sol_cost = random_search(G,K, MAX_ITR)
        avg_cost += K_sol_cost
        if K_sol_cost > max_sol_cost:
          max_sol_cost = K_sol_cost
        time_avg += (time.time() - start)
        if K_sol_cost == K_sol_len:
          MIS_avg += 1
      fp.write(f'{G_nodes},{time_avg/num_runs},{max_sol_cost},{avg_cost/num_runs},{MIS_avg/num_runs},{dp}\n')
      fp.flush()
      print(f'done with run {dp}, {G_nodes}, {dt.datetime.now()}')
  fp.close()



  # SIMULATED ANNEALING
  print(dt.datetime.now())
  alpha_list = [.99, .9, .8]
  density_list = [.1, .7]
  header = 'Num_G_nodes, Avg_time, Max_cost, Avg_cost, Netx_avg_cost, Percent_MIS, Density'
  num_runs = 10 
  X_avg = 0
  for a in alpha_list:
    fp = open(f'./Data/SA_data_a_{a}.txt', 'w')
    fp.write(f'{header}\n')
    for dp in density_list:
      for G_nodes in range(10,20):
        numerator = 0       # for calculating average
        edges = math.floor((G_nodes * (G_nodes - 1)/2) * dp) 
        e = nx.dense_gnm_random_graph(G_nodes, edges) 
        G = GraphAL(G_nodes, e.edges)
        X_avg = 0
        max_K_set = None
        max_sol_cost = float("-inf")
        avg_cost = 0
        MIS_avg = 0
        for _ in range(num_runs):
          X = nx.maximal_independent_set(e)
          X_avg += len(X)
          K = Subgraph(G, random_subset=True)
          start = time.time() 
          K_sol, K_sol_cost = simulated_annealing(G, a, ITR_PER_T, MAX_ITR, FREEZE)
          avg_cost += K_sol_cost
          if K_sol_cost > max_sol_cost:
            max_sol_cost = K_sol_cost
          numerator += (time.time() - start)
          if K_sol_cost == len(K_sol):
            MIS_avg += 1
        fp.write(f'{G_nodes},{numerator/num_runs},{max_sol_cost},{avg_cost/num_runs},{X_avg/num_runs},{MIS_avg/num_runs},{dp}\n')
        fp.flush()
        print(f'done with run {dp}, {G_nodes}, {dt.datetime.now()}')
  fp.close()


