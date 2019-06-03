import networkx as nx
from graph import * 
import math

#Global Parameters
T = 1
LAMBDA = 20          #should be a pos integer constant (Leo, p 870)
ALPHA = .99        # temperature reducing coefficient, try 0.8 <= ALPHA <= 0.99 (higher alphas reduce temp more slowly)
TRIALS_PER_T = 10
ITR_PER_T = 5000    # num of iterations before temp change
MAX_ITR = 500     # relate to itr per t?
#FIXME do I want to include this parameter? Need to verify the math...
K_CONS = 1.455 #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256 (tried range with k = (1...5), higher K values increase likelyhood of accepting worse solution. tested for T = 1 only)
FREEZE = 600     # how many cycles with no change in solution before returning

# cost function notation from Feo paper (greedy alg vs SA)
# makes more sense to MAXIMIZE, since solution is num nodes in independent set 
def cost_dense(num_nodes_K, num_edges_K, t=None):
  """Calculate cost of solution, given |K| and |E_K| where K is a (proper? why proper) subset of V and E_K is set of edges induced by K"""
  if t:
    return num_nodes_K - (LAMBDA * num_edges_K)/t #Skiena's cost function suggestion
  else:
    return num_nodes_K - (LAMBDA * num_edges_K) #Feo paper version

# skiena pg 258
def accept_neighbor_solution(cost_K, cost_K_prime, T):
  #FIXME should verify and think through the math here....
  """Probability of accepting a solution with a worse cost""" 
  delta = cost_K - cost_K_prime
  exp_term = -delta / (K_CONS * T)
  #FIXME optimize exp function
  if math.exp(exp_term) >= random.random(): 
    return True
  else: return False

# this is the version that accepts non-MIS neighbor solutions
def neighbor_union_subtract(G, K):
  # get index of random node in V
  v = random.randint(0, G.nvertices - 1)
  if v in K.node_set:
    return ((K.node_set).difference({v}), K.nvertices - 1)
  else:
    return ((K.node_set).union({v}), K.nvertices + 1)

# FIXME need to make sure that temp doesn't get too small
# throws exp error math error
def simulated_annealing(G):
  t = T
  freeze = FREEZE
  itr = 0
  max_cost = float('-inf')        #want to minimize cost
  max_node_set = None
  max_last_update = 0
  K = Subgraph(G, random_subset=True)
  count = 0
  while(itr < MAX_ITR): #FIXME
    for _ in range(ITR_PER_T):
      # construct neighbor via union/subtract
      K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
      # calculate cost of K and K_prime,
      K_prime_nedges = count_edges(G, K_prime)
      #FIXME should use a flag here, to check if already calculated last solution for K_cost
      K_cost = cost_dense(K.nvertices, K.nedges, t)
      K_prime_cost = cost_dense(K_prime_nvertices, K_prime_nedges, t)
      # See if K or K' cost is better than current maximum
      cur_max_set, cur_max_cost, cur_nedges = ((K.node_set, K_cost, K.nedges) if K_cost > K_prime_cost else (K_prime, K_prime_cost, K_prime_nedges))
      if cur_max_cost > max_cost:
        max_node_set = cur_max_set
        max_nedges = cur_nedges
        max_cost = cur_max_cost
        max_last_update = count
      # relies on MAXIMIZING cost function. 
      # FIXME I think this is wrong... it should accept worse solution either way, not just if the cost isnt better?
      if K_prime_cost >= K_cost:
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
        K.nedges = K_prime_nedges
      # accept worse solution by probability, a_n_s returns t/f
      elif accept_neighbor_solution(K_cost, K_prime_cost, t):
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
        K.nedges = K_prime_nedges
      if (count - max_last_update) > freeze:
        #FIXME 
        if density_ratio(len(max_node_set), max_nedges) == 0: 
          print(f'final temp: {t}')
          return (max_node_set, max_cost, max_nedges, count) 
        else:
          freeze = 0
      count += 1
    # reduce temp
    t = ALPHA * t
    itr += 1
  print(f"Best solution found in {itr*MAX_ITR} iterations. May not be maximal IS.")
  return (max_node_set, max_cost, max_nedges, count)

# G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
nodes = 50 
edges = 200
e = nx.dense_gnm_random_graph(nodes, edges)
G = GraphAL(nodes, e.edges)
for i in range(11):
  best, cost, nedges, itr_stop = simulated_annealing(G)
  print(f'cost {cost}, itr stop: {itr_stop}, cardinality {len(best)}')
  print(f'Density of final solution: {density_ratio(len(best), nedges)}')
  
# for _ in range(10):
#   best, cost, itr_stop = simulated_annealing(G)
#   print(f'best {best} cost {cost}, itr stop: {itr_stop}')