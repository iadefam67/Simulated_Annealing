from graph import * 
import math

#Global Parameters
LAMBDA = 1          #should be a pos integer constant (Leo, p 870)
ALPHA = 0.8         # temperature reducing coefficient
TRIALS_PER_T = 10
ITR_PER_T = 10      # num of iterations before temp change
MAX_ITR = 100       # relate to itr per t?

K = 1 #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256
T = 1

# cost function notation from Feo paper (greedy alg vs SA)
# MIN cost 
def cost(num_nodes_K, num_edges_K):
  """Calculate cost of solution, given |K| and |E_K| where K is a (proper? why proper) subset of V and E_K is set of edges induced by K"""
  return -num_nodes_K + (LAMBDA * num_edges_K) #Feo paper version
  #return num_K - (LAMBDA * E_K)/T #Skiena's cost function suggestion

# skiena pg 258
# does this rely on min or max of cost function? Right now I'm maximizing (easy to fix)
def accept_neighbor_solution(cost_K, cost_K_prime, T):
  """Probability of accepting a solution with a worse cost""" 
  delta = cost_K - cost_K_prime
  exp_term = -delta / (K * T)
  # how should I calculate e^? ok to use a fast function here?
  if math.exp(exp_term) >= random.random(): 
    return True
  else: return False

# this is the version that accepts non-MIS neighbor solutions
def neighbor_union_subtract(G, K):
  # get index of random node in V
  v = random.randint(0,  Graph.nvertices - 1)
  # print(f'v:{v}')
  # union or subtract v with current subset 
  if v in K.node_set:
    K_prime = (K.node_set).difference({v})
    K_prime_nvertices = K.nvertices - 1
  else:
    K_prime = (K.node_set).union({v})
    K_prime_nvertices = K.nvertices + 1
  return (K_prime, K_prime_nvertices)

def simulated_annealing(G, num_itrs, T, ALPHA):
  t = T
  total_itr = 0
  while(total_itr < MAX_ITR): #FIXME
  #until temp doesn't change OR no solution change OR maxitr  
    # create initial solution, with random subset of Graph G nodes
    K = Subgraph(G, random_subset=True)
    for i in range(ITR_PER_T):
      # construct neighbor via union/subtract
      K_prime, K_prime_nvertices = neighbor_union_subtract(G, K)
      # calculate cost of K and K_prime,
      K_cost = cost(len(K.nvertices, count_edges(G, K.node_set)) 
      K_prime_cost = cost(K_prime_nvertices, count_edges(G, K_prime))
      #FIXME is this 100%?
      if K_cost < best_cost:
        best_cost = K_cost
        best_node_set = K.node_set
      elif K_prime_cost < best_cost:
        best_cost = K_prime_cost
        best_node_set = K_prime_cost
      # relies on MINIMIZING cost function. 
      if K_cost >= K_prime_cost:
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
      # accept worse solution by probability, a_n_s returns t/f
      elif accept_neighbor_solution(K_cost, K_prime_cost, t):
        K.node_set = K_prime
        K.nvertices = K_prime_nvertices
    # reduce temp
    t = ALPHA * t
  return (best_node_set, best_cost, total_itr)

G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
K = Subgraph(G, random_subset=True)
print(G,K)

# note that this section relies on MAXIMIZING cost function
# hardcoded random walk solution, ish
# had a problem -- was staying in the K state unless a neighbor solution was strictly better -- this was preventing exploration. Once I fixed that, things work as expected. 
best_cost = float('-inf')
best_node_set = None

K_cost = cost(K.nvertices, count_edges(G, K.node_set))
for i in range(50):
  print(f'k {K.node_set}, k cost {K_cost}')
  K_prime = get_neighbor_solution(G,K)
  K_prime_cost = cost(len(K_prime), count_edges(G, K_prime))
  
  print(f'k prime {K_prime}, cost {K_prime_cost}')
  if K_prime_cost > best_cost and K_prime_cost > K_cost:
    K.node_set = K_prime
    K.nvertices = len(K_prime)
    best_cost = K_prime_cost
    best_node_set = K_prime
    K_cost = cost(K.nvertices, count_edges(G, K.node_set))
  elif K_cost > best_cost:
    best_cost = K_cost
    best_node_set = K.node_set
  if K_cost == K_prime_cost:
    K.node_set = K_prime
    K.nvertices = len(K_prime)

  #FIXME where does this call go
  K_cost = cost(K.nvertices, count_edges(G, K.node_set))
print(f'best cost {best_cost}')
print(f'set final {K.node_set}')

