from graph import * 
import math

#Global Parameters
LAMBDA = 1          #should be a pos integer constant (Leo, p 870)
TRIALS_PER_T = 10
K = 1 #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256
T = 1

# cost function notation from Feo paper (greedy alg vs SA)
# this cost should be maximized (objective function)
def cost(num_nodes_K, num_edges_K):
  """Calculate cost of solution, given |K| and |E_K| where K is a (proper? why proper) subset of V and E_K is set of edges induced by K"""
  return num_nodes_K - (LAMBDA * num_edges_K) #Feo paper version
  #return num_K - (LAMBDA * E_K)/T #Skiena's cost function suggestion

# skiena pg 258
def accept_neighbor_solution(cost_K, cost_K_prime, T):
  """Probability of accepting a solution with a worse cost""" 
  delta = cost_K - cost_K_prime
  exp_term = -delta / (K * T)
  # how should I calculate e^? ok to use a fast function here?
  if math.exp(exp_term) >= random.random(): 
    return True
  else: return False

# this is the version that accepts non-MIS neighbor solutions
def get_neighbor_solution(Subgraph):
  # get index of random node in V
  v = random.randint(0, Subgraph.nvertices - 1)
  # union or subtract v with current subset 
  # FIXME make version that doesn't copy sets
  if v in Subgraph.node_set:
    K_prime = (Subgraph.node_set).difference({v})
  else: K_prime = (Subgraph.node_set).union({v})
  return K_prime

G = GraphAL(5, [(1,2),(2,3),(3,4),(0,4),(1,4)])
K = Subgraph(G, random_subset=True)
print(G,K)


# hardcoded random walk solution

best_cost = float('-inf')
best_node_set = None

K_cost = cost(K.nvertices, count_edges(G, K.node_set))
for i in range(5):
  print(f'k {K.node_set}, k cost {K_cost}')
  K_prime = get_neighbor_solution(K)
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
  #FIXME where does this call go
  K_cost = cost(K.nvertices, count_edges(G, K.node_set))
print(f'best cost {best_cost}')
print(f'set final {K.node_set}')

