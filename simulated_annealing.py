from graph import * 
import math

#Global Parameters
LAMBDA = 1          #should be a pos integer constant (Leo, p 870)
TRIALS_PER_T = 10
K = 1 #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256

G = GraphAL(4, [(2,3)])

cost = float('inf')
K = None
nvertices_K = None
nedges_K = None
print(G)

# cost function notation from Feo paper (greedy alg vs SA)
def cost(num_K, num_E_K, LAMBDA, T):
  """Calculate cost of solution, given |K| and |E_K| where K is a (proper? why proper) subset of V and E_K is set of edges induced by K"""
  return num_K - LAMBDA * num_E #Feo paper version
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

def get_neighbor_solution(Graph):
  # get index of random node in V
  v = random.randint(0, Graph.nvertices - 1)
  # make a copy of the subgraph nodes
  K_prime = {node for node in Graph.subgraph_K}
  # union or subtract it with current subset 
  if v in K_prime:
    K_prime.remove(v)
  else K_prime.add(v)
  # calculate the cost of neighbor K' solution:
  




