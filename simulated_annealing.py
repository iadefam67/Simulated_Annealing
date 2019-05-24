from graph import * 
import math
#Global Parameters
LAMBDA = 1
TRIALS_PER_T = 10
K = 1 #constant used by ap function to "to normalize the cost function so that almost all transitions are accepted at the starting temp", skiena p 256


G = GraphAL(4, [(2,3)])
print(G)

# cost function notation from Feo paper (greedy alg vs SA)
def cost(num_K, num_E_K, LAMBDA, T):
  """Calculate cost of solution, given |K| and |E_K| where K is a proper subset of V and E_K is set of edges induced by K"""
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


