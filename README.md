# SemiClassicalMC
Semi-classical Monte Carlo solver for spin-valley models of the form
$$ H = \sum_{ij} \left(\sum_{\mu,\nu,\kappa,\lambda = x, z, y}(\sigma_i^\mu \otimes \tau_i^\kappa)\ J_{s, ij}^{\mu\nu}\ J^{\kappa\lambda}_{v, ij}\ (\sigma_j^\nu \otimes \tau_j^\lambda) 
+ \sum_{\mu,\nu = x, y, z} \sigma^\mu_i\ K^{\mu\nu}_{s, ij}\ \sigma^\nu_j
+ \sum_{\kappa,\lambda = x, y, z} \tau^\kappa_i\ K^{\kappa\lambda}_{v, ij}\ \tau^\lambda_j
 \right)$$
 in the fundamental (quarter-filling, local Hilberspace dimension d = 4) and self-conjugate (half-filling, d = 6) representation.  The on-site exchange matrices $(i = j)$ need to be diagonal, all other interactions can be completely off-diagonal.
Implemented are both simulated annealing and finite temperature Monte Carlo.
## Example
See Example notebook in folder "example"