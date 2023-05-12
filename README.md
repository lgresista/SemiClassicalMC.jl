# SemiClassicalMC
This package allows the semi-classical analysis of $\mathfrak{su}(4)$ spin models using Metropolis Monte Carlo and simulated annealing algorithms. It can treat models with Hamiltonians of the form

$$H = \sum_{ij} \sum_{\mu,\nu,\kappa,\lambda = 0}^3(\sigma_i^\mu \otimes \tau_i^\kappa)\ J_{ij}^{\mu\nu\kappa\lambda} (\sigma_j^\nu \otimes \tau_j^\lambda) \equiv  \sum_{ij} \sum_{a, b = 0}^{15} T^a_i\ J_{ij}^{ab}\ T^b_j.$$

The spin-valley operators $T^a \equiv \sigma_i^\mu \otimes \tau_i^\kappa$ are generators of the Lie-Group SU(4) and defined via a fermionic parton decomposition as

$$\sigma^\mu_i \otimes \tau^\nu_i = P_n f^\dagger_{isl} \theta^{\mu}_ {ss'} \theta^\nu _{ll'} f _{is'l'},$$

where $\theta^\mu$ are the Pauli matrices (with $\theta^0 = 1$) and $P_n$ is the projector on the subspace of $n$ partons per site, i.e. exactly enforcing $f^\dagger_{isl}f_{isl} = n = 1, 2$. The coupling matrix $J$ can be completely off-diagonal, except for on-site interactions ($i = j$), where we demand it to be diagonal in the sense that $J_{ii}^{\mu\nu\kappa\lambda} \sim \delta_{\mu,\nu} \delta_{\kappa,\lambda}$ or $J_{ii}^{a,b}\sim \delta_{a,b}$ The Hamiltonian can be defined on arbitrary lattice graphs using the [LatticePhysics.jl](https://github.com/janattig/LatticePhysics.jl) package. 

By "semi-classical" we refer to the limit of no entanglement between different lattice sites, which is enforced by making a product-state ansatz for the wave function. In this limit, finite-temperature observables can by efficiently calculated by sampling the state of product-state wave functions via Monte Carlo algorithms. The package implements a standard Metropolis Monte Carlo with adaptive step-size for finite temperature calculations, and a minimization scheme using simulated annealing with subsequent stochastic gradient descent to efficiently obtain the semi-classical ground state. Details on the algorithms are described in the appendix of [our preprint](https://doi.org/10.48550/arXiv.2303.01244).

# Instalation
Before installing the package, first the unregistered dependency [LatticePhysics.jl](https://github.com/janattig/LatticePhysics.jl) has to be installed by following the instructions on its github page. Afterwards, the SemiClassicalMC package can be installed using the Julia package manager by switching to package mode in the REPL (with `]`) and writing
```julia
pkg> add https://github.com/lgresista/SemiClassicalMC.git
```
# Usage
For several examples of the packages capabilities and how to use them, see the Jupyter Notebook in `docs/examples.ipynb`.