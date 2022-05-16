EAMPA Program
#####################################################
Embedded Atom Method Potential Analyser

Written in Python and Fortran (uses F2PY).



Install on Standard Linux Desktop
#####################################################

CD to where you want the python code to be installed.

git clone https://github.com/BenPalmer1983/eampa_v3.git

cd eampa_v3

./install.sh



What does it do?
#####################################################

Can be used to calculate:
· energy
· energy, forces 
· energy, forces, stress
of a collection of atoms

Calculates some other properties:
· a0
· v0
· e0
· b0
· 3 independent elastic constants (Cubic - MSKP)
· 9 independent elastic constants (Orthorhombic - RFKJ)
· shear modulus
· young's modulus
· estimated melting temperature

Can be used to fit potentials to data:
· random only fitting
· simulated annealing
· genetic algorithm
· bin hopping (scipy)
· nelder-mead (scipy)
· conjugate-gradient (scipy)
· bfgs (scipy)


What configuration files can be used?
#####################################################

There is a native format (see examples) similar to an xyz type file.





What potentials?
#####################################################

EAM potentials.

Works with potentials in the following form:


Pair:

V(r)

Density - one or multiple distinct functions:

rho_1(r), rho_2(r), rho_3(r)

Embedding term - one or multiple terms:

F_1(rho_1) + F_2(rho_2) + F_3(rho_3)


The potentials can be tabulated or analytic.  It's easy to add new, custom functions, either in Python or the F2PY module.





RSS Calculation
#####################################################

EFS

Four weights: 
1. config weight - multiplier of overall rss
2. energy weight - multiplier of energy rss
3. force weight - multiplier of force rss
4. stress weight - multiplier of stress rss

How each are calculated
1. energy rss - rss between overall energy of calculated config vs known added up over all configs
2. force rss - rss for all forces, all 3 directions added up over all configurations










