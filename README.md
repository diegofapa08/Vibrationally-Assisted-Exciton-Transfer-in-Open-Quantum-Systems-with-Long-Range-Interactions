# Vibrationally-Assisted-Exciton-Transfer-in-Open-Quantum-Systems-with-Long-Range-Interactions
Sample codes for numerical simulations performed in the paper: Vibrationally Assisted Exciton Transfer in Open Quantum Systems with Long-Range Interactions

The code DimerEvolution.m can be used to simulate the dynamics of the two-monomer (dimer) case. This code can be used to generate Figs 2, 3, and 4 of the paper.

The code threeqMonomerEvolution.m is an extension of the dimer code with increases the number of qubits per site from 2 to 3. Similar codes can be created to further increase the number of qubits per site to generate plots similar to those in Fig. 5 (a) and (b).

The code TrimerEvolution.m is an extension of the dimer code with increases the number of sites from 2 to 3 while keeping two qubits per site. Similar codes can be created to further increase the number of sites to generate plots similar to those in Fig. 5 (c) and (d).

The code TimeEvolDimerSingleTrjv2.m is a different version of the code DimerEvolution.m but using matrix exponentiation instead of ODE solvers to solve for the dynamics. It requires the use of a function expv(t,A,V) that computes numerically exp(tA)v, where t is time, A is a matrix, and vis  a vector. I highly recommend using the EXPOKIT expv.m code that can be found at https://www.maths.uq.edu.au/expokit/
