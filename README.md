# Decoupled SDP for 2-D semidifinite programming

Paper:

* Z. Zhang, Y. Wang, and Z. Tian, "Efficient Two-Dimensional Line Spectrum Estimation Based on Decoupled Atomic Norm Minimization," Signal Processing.
* Z. Tian, Z. Zhang and Y. Wang, Low-complexity optimization for Two-Dimensional Direction-of-arrival Estimation via Decoupled Atomic Norm Minimizationg'', in 42th International Conference on Acoustics, Speech, and Signal Processing (ICASSP 2017), 2017.

This is an efficient optimization technique for super-resolution 2-D harmonic retrieval, named decoupled atomic norm minimization (D-ANM). 

State-of-the-art 2-D ANM approach vectorizes the 2-D measurement to a 1-D equivalence, which incurs huge computational cost in 2-D processing, making it infeasible in real applications. We develop a novel decoupled approach of 2-D ANM via semi-definite programming (SDP), which introduces a new atom set to naturally decouple the joint observation in both dimensions without loss of optimality. Accordingly, the original large-scale 2-D problem is equivalently reformulated via two decoupled one-level Toeplitz matrices, which can be solved by simple 1-D frequency estimation with automatic pairing. 

Compared with the conventional vectorized approach, the proposed technique reduces the computational complexity from the order of O(N^3.5 M^3.5) to O((N +M)^3.5), where N and M are the problem size on the two dimensions respectively. 

It also retains the benefits of ANM in terms of super-resolution, small number of required measurements and robustness to source correlation. 

# Files

* d2dsdp.m

The proposed decoupled SDP.

* vsdp.m

The conventional vectorized SDP.

