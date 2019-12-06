# CoBic-and-PeBic-SMLM
Single Molecule localization microscopy code based on a deconvolution and a sparsity parameter L0. The sparsity term is reformulated using an auxiliary variable, and the final costfunction is biconvex. We propose two algorithms  CoBic (Constrained Biconvex) and PeBic (Penalized Biconvex) algorithm on SMLM. CoBic is the reformulated constrained L2_L0 problem and PeBic is the reformulated penalized L2_L0 problem. Both with the positivity constraint. 
More information is in an article that will soon be published. 

The content is two main functions, CoBic.m and PeBic.m. Both of them call the Fista.m algorithm to do the minimization with respect to x. Cobic calls also on BreakSearchPoint.m which solves the quadratic knapsack problem using the algorithm of Stefanov. 
Testfile.m showes an example of CoBic and PeBic, using ParamsAcq to create the deconvolution matrix and reduction matrix. The file uses a simulated acquistion, Testimage. 
