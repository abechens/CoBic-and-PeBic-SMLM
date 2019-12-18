%% Demonstration of algorithm
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%% Parameters to create the observation matrix 
% (We will load it instead of creating it. 

% paramsAcq.ech = 4; % The L factor. How much upsampling
% paramsAcq.N = 64*paramsAcq.ech; %Size of ground truth grid
% paramsAcq.sig = 109.70; % PSF estiamtion (sigma)
% paramsAcq.c = 0; % center of PSF
% paramsAcq.I0 = 1; % intensisty PSF
% paramsAcq.taillePixel = 100/paramsAcq.ech ; % size of each pixel in nm
% 
% %% Creating the Matrix. 
% 
% matriceConv = ParamAcq(paramsAcq);
load('MatriceA.mat') %The matrix A, containing A and M, A is the
%                       convolution matrix, and M ...
%                       the downsampling matrix. (Downsampling done by
%                       M*x*M'). 

%% Defining the parameters to solve to problem
% Converging parameters
paramsResol.itermaxFista= 4000;
paramsResol.stopCritFista = 10^(-5);
paramsResol.itermax = 4000;
paramsResol.stopCrit = 10^(-5);

%Initial parameters
paramsResol.rhot = 10^0;
paramsResol.ck = 10^4;
paramsResol.dk = 10^4 ;

%Sparsity parameters
paramsResol.kmax = 220; % The number of non-zeros pixels reconstructed
paramsResol.lambda=0.05; % penalty parameter

load('testimage.mat') %The test image is the first acquisition of a 361 High
%                       -density acquisition. Simulated, and found at the
%                       page: http://bigwww.epfl.ch/smlm/datasets/index.html

%% Constrained Biconvex algorithm
%[xCoBic, ~] = CoBic( testimage/max(testimage(:)),matriceConv, paramsResol);  
xCoBic= xCoBic*max(testimage(:));
%% Penalized Biconvex algorithm
[xPeBic, ~] = PeBic( testimage/max(testimage(:)),matriceConv, paramsResol);  
xPeBic= xPeBic*max(testimage(:));
