function [ x, costFinal ] = CoBic( Iaccq, paramsConv, paramsResol )
%CoBic (Constrained Biconvex algortihm) minimizes the function the square
%datafitting function under an l_0 constrait. The algorithm is presented in
%the article "New L_2-L_0 algorithm for single-molecule localization
%microscopy" By Arne Bechensteen, Laure Blanc-Féraud and Gilles Aubert.
%   Inputs:
%   Iaccq= the observation
%   ParamsConv = parameters connected to the modelisation of the matrix A.
%   Here A is the compositopn reduction of size matrix (paramsConv.M) and a convolution with the
%   point spread function (paramsConv.A)
%   paramsResol contains the parameter connected to the solving of the
%   problem. 
%   paramsResol.itermax is the maximum of iterations ofthe algorithm
%   paramsResol.stopCrit is the stopping criteria for the 
%   algorithm. This is in the sens of difference of cost
%   paramsResol.itermaxFista is how maximum of iterations of the inner -
%   algorithm FISTA
%   paramsResol.stopCritFista is the stopping criteria for the inner
%   algorithm FISTA
%   paramsResol.kmax is maximum number of non-zero pixels we
%   want to reconstruct.
%   paramsResol.ck and
%   paramsResol.dk are two parameters connected to PAM, to enforce strict
%   convexity. 
%   paramResol.rhot is the intitial rho
%   paramResol.rhotStop is the max value of rho (for which rh)
% Output:
% The reconstructed sparse image x
% Costfunction for each time we increase rho.
% Note that we normalize the initial problem, so all the norms of the
% colummns are 1. 
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


% === Varible definition

dk = paramsResol.dk;
ck= paramsResol.ck;
rhotStop= 2*4*norm(Iaccq(:)); %is the max value of rho will take. 
rhot = paramsResol.rhot;
k=paramsResol.kmax;
% === inititalization of solutions ===
[sz1, sz2] = size(paramsConv.A);
xt = zeros(sz1, sz2);
ut = zeros(sz1, sz2);

% === Mech is sucht that Mech*x*MechT reduces the dimension of x. 
Mech = paramsConv.M; 
MechT = paramsConv.M';
ech = size(paramsConv.A,1)/size(Iaccq,1);

% === Calculate the norm of each column of the matrix A ===
for i = 1:ech
    for j = 1:ech
        matr = Mech * circshift(paramsConv.A, [j-1,i-1]) *MechT;
        normai(i,j) = sqrt(sum(matr(:).^2));
    end
end

normai = repmat(normai, size(Mech,1));
norminv = 1./normai;


Afft = fft2(fftshift(paramsConv.A)); % Gaussian kernel by FFT2

t=1; % count for the number of big iterations (big, each time we update \rho)
notConv = 1;

while notConv
    iter=0;
    cost=[];
    notConvS = 1;
    while notConvS
        iter=iter+1;
        xold=xt;
        % Minimization with respect to x
        [xt, ~] = Fista( Iaccq, paramsConv, paramsResol, xt, ut,rhot); 
        
        %Minimization with respect to u
            
        z =  reshape(ut + (rhot *dk) * xt,1,sz1*sz2);
            %Applying a function BreakPointSeach by Stavanov in "Convec
            %Quadratic minimization subject to a linearr constraint and box
            %constraints"
        ut= reshape( sign(z) .* BreakPointSearch( 0.5*ones(size(z)) ,abs(z), k,zeros(size(z)) ,ones(size(z)) ),sz1,sz2 );
       
        % Calculating the cost
        xitfft = fft2(norminv.*xt);
        cost(iter) = 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iaccq,2)^2 + rhot * (norm(xt(:),1) - (xt(:)'*ut(:)));
        
        %Testing diverse convergence criteria
        if iter > paramsResol.itermax
            notConvS = 0;
        elseif  iter ~= 1
            diffCout = abs(cost(iter) - cost(iter-1))/abs(cost(iter)+eps);
            
            diffX = norm(xt(:)-xold(:))/(norm(xold(:))+eps);
            
            if (diffCout < paramsResol.stopCrit  && diffX < paramsResol.stopCrit ) || iter == paramsResol.itermax
                
                notConvS = 0;
                
            end
        end
    end
    % Update Rho
    rhot = min (rhotStop, 2*rhot);
    
    if rhot==rhotStop
        notConv=0;
    end
    
    
    
    
    % === Cost for each Rho update ===
        xitfft = fft2(norminv.*xt);
        costFinal(t) = 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iaccq,2)^2 + rhot * (norm(xt(:),1) - (xt(:)'*ut(:)));
    t = t + 1;
    
end

%Final x!
x = norminv.*xt;

end

