function [ x, costFinal ] = PeBic( Iaccq, paramsConv, paramsResol )
%PeBic (Penalized Biconvex algortihm) minimizes the function the square
%datafitting function penalized with the l_0 norm. The algorithm is presented in
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
%   paramsResol.lambda is the penalization paramater linked to L_0 norm
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

lambdat = paramsResol.lambda;
rhotStop= 2*4*norm(Iaccq(:)); %is the max value of rho will take. 
rhot = paramsResol.rhot;
dk = paramsResol.dk;


% === init of solutions ===
[sz1 sz2] = size(paramsConv.A);
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



Afft = fft2(fftshift(paramsConv.A)); % THe Psf in fourier.
Mech = paramsConv.M;
MechT = paramsConv.M';


% count for the number of big iterations (big, each time we update \rho)


t=1;
notConv = 1;
while notConv %This loop is for each time we update rho 
    
    iter=0;
    notConvS = 1;
    
    while notConvS % This is the algorithm Pam
        iter=iter+1;
        xold=xt;
        %Minimization with respect to x
        [xt, ~] = Fista( Iaccq, paramsConv, paramsResol, xt, ut, rhot);
        %Minimization with respect to u.
        yt = ut + rhot *dk .* xt ;
        ut =sign(yt) .* (min(max( (abs(yt) - lambdat*dk),0),1));
        
        %cost Function
        xitfft = fft2(norminv.*xt);
        cost(iter) = 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iaccq,2)^2 +lambdat * norm(ut(:),1)+ rhot * (norm(xt(:),1) - (xt(:)'*ut(:)));
        
        %Testing diverse convergence criteria
        if iter > paramsResol.itermax
            notConvS = 0;
        elseif  iter ~= 1
            diffCout = abs(cost(iter) - cost(iter-1))/abs(cost(iter)+eps);
            
            diffX = norm(xt(:)-xold(:))/(norm(xold(:))+eps);
            
            if (diffCout < paramsResol.stopCrit  && diffX < paramsResol.stopCrit) || iter == paramsResol.itermax
                
                notConvS = 0;
                
            end
        end
    end
    %Updating Rho
    rhot = min (rhotStop, 2*rhot);
    
    if rhot==rhotStop
        notConv=0;
    end
    
    
    
    
    % === Global Cost  ===
    if all(sum(abs(ut)>1))==0
        xitfft = fft2(norminv.*xt);
        costFinal(t) = 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iaccq,2)^2 +lambdat * norm(ut(:),1) + rhot * (norm(xt(:),1) - (xt(:)'*ut(:)));
    else
        cost(t) = inf;
    end
   
    t = t + 1;

end
%Reconstructing x
x = norminv.*xt;
end

