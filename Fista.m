function [xopt, cost ] = Fista( Iacq, paramsConv, paramsResol, xit, ut,rhot)
%Fista minimizes, for a given rhos and u, the functional
%G_rho(x,u)+1/(2ck)||x-xk||^2 using the FISTA algorithm
%   Inputs: 
%   Iacq  is the observation 
%   ParamsConv = parameters connected to the modelisation of the matrix A.
%   Here A is the compositopn reduction of size matrix (paramsConv.M) and a convolution with the
%   point spread function (paramsConv.A)
%   paramsResol contains the parameter connected to the solving of the
%   problem. 
%   paramsResol.itermaxFista is how maximum of iterations of the inner -
%   algorithm FISTA
%   paramsResol.stopCritFista is the stopping criteria for the inner
%   algorithm FISTA
%   Note that we are not using rho defined from paramsResol, since that is
%   the initialization.
%Output: 
%   xopt is the the argmin
% cost is the costfunction evolution
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


% -- Set default values

if ~isfield(paramsResol,'itermaxFista'), paramsResol.itermaxFista=10000; end
if ~isfield(paramsResol,'stopCritFista '),paramsResol.stopCritFista  = 1e-5; end
if ~isfield(paramsConv,'A') || ~isfield(paramsConv,'M') , error('No observation matrix defined'); end

% === Varible definition

ck=paramsResol.ck;

Afft = fft2(fftshift(paramsConv.A)); 
Aconj = conj(Afft);
Mech = paramsConv.M;
MechT = paramsConv.M';
zt = xit + rhot * ck * ut;
ech = size(Afft,1)/size(Iacq,1);

% === norm of each column of A ===
for i = 1:ech
    for j = 1:ech
        matr = Mech * circshift(paramsConv.A, [j-1,i-1]) *MechT;
        normai(i,j) = sqrt(sum(matr(:).^2));
    end
end

normai = repmat(normai, size(Mech,1));
norminv = 1./normai;

% === Lipchitz constant  ===
normA= max(norm(Mech*MechT)*abs(Afft(:)))*max(norminv(:));
lipCst =normA^2 + 1/ck +1;
gam = 1/lipCst;

% === Fista ===
nonConv=1;
y = xit;
xt=xit;

cost = [];

tk = 1;
iter = 1;

while nonConv
    xold = xt;
    told = tk;
    %Gradient calculations
    yfft = fft2(norminv.*y);
    gradientF = norminv.*( real(ifft2(Aconj.*fft2((MechT*(Mech * real(ifft2(Afft .* yfft)) * MechT - Iacq)*Mech)))) ) + 1/ck * (y-zt);
    yt = y - gam*gradientF;
    %Proximal step
    xt = max( ( yt - rhot * gam),0);
    %Fista steps
    tk = 0.5*(1+sqrt(1+4*told^2));
    temp1 = (told-1)/tk;
    y =( xt + temp1 * (xt - xold));
    
    
    % === Cost calculations  ===
    
    xitfft = fft2(norminv.*xt);
    cost(iter) = 0.5*norm(Mech*real(ifft2(Afft.*xitfft))*MechT-Iacq,'fro')^2 + rhot * norm(xt(:),1) + 1/ck * 0.5 *norm(xt(:)-zt(:))^2 ;
   
    
    % === Convergence  ===
    
    if iter > paramsResol.itermaxFista
        nonConv = 0;
    elseif  iter ~= 1
        diffCout = abs(cost(iter) - cost(iter-1))/abs(cost(iter)+eps);
       
        diffX = norm(xt(:)-xold(:))/(norm(xold(:))+eps);
        
        if (diffCout < paramsResol.stopCritFista  && diffX < paramsResol.stopCritFista ) || iter == paramsResol.itermaxFista
         
            nonConv = 0;
           
        end
    end
    
    iter = iter + 1;
end

xopt=xt;

end

