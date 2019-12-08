function [ xfinal ] = BreakPointSearch( mj,s,alp,aj, bj)
% This program solves the following QP:
% min_x x'*diag(mj)*x+s'x
% s.t. sum(d*x) =< alp, aj<=x<=bj
% The algorithm is detailed in the article "Convec Quadratic minimization
% subject to a linear constraint and box constraints" by Stefan Stefanov. 

%Step 0. initialisations 
d=ones(1,length(mj)); %in our case d=1

J = logical(ones(1,length(mj)));
k = 0;
Ja= logical(zeros(1,length(mj)));
Jb= logical(zeros(1,length(mj)));
indiceG= logical(zeros(1,length(mj)));
temp1 = (-2*mj.*aj+s)./d;
temp2 = (-2*mj.*bj+s)./d;
JLa =  0>=temp1;
JLb = 0<=temp2;
JL= logical(ones(1,length(mj))-JLa - JLb);

alpk = alp;
n=length(mj);

%Step 1. Calculate delta 0
dzero = sum(d(JLa).*aj(JLa))+sum(d(JLb).*bj(JLb))+ 0.5 * sum(d(JL) .* s(JL) ./ mj(JL))-alpk;


if dzero <= 0 % then step 8
    x(JLa) = aj(JLa);
    x(JLb) = bj(JLb);
    x(JL) = 0.5 * d(JL).* s(JL)./mj(JL);
else %Step 2
    
    JLa = Ja;
    JLb = Jb;
    while 1 %Step3
        
        JLs = logical(max(J-indiceG,0));
        lambdat = sum (d(JLs).^2./mj(JLs))^-1 * ( 2 * sum(d(JLa) .* aj(JLa)) + 2 * sum (d(JLb) .* bj(JLb)) + sum(d(JLs) .* s(JLs) ./ mj(JLs))-2* alpk );
        temp1 = (-2*mj.*aj+s)./d;
        temp2 = (-2*mj.*bj+s)./d;
        JLa =  logical(max((lambdat>=temp1) - indiceG,0));
        JLb = logical(max((lambdat<=temp2) - indiceG,0));
        JL= logical(ones(1,length(mj))-JLa - JLb - indiceG);
        %Step 4
        dk = sum(d(JLa).*aj(JLa))+sum(d(JLb).*bj(JLb))+ 0.5 * sum(d(JL) .* s(JL) ./ mj(JL)) -0.5 * lambdat * sum (d(JL).^2./mj(JL))-alpk;
        % Step 5
        if dk== 0 || sum(JLs)==0 || sum(JLa)+sum(JLb)==0
            Ja=logical(min(JLa+Ja,1));
            Jb=logical(min(JLb+Jb,1));
            J=JL;
            x(Ja) = aj(Ja);
            x(Jb) = bj(Jb);
            x(J) = 0.5 *(s(J)-lambdat*d(J))./mj(J);
            break % Step 8 
        elseif dk>0 %Step 6
            x(JLa) = aj(JLa);
            alpk = alpk - sum (d(JLa).*aj(JLa));
            indiceG = indiceG +JLa;
            J=logical(max(J-indiceG,0));
            n = n  - sum(JLa);
            Ja=logical(min(JLa+Ja,1));
        else % Step 7
            x(JLb) = bj(JLb);
            alpk = alpk - sum (d(JLb).*bj(JLb));
            indiceG = indiceG +JLb;
            J=logical(max(J-indiceG,0));
            n = n  - sum(JLb);
            Jb=logical(min(JLb+Jb,1));
        end
    end
end
xfinal= x; %Step 10
end
