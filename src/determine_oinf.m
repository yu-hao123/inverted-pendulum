function [Sf,bf] = determine_oinf(Af,phi,Gamma,h,Spsi,bpsi,max_iter,tol)
    if nargin < 8, tol = 0; end
    r = length(bpsi);
    SpsiGamma = Spsi*Gamma;
    Spsih = Spsi*h;
    S = SpsiGamma; b = bpsi - Spsih;
    n = size(Af,1);
    flag_redund = 0; i = 1;
    
    while ( (i <= max_iter) && (flag_redund == 0) )
        disp(['Iteration' num2str(i) ' / ' num2str(max_iter)])
        flag_redund = 1;
        SGAi = SpsiGamma*(Af^i);
        aux = eye(n);
        for l = 1:i-1
            aux = aux + (Af^l);
        end
        right_side = bpsi - SpsiGamma*aux*phi - Spsih;
        for j = 1:r
            c = SGAi(j, :)';
            d = right_side(j);
            
            t(j) = redundancy_test(S, b, c, d);
            if t(j) > tol
                flag_redund = 0;
                S = [S; c']; b = [b; d];
            end
        end
        i = i + 1;
    end
    
    if (flag_redund == 0)
        disp('Could not find the larger invariant set')
        Sf = []; bf = [];
    else
        Sf = S; bf = b;
    end