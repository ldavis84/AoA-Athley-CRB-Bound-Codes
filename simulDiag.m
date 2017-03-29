function [T,sigQ,sigP,iinull] = simulDiag(Q,P,toler)
% [T,sigQ,sigP,iinull] = simulDiag(Q,P,toler);
% Simultaneously diagonalizes the non-negative definite matrices Q, and P,
% robustly, such that T'*Q*T = diag(sigQ), T'*P*T = diag(sigP), and
% arranged in descending ratio of sigQ ./ sigP, with values (iinull index)
% for both sigQ and sigP zero at the end, if any.
% Q      -- n x n input matrix Q = Q' >= 0 in terms of definiteness
% P      -- n x n input matrix P = P' >= 0
% toler  -- tolerance for small eigenvalues
% T      -- diagonalizing matrix, invertible but not orthogonal
% sigQ   -- numerator diagonalized values Q, order in descending ratio of 
%           sigQ ./ sigP, with both zeros at end.  
% sigP   -- denominator diagonalized values of P, same order as Q
% iinull -- indices of columns of T, and values in Q and P which are both
%           zero.  T(:,iinull) is basis for a common null space

n = size(Q,1);

if (~exist('toler','var'))
    toler = 4*eps;
end

[U,sigQ] = svd(Q); 
sigQ = diag(sigQ);

iigood = sigQ >= sigQ(1)*toler;
sigQ(~iigood) = 0;

U1 = U(:,iigood) ./ (ones(n,1)*sqrt(sigQ(iigood)'));

if (all(iigood))     % Q has no null space.  Easy case.
    
    Pbar = U1'*P*U1;
    sigQ = ones(n,1);
    [V,sigP] = svd(Pbar);
    sigP = diag(sigP);
    
    T = U1 * V;
    
    iiorder = n:-1:1;
    sigP = sigP(iiorder);
    T = T(:,iiorder);
    iinull = [];
    
else  % Q has a null space, P might also
    
    sigQ(iigood) = 1;
    
    ngood = sum(iigood);
    nbad = n - ngood;
    U2 = U(:,~iigood);
    
    T1 = U1;
    T2 = U2;
    
    Pb11 = U1'*P*U1;
    Pb12 = U1'*P*U2;
    Pb22 = U2'*P*U2;
    
    Ptoler = max(real(diag(P)))*toler;
    
    % diagonalize Pb22, pick out/zero nulls
    
    [V2,gam2] = svd(Pb22);   % lies in Q null space
    gam2 = diag(gam2);
    
    jjgood = gam2 >= Ptoler;
    gam2(~jjgood) = 0;
    mgood = sum(jjgood);
    mbad = nbad - mgood;
    jjorder = [(mgood:-1:1) (mgood+1):nbad];
    
    V2 = V2(:,jjorder);
    gam2 = gam2(jjorder);
    T2 = T2 * V2;
    
    iinull = ~iigood;
    iinull(~iigood) = ~jjgood;
    
    % off diagonal stuff
    
    if (any(jjgood))
        
        Pb12 = Pb12 * V2;
        Pb12(:,~jjgood) = 0;
        M12 = Pb12(:,jjgood) ./ (ones(ngood,1)*sqrt(gam2(jjgood)'));
        T12tmp = [-Pb12(:,jjgood) ./ (ones(ngood,1)*gam2(jjgood)') ...
            zeros(ngood,mbad)];
        T1 = T1 + T2 * T12tmp';
        Pb11 = Pb11 - M12*M12';
    end
    
    % Now diagonalize the top
    
    [V1,gam1] = svd(Pb11);
    gam1 = diag(gam1);
    
    kkgood = gam1 >= gam1(1)*toler;
    gam1(~kkgood) = 0;
    
    iiorder = length(gam1):-1:1;
    gam1 = gam1(iiorder);
    V1 = V1(:,iiorder);
    
    T1 = T1 * V1;
    
    T = [T1 T2];
    sigP = [gam1; gam2];
    
end

Tmax = max(abs(T),[],1);    % get decent scaling for T
T = T ./ repmat(Tmax,n,1);

sigQ = sigQ ./ Tmax(:).^2;
sigP = sigP ./ Tmax(:).^2;

end