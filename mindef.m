function minQP = mindef(Q,P)
% minQP = mindef(Q,P);
% Finds the largest matrix minQP in sense of definiteness such that
% minQP <= Q, minQP <= P.  Does so via simultaneous diagonalization.
% Q          -- n x n non-negative definite matrix
% P          -- n x n second matrix P >= 0
% minQP      -- desired minimum

[T,sigQ,sigP] = simulDiag(Q,P);
Tinv = inv(T);
minQP = Tinv' * diag(min(sigQ,sigP)) * Tinv;

end