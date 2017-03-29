function maxQP = maxdef(Q,P)
% maxQP = maxdef(Q,P);
% Finds the smallest matrix minQP in sense of definiteness such that
% minQP >= Q, minQP >= P.  Does so via simultaneous diagonalization.
% Q          -- n x n non-negative definite matrix
% P          -- n x n second matrix P >= 0
% maxQP      -- desired minimum

[T,sigQ,sigP] = simulDiag(Q,P);
Tinv = inv(T);
maxQP = Tinv' * diag(max(sigQ,sigP)) * Tinv;

end