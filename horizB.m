function G = horizB(i,Azvals,Elvals,lam,ipol)
% G = horizB(i,Azvals,Elvals,lam,ipol);

if (nargin < 5)
    ipol = 1;
end

if (ipol == 1)    % vertical polarization
    G = cosd(Azvals);
    G(Elvals < 0) = 0;
else              % horizontoal polarization--no response
    G = 0;
end

end