function G = vertdipole(i,Azvals,Elvals,lam,ipol)
% G = vertdipole(i,Azvals,Elvals,lam,ipol);

if (nargin < 5)
    ipol = 1;
end

if (ipol == 1)    % vertical polarization
    G = cosd(Elvals);
%     G(Elvals < 0) = 0;
else              % horizontoal polarization--no response
    G = eps;
end

end