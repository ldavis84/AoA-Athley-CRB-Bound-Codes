function [i,xsi] = interpHelper(x,xvec)
% [i,xsi] = interpHelper(x,xvec);
% An interpolation helper.  If x is a set of request values, and xvec a 
% vector of monotonically increasing values, it finds i and xsi such that
% x = xvec(i)*(1-xsi) + xvec(i+1)*xsi;
% xsi is in [0 1] unless an extrapolation is being performed (x not in the
% range).
% x     -- n x 1 values for which to obtain the reference
% xvec  -- vector of reference values

xvec = xvec(:);

n = length(xvec);
x = x(:);

nx = length(x);
i = zeros(nx,1);

extrapleft = x < xvec(1);
i(extrapleft) = 1;

extrapright = x >= xvec(n);
i(extrapright) = n-1;

interp = ~extrapleft & ~extrapright;
xinterp = x(interp);
ninterp = length(xinterp);

kinterp = find(interp);

xspan = xvec(n) - xvec(1);

if (ninterp > 0)
    
    notdone = true(ninterp,1);
    
    ilow = ones(ninterp,1);
    ihi = n*ilow;
    
    while (any(notdone))
        
        knotdone = kinterp(notdone);
        
        xi = xinterp(notdone);
        il = ilow(notdone);
        ih = ihi(notdone);
        
        nu = (xi - xvec(il)) ./ (xvec(ih) - xvec(il));
        
        iest = il + nu.*(ih - il);  % MUST be between il and ih
        
        % check very close to a reference point
        
        itmp = round(iest);
        equal = abs(xi - xvec(itmp)) <= 2*eps*xspan;
        kassign = knotdone(equal);
        i(kassign) = itmp(equal);
        
        % check other cases
        
        itmp = floor(iest);
        
        % This interval is ok.
        
        thisinterval = (xi > xvec(itmp) & xi < xvec(itmp+1)) & ~equal;
        kassign = knotdone(thisinterval);
        i(kassign) = itmp(thisinterval);
        
        % The true interval is below this one
        
        below = xi < xvec(itmp) & ~equal;
        ih(below) = itmp(below);
        ihi(notdone) = ih;
        
        % The true interval is above this one
        
        above = xi >= xvec(itmp+1) & ~equal;
        il(above) = itmp(above) + 1;
        ilow(notdone) = il;
       
        notdone(thisinterval | equal) = false;

    end
end

xsi = (x - xvec(i)) ./ (xvec(i+1) - xvec(i));

end
