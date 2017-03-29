function A = arrayManifold(egain,pol,lam,Rarray,Euler,Azvals,Elvals,...
    normalize)
% A = ...
%    arrayManifold(elemgain,pol,lam,Rarray,Euler,Azvals,Elvals,normalize);
% Computes the total array manifold at specified azimuth and elevation
% directions for given polarization:  vertical (E local vertical), and
% horizontal (E local u, which is cross(k vector, vertical)).
%   El rotates up to z from xy, Az about z right-handed (CCW).
% Source direction is d = [cosd(El)*cosd(Az); cosd(El)*sind(Az); sind(El)];
%   The user provides a function G = elemgain(i,Azvals,Elvals,lam); This 
% function provides the complex gain of element i relative to the gain
% assumed in the SNR in a matrix conformal to Azvals, Elvals.  
% egain     -- Handle for an element gain function @g such that
%              G = g(i,Azvals,Elvals,lam,ipol) must return the 
%              amplitude gains of element i in matrix G with same dims
%              as Azvals, Elvals, when these are in degrees.
%              ipol = 1 --> vertical, 2 --> Horizontal
% pol       -- 2 x 1 complex vector indicating the incident polarization 
%              in terms of vertical and horizontal.  E.G., 
%              pol = [1; 0] for vert, [0; 1] for horizontal, pol = [1; 1i]
%              for RHCP.
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about 
%              +z (Azimuth), then about rotated -y (Elev), then about 
%              rotated-rotated +x (Roll). 
% Azvals    -- nEl x nAz azimuth search angles, constant down each column,
%              angles in degrees
% Elvals    -- nEl x nAz elevation search angles, constant
%              along a row, degrees.  El = pi/2 is straight up along z.
% normalize -- OPTIONAL argument, true ==> produce unit vectors
% A         -- nelem x nEL x nAz complex antenna manifold at each elevation
%              and azimuth, for the given incident polarization

pol = pol / norm(pol);   

[Av,Ah] = ...
    arrayManifoldPol(egain,lam,Rarray,Euler,Azvals,Elvals,normalize);

A = Av*pol(1) + Ah*pol(2);

end

