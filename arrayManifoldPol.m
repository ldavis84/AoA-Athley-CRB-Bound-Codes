function [Av,Ah] = ...
    arrayManifoldPol(egain,lam,Rarray,Euler,Azvals,Elvals,normalize)
% [Av,Ah] = ...
%    arrayManifoldPol(eain,lam,Rarray,Euler,Azvals,Elvals,normalize);
% Computes the total array manifold at specified azimuth and elevation
% directions for both polarizations:  vertical (E local vertical), and
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
%              ipol = 1 --> vertical (Z Electric Field), 2 --> Horizontal
%              polarization (+ electric field in cross(k, vert) direction 
%              where k is the direction of propagation
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
% [Av,Ah]   -- nelem x nEL x nAz complex antenna manifold at each elevation
%              and azimuth, Av for vertical pol, Ah for horizontal

normalize = exist('normalize','var') && normalize;

nelem = size(Rarray,2);
[nEl,nAz] = size(Azvals);
Azvals = Azvals(:)';
Elvals = Elvals(:)';

Av = zeros(nelem,nEl*nAz);
Ah = Av;

% Build direction vector, local vertical v, and local horizontal u in 
% Global coordinates.

d = [cosd(Elvals) .* cosd(Azvals); cosd(Elvals) .* sind(Azvals); ...
    sind(Elvals)];

u = [-sind(Azvals); cosd(Azvals); zeros(1,nEl*nAz)];

v = [-sind(Elvals) .* cosd(Azvals); -sind(Elvals) .* sind(Azvals); ...
    cosd(Elvals)];

kwav = 2*pi/lam;

for i = 1:nelem
    
    % Compute local azimuth and elevation values for this aperture.
    
    Qga = euler2Q(Euler(1,i),Euler(2,i),Euler(3,i));
    da = Qga' * d;
    Aza = atan2d(da(2,:),da(1,:));
    Ela = atan2d(da(3,:),sqrt(da(1,:).^2 + da(2,:).^2));
    
    % Compute the local horizontal reference u and vertical v in global
    % coordinates.
    
    horiz = Qga * [-sind(Aza); cosd(Aza); zeros(1,nAz*nEl)];
    vert = Qga * [-sind(Ela).*cosd(Aza); -sind(Ela).*sind(Aza); cosd(Ela)];
    
    % Compute the response of the aperture to its own local vertical and
    % horizontal polarizations, including array factor phasor.
    
    phasor = exp(1i*kwav*(Rarray(:,i)'*d));
    
    Ava = phasor .* egain(i,Aza,Ela,lam,1);
    Aha = phasor .* egain(i,Aza,Ela,lam,2);
    
    Av(i,:) = sum(vert.*v,1) .* Ava + sum(horiz.*v,1) .* Aha;
    Ah(i,:) = sum(vert.*u,1) .* Ava + sum(horiz.*u,1) .* Aha;
end

if (normalize)
    
    magv = sqrt(sum(abs(Av).^2,1));
    magh = sqrt(sum(abs(Ah).^2,1));

    Av = Av ./ repmat(magv,nelem,1);
    Ah = Ah ./ repmat(magh,nelem,1);
    
    iz = magv <= 100*eps*(magv + magh);     % zero out null elements
    Av(:,iz) = zeros(size(Av(:,iz)));
    
    iz = magh <= 100*eps*(magv + magh);
    Ah(:,iz) = zeros(size(Ah(:,iz)));
   
end

Av = reshape(Av,nelem,nEl,nAz);
Ah = reshape(Ah,nelem,nEl,nAz);

end

