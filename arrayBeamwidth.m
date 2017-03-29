function [BWAz,BWEl] = arrayBeamwidth(egain,lam,Rarray,Euler,Az,El,Pol,SrchRange)
% [BWAz,BWEl] = arrayBeamwidth(egain,lam,Rarray,Euler,Az,El,Pol,SrchRange);
% Computes array 3 dB beamwidth using manifold and second derivatives.
%    El rotates up to z from xy, Az about z right-handed (CCW).
% Source direction is d = [cos(El)*cos(Az); cos(El)*sin(Az); sin(El)];
%   The user provides a function G = elemgain(i,Azvals,Elvals,lam); This 
% function provides the complex gain of element i relative to the gain
% assumed in the SNR in a matrix conformal to Azvals, Elvals.  A single 
% polarization is  assumed.  
% LINEAR arrays are supported by providing Rarray as a 1-dimensional row
% vector of the X positions of the elements, with Az = 0 down the X axis
% and Az = pi/2 as broadside to the array.
% elemgain  -- function returning a single element's complex gain for each
%              value of azimuth and elevation.  If g() is
%              the function, one passes @g in the call to this function and
%              G = g(i,Azvals,Elvals,lam) must return the complex
%              amplitude gain of element i in matrix G with same dimensions
%              as Azvals, Elvals, when these are in radians.
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements.  For LINEAR
%              array, 1 x nelem X coordinates of elements  
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about
%              +z (Azimuth), then about rotated +x (Elev), then
%              about rotated-rotated +y (Roll). If empty, then all
%              rotations are assumed zero.  If only 1 angle provided, it is
%              assumed to be an azimuth rotation, and the rest are zeros.
% Az,El     -- nominal pointing direction, deg, for LINEAR array the EL
%              argument is ignored and El = 0 is assumed
% Pol       -- 1 x 2 true polarization parameters p1 and p2, deg
% SrchRange -- [Azmin Azmax Elmin Elmax] OPTIONAL limitation on search
%              range for points, pi * [-2 2 -1/2 1/2] default.  For LINEAR
%              only the first two need be provided or are relevant.  In
%              searching for 3 dB points search will not proceed past these
%              boundaries so that beamwdiths may be appropriately truncated
% BWAz      -- azimuth full width at half maximum (3 dB) beamwdith, deg
% BWEl      -- elevation full width at half maximum (3 dB) beamwdith, deg

[ndim,nelem] = size(Rarray);

ERROR = ndim ~= 1 && ndim ~= 3;

if ERROR
    fprintf('\narrayBeamwidth:  Rarray dimensions are %d x %d!\n',...
        ndim,nelem);
    BWAz = [];
    BWEl = [];
    return;
end

LINEAR = ndim == 1;

if LINEAR
    El = 0;
    Rarray = [Rarray; zeros(2,nelem)];
    BWEl = [];
end

z = pol2jones(Pol);

Ao = arrayManifold(egain,z,lam,Rarray,Euler,Az,El,true);   % normalized
U = Ao*Ao';

[~,~,hess] = AoAderiv(egain,z,lam,Rarray,Euler,Az,El,U);

H = -hess;     % gain = 1 - 0.5*u'*H*u

BWAz = 2 / sqrt(H(1,1));  % starting point of search

toler = 1e-2;   % tolerance on 3 dB gain

godB = -10*log10(2);   % desired gain

% Search limits

if (nargin >= 6)
    
    if (LINEAR)
        Azmin = SrchRange(1);
        Azmax = SrchRange(2);
    else
        Azmin = SrchRange(1);
        Azmax = SrchRange(2);
        Elmin = SrchRange(3);
        Elmax = SrchRange(4);
    end
    
else
    Azmin = -360;
    Azmax = 360;
    Elmin = -90;
    Elmax = 90;
end

%--------------------------------------------------------------------------
% Find 3 dB gain points to the left
%--------------------------------------------------------------------------

Azhigh = Az;
gdB_high = 0;

dAz = BWAz / 10;
Aztry = max(Az - BWAz/2,Azmin);

while (1)   % find a point < 3 dB gain
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Aztry,El,true);
    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (Aztry == Azmin)
        gdB = godB;
    end
    
    if (gdB <= godB)
        Azlow = max(Aztry,Azmin);
        gdB_low = gdB;
        break;
    end
    
    Azhigh = Aztry;
    gdB_high = gdB;
    
    Aztry = max(Aztry - dAz,Azmin);
    
end

% Now zoom in...

while (gdB ~= godB && abs(Azhigh - Azlow) >= BWAz*1e-3)
    
    Aztry = (godB*(Azhigh - Azlow) ...
        + Azlow*gdB_high - Azhigh*gdB_low) / (gdB_high - gdB_low);
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Aztry,El,true);    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (gdB > godB + toler)
        Azhigh = Aztry;
        gdB_high = gdB;
    elseif (gdB < godB - toler)
        Azlow = Aztry;
        gdB_low = gdB;
    else
        break;
    end

end

Azleft = Aztry;

%--------------------------------------------------------------------------
% Find 3 dB gain points to the right
%--------------------------------------------------------------------------

Azhigh = Az;
gdB_high = 0;

Aztry = min(max(2*Az - Azleft,Az + dAz),Azmax);

while (1)   % find a point < 3 dB gain
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Aztry,El,true);
    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (Aztry == Azmax)
        gdB = godB;
    end
    
    if (gdB <= godB)
        Azlow = Aztry;
        gdB_low = gdB;
        break;
    end
    
    Azhigh = Aztry;
    gdB_high = gdB;
    
    Aztry = min(Aztry + dAz,Azmax);
    
end

% Now zoom in...

while (gdB ~= godB && abs(Azhigh - Azlow) >= BWAz*1e-3)
    
    Aztry = (godB*(Azhigh - Azlow) ...
        + Azlow*gdB_high - Azhigh*gdB_low) / (gdB_high - gdB_low);
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Aztry,El,true);    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (gdB > godB + toler)
        Azhigh = Aztry;
        gdB_high = gdB;
    elseif (gdB < godB - toler)
        Azlow = Aztry;
        gdB_low = gdB;
    else
        break;
    end

end

Azright = Aztry;

BWAz = Azright - Azleft;

if LINEAR    % no Elevation search
    return;
end

%--------------------------------------------------------------------------
% Find 3 dB gain points to the left, Elevation.
%--------------------------------------------------------------------------

Elhigh = El;
gdB_high = 0;

BWEl = min(2 / sqrt(H(2,2)),90);

dEl = BWEl / 10;
Eltry = max(El - BWEl/2,Elmin);

while (1)   % find a point < 3 dB gain
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Az,Eltry,true);
    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (Eltry == Elmin)
        gdB = godB;
    end
    
    if (gdB <= godB)
        Ellow = max(Eltry,Elmin);
        gdB_low = gdB;
        break;
    end
    
    Elhigh = Eltry;
    gdB_high = gdB;
    
    Eltry = max(Eltry - dEl,Elmin);
    
end

% Now zoom in...

while (gdB ~= godB && abs(Elhigh - Ellow) >= BWEl*1e-3)
    
    Eltry = (godB*(Elhigh - Ellow) ...
        + Ellow*gdB_high - Elhigh*gdB_low) / (gdB_high - gdB_low);
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Az,Eltry,true);    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (gdB > godB + toler)
        Elhigh = Eltry;
        gdB_high = gdB;
    elseif (gdB < godB - toler)
        Ellow = Eltry;
        gdB_low = gdB;
    else
        break;
    end

end

Elleft = Eltry;

%--------------------------------------------------------------------------
% Find 3 dB gain points to the right
%--------------------------------------------------------------------------

Elhigh = El;
gdB_high = 0;

Eltry = min(max(2*El - Elleft,El + dEl),Elmax);

while (1)   % find a point < 3 dB gain
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Az,Eltry,true);
    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (Eltry == Elmax)
        gdB = godB;
    end
    
    if (gdB <= godB)
        Ellow = Eltry;
        gdB_low = gdB;
        break;
    end
    
    Elhigh = Eltry;
    gdB_high = gdB;
    
    Eltry = min(Eltry + dEl,Elmax);
    
end

% Now zoom in...

while (gdB ~= godB && abs(Elhigh - Ellow) >= BWEl*1e-3)
    
    Eltry = (godB*(Elhigh - Ellow) ...
        + Ellow*gdB_high - Elhigh*gdB_low) / (gdB_high - gdB_low);
    
    Atry = arrayManifold(egain,z,lam,Rarray,Euler,Az,Eltry,true);    
    gdB = 20*log10(abs(Atry'*Ao));
    
    if (gdB > godB + toler)
        Elhigh = Eltry;
        gdB_high = gdB;
    elseif (gdB < godB - toler)
        Ellow = Eltry;
        gdB_low = gdB;
    else
        break;
    end

end

Elright = Eltry;

BWEl = Elright - Elleft;

end