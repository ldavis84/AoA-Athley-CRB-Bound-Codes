function G = elemGainLookup(k,Azvals,Elvals,lam,ipol)
% G = elemGainLookup(k,Azvals,Elvals,lam,ipol); 
% Interpolating gain lookup function that relies on previously calling
% elemGainLoad.m to operate. Values outside the range of the loaded table
% are extrapolated.
% k      -- which element's gain is to be provided
% Azvals -- n x m matrix of azimuth values, deg
% Elvals -- n x m matrix of elevation values, deg, 0 deg is horizontal
%           and 90 is vertical
% lam    -- wavelength
% ipol   -- OPTIONAL, if present, which polarization to look up, 1 or 2
% G      -- n x m matrix of gain values for element k

%==========================================================================
% global variables loaded:
% Gref         -- nEl x nAz x nelem x nlam stored lookup tables for
%                 elements
% ElRef        -- nEl x 1 reference elevation values, rad
% AzRef        -- nAz x 1 reference azimuth values, rad
% lamRef       -- nlam x 1 reference wavelengths, if more than 1
% modeInterp   -- 0 ==> linear interpolation, 1 ==> logarithmic

global Gref AzRef ElRef lamRef modeInterp AzHi

[nEl,~,~,nlam,~] = size(Gref);

if (~exist('ipol','var'))
    ipol = 1;
end

[n,m] = size(Azvals);

Azvals = Azvals(:);
Elvals = Elvals(:);

%--------------------------------------------------------------------------
% Obtain Az/El lookup matrix for this element.  Interpolate lam if
% required.
%--------------------------------------------------------------------------

if (nlam == 1)  % no interpolation
    Go = Gref(:,:,k,1,ipol);
else
    [ilam,xsilam] = interpHelper(lam,lamRef);
    Go = Gref(:,:,k,ilam,ipol);
    G1 = Gref(:,:,k,ilam+1,ipol);
    
    if (modeInterp == 0)   % Linear
        Go = Go*(1-xsilam) + G1*xsilam;
    else
        Go = Go .* (G1 ./ Go).^xsilam;
    end
end

Go = Go(:);   % vectorize

%--------------------------------------------------------------------------
% Perform lookups in az and el directions.  Find corner lookup values for
% each point, then perform interpolation.
%--------------------------------------------------------------------------

kaz = ceil((Azvals - AzHi)/360);     % maps to closest branch cut
Azvals = Azvals - kaz*360;

[iaz,xsi] = interpHelper(Azvals,AzRef);
[iel,nu] = interpHelper(Elvals,ElRef);

ii00 = iel + nEl*(iaz-1);
G00 = Go(ii00);

ii01 = ii00 + 1;    % increment in El
G01 = Go(ii01);

ii10 = ii00 + nEl;      % increment in Az
G10 = Go(ii10);

ii11 = ii00 + (1 + nEl); % both
G11 = Go(ii11);

if (modeInterp == 0)   % Linear
    
    G = G00.*(1-xsi).*(1-nu) + G01.*(1-xsi).*nu ...
        + G10.*xsi.*(1-nu) + G11.*xsi.*nu;
    
else
    G = G00 .* (G01 ./ G00).^((1-xsi).*nu) ...
        .* (G10 ./ G00).^(xsi.*(1-nu)) .* (G11./G00).^(xsi.*nu);
    
end

G = reshape(G,n,m);

end
