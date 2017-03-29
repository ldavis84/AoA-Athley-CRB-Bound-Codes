function elemGainLoad(Glook,Ellook,Azlook,lamlook,modelook)
% (Glook,Ellook,Azlook,lamlook,modelook);
% Loads variables into the global area for future use by elemGainLookup.
% Supports interpolation of values and wavelengths for either one of two
% polarizations.  After calling this routine, one then uses 
% @elemGainLookup as the function for arrayAthley and other codes.
% Glook        -- nEl x nAz x nelem x nlam x npol stored lookup tables for
%                 antenna elements
% Ellook       -- nEl x 1 reference elevation values, deg
% Azlook       -- nAz x 1 reference azimuth values,deg.  Should include 0
%                 and 360 or -180 and 180 to ensure good circular interp.
% lamlook      -- nlam x 1 reference wavelengths.  If one 1 wavelength or
%                 if Glook is only 2 or 3 dimensions, lamlook is ignored
% modelook     -- use 0 ==> linear interpolation, 1 ==> logarithmic
%                 Logarithmic interpolation works as:
%                 f(x) = f0 .* (f1./f0).^u, where u is a value between 0
%                 and 1.  Linear is f(x) = f0 + (f1-f0).*u.  Log can be 
%                 better for antenna gains because it linearizes phase and
%                 accommodates widely varying gains better.

global Gref AzRef ElRef lamRef modeInterp

% check dimensions

[nEl,nAz,~,nlam] = size(Glook);

Ellook = Ellook(:);

if (length(Ellook) ~= nEl)
    fprintf('elemGainLoad:  lookup table does not match #El values!!\n');
    return;
end

Azlook = Azlook(:);

if (length(Azlook) ~= nAz)
    fprintf('elemGainLoad:  lookup table does not match #Az values!!\n');
    return;
end

if (nlam > 1)
    lamlook = lamlook(:);
    
    if (length(lamlook) ~= nlam)
        fprintf(['elemGainLoad:  lookup table does not ',...
            'match #lam values!!\n']);
        return;
    end
end

Gref = Glook;
AzRef = Azlook;
ElRef = Ellook;
lamRef = lamlook;
modeInterp = modelook;

if (modelook == 1) % ensure no zero values
    
    Ggain = abs(Glook(:));
    toler = max(Ggain)*(10*eps);
    izero = Ggain <= toler;
    
    Gref(izero) = toler;
    
end

end
