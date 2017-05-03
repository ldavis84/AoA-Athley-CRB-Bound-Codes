function [Q,Qcrb,AzThresh,ElThresh,GaindB,iSL] = arrayAthley(...
    egain,SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,Azvals,Elvals,casename)
% [Q,Qcrb,AzThresh,ElThresh,GaindB,iSL] = arrayAthley(...
%  egain,SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,Azvals,Elvals,casename);
% Provides an estimate of the minimum achievable angle-of-arrival (AoA)
% error given that the true signal comes from specified Az and El
% with a given polarization for different SNR values.  
%    El rotates up from xy plane to z, Az about z right-handed (CCW). A
% local coordinate system is defined as follows.  Direction TO the source
% is d = [cosd(El) .* cosd(Az); cosd(El) .* sind(Az); sind(El)];
% Next, a local horizontal vector orthogonal to d:
% u = [-sind(Az); cosd(Az); 0];
% Finally, a local vertical v = cross(d,u),
% v = [-sind(El) .* cosd(Az); -sind(El) .* sind(Az); cosd(El)];
%    Polarization is specified as vertical:  E along v, horizontal:  E is 
% along u.  The total polarized wave is along [v u] * p; where p is the
% complex polarization vector given by
% p = [cosd(TruePol(1)) * exp(1i*TruePol(2)*pi/180); sind(TruePol(1)];
% where TruePol is the Jones parameter vector in degrees.
% Small errors in azimuth are changes in estimated d in the u direction.
% This is not the same as azimuth angle because azimuth angle errors get
% large at the zenith, whereas this quantity does not.  Small elevation
% errors are as expected.
%   The SNR is for the whole array over one snapshot.
%   The signal model is 1 unknown deterministic signal in the presence of
% AWGN of known variance. 
%   The search grid should be sufficiently fine that important sidelobes
% will not be missed.  I.E. dAz or dEl < ~0.1 * lam / D, where D is largest
% array dimension.
% egain     -- Handle for an element gain function @g such that
%              G = g(i,Azvals,Elvals,lam,ipol) must return the 
%              amplitude gains of element i in matrix G with same dims
%              as Azvals, Elvals, when these are in degrees.
%              ipol = 1 --> vertical (Z Electric Field), 2 --> Horizontal
%              polarization (+ electric field in cross(z,d) direction 
%              where d is the direction of propagation
% SNRdB     -- 1 x nsnr the SNR values for analysis.
%              SNR is with respect to a the whole array 1 snapshot,
%              i.e., SNRpower = qx * norm(A)^2 / sig^2, where qx is
%              the signal m.s. at a gain-of-1 antenna, norm(A) is the norm
%              of the vector of gains of all apertures, and sig = noise rms
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements.  If 1 x nelem, 
%              then zeroes are filled in y and z.
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about 
%              +z (Azimuth), then about rotated -y (Elev), then about 
%              rotated-rotated +x (Roll). If empty, then all
%              rotations are assumed zero.  If only 1 angle provided, it is
%              assumed to be an azimuth rotation, and the rest are zeros.
% M         -- number of sample points assumed
% TrueAzEl  -- 1 x 2 true source azimuth, elevation in degrees.  
% TruePol   -- 1 x 2 true polarization parameters p1 and p2, deg
% Azvals    -- nEl x nAz azimuth search angles, constant down each column,
%              angles in degrees.  If nEl == 1, then only azimuth 
%              performance is examined, and elevation is fixed.
% Elvals    -- nEl x nAz elevation search angles, constant
%              along a row, degrees.  El = pi/2 is straight up along z.
%              If nEl == 1, then only azimuth 
%              performance is examined, and elevation is fixed at 1st row.
% casename  -- OPTIONAL casename.  If present, figures are plotted and
%              saved to files with this prefix
% Q         -- 2 x 2 x nsnr minimum MSE for each SNR in deg^2:  
%              E{ [du;dEl] * [du;dEl]' | >= Q(:,:,i) for SNRdB(i).  For
%              LINEAR array, both Q and Qcrb are 1 x nsnr Azimuth angle
%              errors only; du = cosd(El) * dAz
% Qcrb      -- 2 x 2 x nsnr. Just the Cramer-Rao portion of Q
% AzThresh  -- 1 x 2 Gain in dB followed by RMS Az error in Deg for the
%              point at which Athley bound is 1 dB higher than CRB.  Empty
%              if no threshold detected in the range of SNR
% ElThresh  -- 1 x 2 Gain in dB followed by RMS El error in Deg for the
%              point at which Athley bound is 1 dB higher than CRB.  Empty
%              if no threshold detected in the range of SNR. 
%              Empty for LINEAR arrays
% GaindB    -- nEl x nAz array gain relative to boresight, dB.  I.E.
%              boresight gain is 0 dB, using unit array manifold vector.
%              For LINEAR array, 1 x nAz only
% iSL       -- 1 x nsidelobes indices of detected sidelobes in the Azvals
%              Elvals and GaindB matrices

pol = pol2jones(TruePol);   % complex polarization vector

AzFirst = Azvals(1,1);
TrueAzEl(1) = TrueAzEl(1) - 360*floor((TrueAzEl(1) - AzFirst)/360);

[dim,nelem] = size(Rarray);
if (dim == 1)
    Rarray = [Rarray; zeros(2,nelem)];
end

nsnr = length(SNRdB);
rad2deg = 180/pi;

if (size(Euler,1) == 1)
    Euler = [Euler; zeros(2,nelem)];
end

ERROR = dim ~= 1 && dim ~= 3;

if ERROR
    fprintf('\narrayAthley:  Rarray dimensions are %d x %d!\n',...
        dim,nelem);
    return;
end

[nEl,nAz] = size(Azvals);

azONLY = (nEl == 1);   % This is the linear array case.

Az = TrueAzEl(1);
El = TrueAzEl(2);

AzThresh = [];    % defaults in case no threshold is found
ElThresh = [];  

YesPlot = exist('casename','var');

%--------------------------------------------------------------------------
% Array geometry.
%--------------------------------------------------------------------------

if YesPlot
    figure(1);
    clf;
    subplot(2,2,1);
    
    plot(Rarray(1,:),Rarray(2,:),'ro');
    xlabel('X Position');
    ylabel('Y Position');
    title(sprintf('Array Plan View, Wavelength = %.3f',lam));
    legend('Phase Centers');
    grid;
    
    xmin = min(Rarray(1,:));
    xmax = max(Rarray(1,:));
    ymin = min(Rarray(2,:));
    ymax = max(Rarray(2,:));
    
    if (xmax > xmin)
        dx = 0.1*(xmax - xmin);
        xlim([xmin-dx xmax+dx]);
    end
    
    if (ymax > ymin)
        dy = 0.1*(ymax - ymin);
        ylim([ymin-dy ymax+dy]);
    end
    
    axis equal;
end

%--------------------------------------------------------------------------
% Compute the CRB.
%--------------------------------------------------------------------------

Qcrb = CRBAoA(egain,SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,azONLY);

%--------------------------------------------------------------------------
% Find beamformed array gain over the search grid.
%--------------------------------------------------------------------------

normalize = true;  % Produce unit vectors

Ao = arrayManifold(egain,pol,lam,Rarray,Euler,Az,El,normalize); 
A = arrayManifold(egain,pol,lam,Rarray,Euler,Azvals,Elvals,normalize);

G = Ao'*reshape(A,nelem,nEl*nAz);
G = reshape(G,nEl,nAz);

GaindB = db(G);

%--------------------------------------------------------------------------
% Find sidelobes.  They are inner points not close to the boresight that 
% are local maxima of the gain.
%--------------------------------------------------------------------------

if azONLY
    iiprev = [1 1:nAz-1];
    iinext = [2:nAz nAz];
    
    dAz = Azvals(1,2) - Azvals(1,1);
    
    notBore = abs(angle(exp(1i*(Azvals - Az)/rad2deg))) > 4*dAz;
    
    iSL = find(GaindB >= GaindB(iiprev) & GaindB > GaindB(iinext) ...
        & notBore);
    
    AzSL = Azvals(iSL);
    nSL = length(iSL);
    
    SLLdB = GaindB(iSL);
    
    if (YesPlot)
        subplot(2,2,3);
        
        if (isempty(iSL))
            plot(Az,0,'bx',Azvals,GaindB,'b-');
            legend('Main Lobe');
        else
            plot(Az,0,'bx',AzSL,GaindB(iSL),'ro',...
                Azvals,GaindB,'b-');
            legend('Main Lobe','Sidelobes');
        end
        
        grid;
        xlabel('Az (Deg)');
        ylabel('Relative Gain (dB)');
        title('Directed Beam Pattern with Sidelobes Identified');
    end

else
    
    dAz = Azvals(1,2) - Azvals(1,1);
    dEl = Elvals(2,1) - Elvals(1,1);
    
    notBore = abs(angle(exp(1i*(Azvals - Az)/rad2deg)))*rad2deg > 4*dAz ...
        | abs(Elvals - El) > 4*dEl;
    
    iel = (1:nEl)';
    ielup = [(2:nEl)'; nEl];
    ieldn = [1; (1:nEl-1)'];
    
    iaz = (1:nAz) - 1;
    iazlf = [1 1:nAz-1] - 1;
    iazrt = [2:nAz nAz] - 1;
    
    iiul = repmat(ielup,1,nAz) + repmat(nEl*iazlf,nEl,1);
    iiu = repmat(ielup,1,nAz) + repmat(nEl*iaz,nEl,1);
    iiur = repmat(ielup,1,nAz) + repmat(nEl*iazrt,nEl,1);
    iil = repmat(iel,1,nAz) + repmat(nEl*iazlf,nEl,1);
    iir = repmat(iel,1,nAz) + repmat(nEl*iazrt,nEl,1);
    iidl = repmat(ieldn,1,nAz) + repmat(nEl*iazlf,nEl,1);
    iid = repmat(ieldn,1,nAz) + repmat(nEl*iaz,nEl,1);
    iidr = repmat(ieldn,1,nAz) + repmat(nEl*iazrt,nEl,1);
    
    iSL = GaindB >= GaindB(iil) & GaindB > GaindB(iir) ...
        & GaindB >= GaindB(iid) & GaindB >= GaindB(iiu) ...
        & GaindB >= GaindB(iidl) & GaindB > GaindB(iidr) ...
        & GaindB >= GaindB(iiul) & GaindB > GaindB(iiur) ...
        & notBore;
    
    AzSL = Azvals(iSL);
    ElSL = Elvals(iSL);
    nSL = length(iSL);
    
    SLLdB = GaindB(iSL);
    
    if (YesPlot)
        subplot(2,2,3);
        imagequick(Azvals(1,:),Elvals(:,1),GaindB,[-20 0]);
        hold on;
        
        if (~any(iSL(:)))
            plot(Az,El,'wx');
            legend('Main Lobe');
        else
            plot(Az,El,'cx',AzSL,ElSL,'bo');
            legend('Main Lobe','Sidelobes');
        end
        
        hold off;
        xlabel('Az (Deg)');
        ylabel('El (Deg)');
        title('Directed Beam Pattern with Sidelobes Identified');
    end
end

%--------------------------------------------------------------------------
% Now assess Athley bound for each SNR case.
%--------------------------------------------------------------------------

Q = Qcrb;

if (nSL > 0)
    
    if azONLY
        uvSL = AzSL - Az;
        qmax = 180^2/3;    % worst it could ever be--know nothing
    else
        % determine direction vectors to sidelobes, measure rotation axis
        % and angle and project onto u an v
        
        ElSL = ElSL(:).';
        AzSL = AzSL(:).';
        
        d = [cosd(El)*cosd(Az); cosd(El)*sind(Az); sind(El)];
        u = [-sind(Az); cosd(Az); 0];
        v = [-sind(El)*cosd(Az); -sind(El)*sind(Az); cosd(El)];
        
        dSL = [cosd(ElSL).*cosd(AzSL); cosd(ElSL).*sind(AzSL); sind(ElSL)];
        angSL = acos(d'*dSL);
        axisSL = [0 -d(3) d(2); d(3) 0 -d(1); -d(2) d(1) 0]*dSL;
        
        ii = abs(exp(1i*angSL) + 1) >= 10*eps; % not opposite poles
        
        axisSL(:,ii) = axisSL(:,ii) ...
            ./ repmat(sqrt(sum(axisSL(:,ii).^2)),3,1);
        uvSL = ([-v u]'*axisSL) .* repmat(angSL,2,1) * rad2deg;
        uvSL(:,~ii) = repmat([180/2; 180/2],1,sum(~ii));

        qmax = 180^2/3 * [1 0; 0 0.25];
    end
    
    for k = 1:nsnr
        
        if azONLY
            qcrb = Qcrb(k);
        else
            qcrb = Qcrb(:,:,k);
        end
        
        q = athleyunion(SNRdB(k),M,SLLdB,uvSL,qcrb);
        
        if azONLY
            Q(k) = max(qcrb,min(qmax,q));
        else
            Q(:,:,k) = maxdef(qcrb,mindef(qmax,q));
        end
    end
end

if (azONLY)
    rms_deg = sqrt(Q(:));
    rmscrb_deg = sqrt(Qcrb(:));
    
    dB1 = 10^(1/20);
    
    iabove = find(rms_deg >= rmscrb_deg * dB1);
    if (~isempty(iabove) && iabove(end)+1 <= nsnr)
        ithresh = iabove(end)+1;
        SNRthresh = SNRdB(ithresh);
        rmsthresh = rms_deg(ithresh);
        
        AzThresh = [SNRthresh, rmsthresh];
    else
        ithresh = [];
        AzThresh = [];
    end
    
    if (YesPlot)   % Just plot the azimuth
        
        subplot(2,2,2);
        
        if (~isempty(ithresh))
            semilogy(SNRdB(:),rms_deg,'b',SNRthresh,rmsthresh,'bx',...
                SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('Az Error (deg)')
            
            threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
                SNRthresh,rmsthresh);
            
            legend('Athley',threshstr,'CRB');
        else
            semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('Az Error (deg)')
            
            legend('Athley','CRB');
        end
        
        title(sprintf('Azimuth Error vs. SNR, %d Snapshots',M));
    end
end

if (~azONLY)
    
    rms_deg = sqrt(squeeze(Q(1,1,:)));
    rmscrb_deg = sqrt(squeeze(Qcrb(1,1,:)));
    
    dB1 = 10^(1/20);
    
    iabove = find(rms_deg >= rmscrb_deg * dB1);
    if (~isempty(iabove) && iabove(end)+1 <= nsnr)
        ithresh1 = iabove(end)+1;
    else
        ithresh1 = [];
    end
    
    if (~isempty(ithresh1))
        
        SNRthresh = SNRdB(ithresh1);
        rmsthresh = rms_deg(ithresh1);
        
        AzThresh = [SNRthresh, rmsthresh];
    end
    
    if (YesPlot)
        
        subplot(2,2,2);
        
        if (~isempty(ithresh1))
            
            semilogy(SNRdB(:),rms_deg,'b',SNRthresh,rmsthresh,'bx',...
                SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('Az Error (deg)')
            
            threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
                AzThresh(1),AzThresh(2));
            
            legend('Athley',threshstr,'CRB');
        else
            
            semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('Az Error (deg)')
            
            legend('Athley','CRB');
        end
        
        title(sprintf('Local Azimuth (du) Error vs. SNR, %d Snapshots',M));
        
    end
    
    rms_deg = sqrt(squeeze(Q(2,2,:)));
    rmscrb_deg = sqrt(squeeze(Qcrb(2,2,:)));
        
    iabove = find(rms_deg >= rmscrb_deg * dB1);
    if (~isempty(iabove) && iabove(end)+1 <= nsnr)
        ithresh2 = iabove(end)+1;
    else
        ithresh2 = [];
    end
    
    if (~isempty(ithresh2))
        
        SNRthresh = SNRdB(ithresh2);
        rmsthresh = rms_deg(ithresh2);
        
        ElThresh = [SNRthresh, rmsthresh];
    end
    
    if (YesPlot)
        
        subplot(2,2,4)
        
        if (~isempty(ithresh2))
            
            semilogy(SNRdB(:),rms_deg,'b',SNRthresh,rmsthresh,'bx',...
                SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('El Error (deg)')
            
            threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
                ElThresh(1),ElThresh(2));
            
            legend('Athley',threshstr,'CRB');
        else
            semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
            grid;
            xlabel('Single Snapshot Array SNR (dB)');
            ylabel('El Error (deg)')
            
            legend('Athley','CRB');
        end
        
        title(sprintf('Elevation Error vs. SNR, %d Snapshots',M));
        
    end

end

if (YesPlot)
    figure(1);
    filename = [casename,'_ArrayAthleyAnalysis'];
    print('-dpng',filename);
end

end