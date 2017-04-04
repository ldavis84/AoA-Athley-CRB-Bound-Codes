function [Q,Qcrb,Thresh,GaindB,Pol,iSL] = arrayAthleyPol(egain,...
    SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,Azvals,Elvals,casename)
% [Q,Qcrb,Thresh,GaindB,Pol,iSL] = arrayAthleyPol(egain,...
%    SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,Azvals,Elvals,casename);
% Provides an estimate of the minimum achievable angle-of-arrival (AoA)
% error given that the true signal comes from specified Az and El
% for different SNR values, with polarization sensitivity and accuracy
% of polarization ID also considered.
%    El rotates up from xy plane to z, Az about z right-handed (CCW).
% Source direction is d = [cos(El)*cos(Az); cos(El)*sin(Az); sin(El)];
%    Errors are in terms of local du,dv where uhat is unit direction vector 
% along cross(z,d), where z = [0;0;1], i.e., a local azimuthal rotation.  
% dAz = du / cos(El).  dEl = dv (dv is small rotation up).  So 
% vhat = cross(d,uhat) is local elevation, and [d uhat vhat]
% form a right-handed local coordinate system.
%   The SNR is for the whole array over one snapshot.
%   The signal model is 1 unknown deterministic signal in the presence of
% AWGN of known variance. 
%   The search grid should be sufficiently fine that important sidelobes
% will not be missed.  I.E. dAz or dEl < ~0.1 * lam / D, where D is largest
% array dimension.  The four parameters specifying the the input wave are
% [Az; El; p1; p2] where p1, p2 are the real parameters specifying the
% input polarization (Jones) vector:  
% mu = [cos(p1)*exp(1i*p2); sin(p2)];
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
% Rarray    -- 3 x nelem locations of the array elements.   
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
%              angles in degrees.  
% Elvals    -- nEl x nAz elevation search angles, constant
%              along a row, degrees.  El = pi/2 is straight up along z.
% casename  -- OPTIONAL casename.  If present, figures are plotted and
%              saved to files with this prefix
% Q         -- 4 x 4 x nsnr minimum MSE for each SNR in rad^2:  
%              E{ [du;dEl;dp1,dp2] * [same]' | >= Q(:,:,i) for SNRdB(i).  
%              Note the Jones vector is u = [cos(p1)*exp(1i*p2); sin(p2)];
% Qcrb      -- 4 x 4 x nsnr. Just the Cramer-Rao portion of Q
% Thresh    -- 4 x 2 For each parameter the SNR value and RMS error of 
%              the threshold when the Athley bound is 1 dB to the right of
%              the CRB.  If there Athley is always close to CRB, the values
%              are those for the first SNR value listed, and Inf,Inf if it
%              never gets there.  row 1,2,3,4:  [Az, El, p1, p2].
% GaindB    -- nEl x nAz array gain relative to boresight, dB.  I.E.
%              boresight gain is 0 dB, using unit array manifold vector.
%              For LINEAR array, 1 x nAz only
% Pol       -- nEl x nAz x 2 MLE polarization Jones parameters for each
%              direction
% iSL       -- 1 x nsidelobes indices of detected sidelobes in the Azvals
%              Elvals and GaindB matrices

[dim,nelem] = size(Rarray);
nsnr = length(SNRdB);

if isempty(Euler)
    Euler = zeros(3,nelem);
elseif (size(Euler,1) == 1)
    Euler = [Euler; zeros(2,nelem)];
end

ERROR = dim < 1 || dim > 3;

if ERROR
    fprintf('\narrayAthley:  Rarray dimensions are %d x %d!\n',...
        dim,nelem);
    return;
end

if (dim < 3)    % fill out to three dimensions
    Rarray = [Rarray; zeros(3-dim,nelem)];
end

[nEl,nAz] = size(Azvals);

Az = TrueAzEl(1);
El = TrueAzEl(2);

Thresh = Inf(4,2);

YesPlot = exist('casename','var');

Pol = zeros(nEl,nAz,2);   % optimum polarization choices
G = zeros(nEl,nAz);

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
        xlim(1.1*[xmin xmax]);
    end
    
    if (ymax > ymin)
        ylim(1.1*[ymin ymax]);
    end
    
    axis equal;
end

%--------------------------------------------------------------------------
% Compute the CRB.
%--------------------------------------------------------------------------

Qcrb = CRBAoAPol(egain,SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol);

%--------------------------------------------------------------------------
% Find beamformed array gain over the search grid with polarization.
% Calculate the optimum polarization for each point and the optimum
% inner product.
%--------------------------------------------------------------------------

normalize = true;  

[A1o,A2o] = arrayManifoldPol(egain,lam,Rarray,Euler,Az,El,normalize);

Polo = pol2jones(TruePol);
bo = [A1o A2o]*Polo;
bo = bo / norm(bo);

[A1,A2] = arrayManifoldPol(egain,lam,Rarray,Euler,Azvals,Elvals,normalize);

% Run through all points and get the optimum polarization for each point
% and the optimum gain.

for ia = 1:nAz
    for ie = 1:nEl
        A = [A1(:,ie,ia) A2(:,ie,ia)];
        z = pseudo(A,100*eps)*bo;
        p = jones2pol(z);
        b = A*z;
        b = b / norm(b);
        Pol(ie,ia,:) = p';
        G(ie,ia) = bo'*b;
    end
end

GaindB = db(G);

%--------------------------------------------------------------------------
% Find sidelobes.  They are inner points not close to the boresight that 
% are local maxima of the gain.
%--------------------------------------------------------------------------

dAz = Azvals(1,2) - Azvals(1,1);
dEl = Elvals(2,1) - Elvals(1,1);

notBore = abs(angle(exp(1i*(Azvals - Az))*pi/180))*180/pi > 4*dAz ...
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

iSL = find(iSL);

AzSL = Azvals(iSL);
ElSL = Elvals(iSL);
nSL = length(AzSL(:));

SLLdB = GaindB(iSL);

if (YesPlot)
    subplot(2,2,3);
    imagequick(Azvals(1,:),Elvals(:,1),GaindB,[-20 0]);
    hold on;
    
    if (isempty(iSL))
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
    axis normal;
end

%--------------------------------------------------------------------------
% Now assess Athley bound for each SNR case.
%--------------------------------------------------------------------------

Q = Qcrb;

if (nSL > 0)
    
    errSL = zeros(4,nSL);  % numerical errors for each sidelobe
    
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
    uvSL = ([-v u]'*axisSL) .* repmat(angSL,2,1) * 180/pi;
    uvSL(:,~ii) = repmat([180/2; 180/2],1,sum(~ii));

    errSL(1:2,:) = uvSL;
    
    % put in the polarization errors for each sidelobe
    
    Pol = reshape(Pol,nEl*nAz,2);
    SLP = Pol(iSL,:)';
    errP = 180/pi*angle(exp(1i*pi/180*(SLP - repmat(TruePol',1,nSL))));
    errSL(3:4,:) = errP;
    
    qmax = 180^2 * diag([1/3, 1/6, 1/3, 1/3]);
    
    for k = 1:nsnr
        
        qcrb = Qcrb(:,:,k);
        
        q = athleyunion(SNRdB(k),M,SLLdB,errSL,qcrb);
        
        Q(:,:,k) = maxdef(qcrb,mindef(qmax,q));
    end
end

%--------------------------------------------------------------------------
% Az threshold and plot.
%--------------------------------------------------------------------------

i = 1;
rms_deg = sqrt(squeeze(Q(i,i,:)));
rmscrb_deg = sqrt(squeeze(Qcrb(i,i,:)));

dB1 = 10^(1/20);

plotthresh = false;
iabove = find(rms_deg >= rmscrb_deg * dB1);
if (~isempty(iabove) && iabove(end)+1 <= nsnr)
    ithresh = iabove(end)+1;
    Thresh(i,:) = [SNRdB(ithresh) rms_deg(ithresh)];
    plotthresh = true;
elseif (rms_deg(1) <= rmscrb_deg(1) * dB1)
    Thresh(i,:) = [SNRdB(1) rms_deg(1)];
end

if (YesPlot)
    subplot(2,2,2);
    
    if (plotthresh)
        
        semilogy(SNRdB(:),rms_deg,'b',Thresh(i,1),Thresh(i,2),'bx',...
            SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('Az Error (deg)')
        
        threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
            Thresh(i,1),Thresh(i,2));
        
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

%--------------------------------------------------------------------------
% El plot.
%--------------------------------------------------------------------------

i = 2;
rms_deg = sqrt(squeeze(Q(i,i,:)));
rmscrb_deg = sqrt(squeeze(Qcrb(i,i,:)));

plotthresh = false;
iabove = find(rms_deg >= rmscrb_deg * dB1);
if (~isempty(iabove) && iabove(end)+1 <= nsnr)
    ithresh = iabove(end)+1;
    Thresh(i,:) = [SNRdB(ithresh) rms_deg(ithresh)];
    plotthresh = true;
elseif (rms_deg(1) <= rmscrb_deg(1) * dB1)
    Thresh(i,:) = [SNRdB(1) rms_deg(1)];
end

if (YesPlot)
    
    subplot(2,2,4)
    
    if (plotthresh)
        
        semilogy(SNRdB(:),rms_deg,'b',Thresh(i,1),Thresh(i,2),'bx',...
            SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('El Error (deg)')
        
        threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
            Thresh(i,1),Thresh(i,2));
        
        legend('Athley',threshstr,'CRB');
    else
        
        semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('El Error (deg)')
        
        legend('Athley','CRB');
    end
    
    title(sprintf('Elevation Error vs. SNR, %d Snapshots',M));

    figure(1);
    filename = [casename,'_ArrayAthleyAnalysis'];
    print('-dpng',filename);
end

%--------------------------------------------------------------------------
% Polarization Parameters Plot.
%--------------------------------------------------------------------------

i = 3;
rms_deg = sqrt(squeeze(Q(i,i,:)));
rmscrb_deg = sqrt(squeeze(Qcrb(i,i,:)));

plotthresh = false;
iabove = find(rms_deg >= rmscrb_deg * dB1);
if (~isempty(iabove) && iabove(end)+1 <= nsnr)
    ithresh = iabove(end)+1;
    Thresh(i,:) = [SNRdB(ithresh) rms_deg(ithresh)];
    plotthresh = true;
elseif (rms_deg(1) <= rmscrb_deg(1) * dB1)
    Thresh(i,:) = [SNRdB(1) rms_deg(1)];
end

if (YesPlot)
    figure(2);
    subplot(1,2,1)
    
    if (plotthresh)
        
        semilogy(SNRdB(:),rms_deg,'b',Thresh(i,1),Thresh(i,2),'bx',...
            SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('P1 Error (deg)')
        
        threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
            Thresh(i,1),Thresh(i,2));
        
        legend('Athley',threshstr,'CRB');
    else
        
        semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('Pl Error (deg)')
        
        legend('Athley','CRB');
    end
    
    title(sprintf('Jones Param. 1 Error vs. SNR, %d Snapshots',M));

end

i = 4;
rms_deg = sqrt(squeeze(Q(i,i,:)));
rmscrb_deg = sqrt(squeeze(Qcrb(i,i,:)));

plotthresh = false;
iabove = find(rms_deg >= rmscrb_deg * dB1);
if (~isempty(iabove) && iabove(end)+1 <= nsnr)
    ithresh = iabove(end)+1;
    Thresh(i,:) = [SNRdB(ithresh) rms_deg(ithresh)];
    plotthresh = true;
elseif (rms_deg(1) <= rmscrb_deg(1) * dB1)
    Thresh(i,:) = [SNRdB(1) rms_deg(1)];
end

if (YesPlot)
    
    subplot(1,2,2)
    
    if (plotthresh)
        
        semilogy(SNRdB(:),rms_deg,'b',Thresh(i,1),Thresh(i,2),'bx',...
            SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('P1 Error (deg)')
        
        threshstr = sprintf('Threshold SNR %.1f dB, %.2f Deg RMS',...
            Thresh(i,1),Thresh(i,2));
        
        legend('Athley',threshstr,'CRB');
    else
        
        semilogy(SNRdB(:),rms_deg,'b',SNRdB(:),rmscrb_deg,'r--');
        grid;
        xlabel('Single Snapshot Array SNR (dB)');
        ylabel('P2 Error (deg)')
        
        legend('Athley','CRB');
    end
    
    title(sprintf('Jones Param. 2 Error vs. SNR, %d Snapshots',M));

    figure(2);
    filename = [casename,'_ArrayAthleyAnalysisJones'];
    print('-dpng',filename);
end

end