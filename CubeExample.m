%--------------------------------------------------------------------------
% Try out an array on a cube to demonstrate its properties.
%--------------------------------------------------------------------------

clear;
rng(1);

descriptor = 'CubeExample';

lam = 2;  % wavelength
TrueAzEl = [30,45];
TruePol = [0; 0];    % vertical pol

% Cube Geometry.  Build one face then replicate by displacement and 
% rotation.

L = 5;
npface = 5;   % elements per face

Rproto = L* [rand(2,npface)-0.5; zeros(1,npface)];
Eulproto = [360*rand(1,npface); zeros(2,npface)];

% top--just raise it up

RTop = Rproto + [0;0;L]*ones(1,npface);
EulTop = Eulproto;

% +X face

QgX = euler2Q(0,-90,0);
RXp = QgX*RTop + L/2*[-1;0;1]*ones(1,npface);

EulXp = 0*EulTop;
for i = 1:npface
    Qpi = euler2Q(Eulproto(1,i),Eulproto(2,i),Eulproto(3,i));
    [az,el,roll] = Q2euler(QgX * Qpi);
    EulXp(:,i) = [az; el; roll];
end

% -X face

QgX = euler2Q(0,90,0);
RXm = QgX*RTop + L/2*[1;0;1]*ones(1,npface);

EulXm = 0*EulTop;
for i = 1:npface
    Qpi = euler2Q(Eulproto(1,i),Eulproto(2,i),Eulproto(3,i));
    [az,el,roll] = Q2euler(QgX * Qpi);
    EulXm(:,i) = [az; el; roll];
end

% +Y face

QgY = euler2Q(0,0,-90);
RYp = QgY*RTop + L/2*[0;-1;1]*ones(1,npface);

EulYp = 0*EulTop;
for i = 1:npface
    Qpi = euler2Q(Eulproto(1,i),Eulproto(2,i),Eulproto(3,i));
    [az,el,roll] = Q2euler(QgY * Qpi);
    EulYp(:,i) = [az; el; roll];
end

% -Y face

QgY = euler2Q(0,0,90);
RYm = QgY*RTop + L/2*[0;1;1]*ones(1,npface);

EulYm = 0*EulTop;
for i = 1:npface
    Qpi = euler2Q(Eulproto(1,i),Eulproto(2,i),Eulproto(3,i));
    [az,el,roll] = Q2euler(QgY * Qpi);
    EulYm(:,i) = [az; el; roll];
end

% Assembled Cube

Rarray = [RTop RXp RYp RXm RYm];
Euler = [EulTop EulXp EulYp EulXm EulYm];

% Spacing in azimuth and elevation of the beampattern sample points in
% degrees
daz = 1;
del = 1;

% Computes the azimuth and elevation values to be evaluated
azvals = (-180 : daz : 180 + daz);
nAz = length(azvals);
elvals = (0: del : 90 + del)';
nEl = length(elvals);

% This is the format in which the azimuth and elevation values are passed 
% to the function
Azvals = repmat(azvals,nEl,1);
Elvals = repmat(elvals,1,nAz);

% Number of snapshots averaged together (Pick a number between 5 and 10)
M = 5;

% Signal to noise ratios (including the beamforming gain) to be evaluated
SNRdB = 0:0.1:30;

% Compute the arrayAthley bounds
[Q,Qcrb,AzThresh,ElThresh,GaindB,iSL] = arrayAthley(@horizB, ...
    SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol,Azvals,Elvals,descriptor);
