% Script to demonstrate usage of kappa estimator algorithm
% 
% Requires 3-component accelerometer data with associated P- and S-phase
% times. Example file 'MyData.mat' is simple struct for one such record.
%
% 2019-07-08 Tim Sonnemann: script created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN SETTINGS
% 
VIS = 'off'; % 'on' show figures, 'off' hide figures
NCOMP = 3; % number of components per station to use (3=geom.mean(1,2))
YLAB  = {'N','E','SQRT(N*E)'}; % Label for component 1,2,3 plots

TauS1 = 0.5;  % S-time window starts by this much BEFORE tS
TauS2 = 3.0;  % min. s-wave window length
TauS3 = 10.0; % max. s-wave window length

TauN1 = 0.5;  % time before tP to stop noise window
TauN2 = 3.5;  % shortest possible noise window
TauN3 = 10.0; % max noise window length
TauN4 = 10.0; % minimum time after 95% energy for noise window at end

SSTc  = 2.0; % spectral smoothing corner "time"
MAX_LINF = 50; % max frequ displayed in semilogy plots (for kappa picks)
MIN_FreqBandWidth = 9; % each tested freq.band has at least this width [Hz]
SNR_LIMIT = 5; % signal to noise limit for accepting Kappa measurement
%
% frequency limits for Kappa estimation:
dfE =  1; dfX = 2; % f increments lower & upper limits
fE  =  1:dfE:15; % lower f limits
fX  = 10:dfX:50; % upper f limits
[FE,FX] = meshgrid(fE,fX);
II = FX-FE >= MIN_FreqBandWidth;
f1 = FE(II); f2 = FX(II);
BP = [f1 f2]; % all used frequency bands
    % trade-off scaling parameter by freq.band:
FB_CON = 1./sqrt(BP(:,2)-BP(:,1));
N_BP = size(BP,1); % number of freq.bands to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STORE SETTINGS FOR FUNCTIONS
% Fill in the algorithm's parameter values:
ALGPAR = struct();
ALGPAR.TauS1 = TauS1; % S-window starts this much BEFORE tS
ALGPAR.TauS2 = TauS2; % min. s-wave window length
ALGPAR.TauS3 = TauS3; % max. s-wave window length
ALGPAR.TauN1 = TauN1; % time before tP to stop noise window
ALGPAR.TauN2 = TauN2; % shortest possible noise window
ALGPAR.TauN3 = TauN3; % max noise window length
ALGPAR.TauN4 = TauN4; % min time after 95% energy for noise
ALGPAR.BP = BP; % frequency band limits
ALGPAR.N_BP = N_BP; % length(BP)
ALGPAR.FB_CON = FB_CON; % freq.band trade-off const
ALGPAR.SSTc = SSTc; % spectral smoothing corner "time"
ALGPAR.SNR_LIMIT = SNR_LIMIT; % SNR lim for accepting spec
ALGPAR.YLAB = YLAB; % ylabels for each component (cell of strings)
ALGPAR.MAX_LINF = MAX_LINF; % linear freq. plot upper limit
ALGPAR.msg_id = ''; % some identifying info (see below in loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DATA AND RUN KAPPA ESTIMATOR
load('MyData.mat');
%
% Get 3comp. accelerograms, timevector, sample length, P and S phase times:
A = MyData.A;
ta = MyData.ta;
dt = MyData.dt;
Nt = MyData.Nt;
tP = MyData.tP;
tS = MyData.tS;
A = detrend(A);
%
% Include manual kappa value info in plot (if plotting):
% mkap = row cell, 20 elements, using 8:20
mkap = {}; % don't add other kappa measurements to plot
%
% Kappa measurement algorithm:
[vKappa,mFLU,tStru,SNR,RMS,TRF,KAP,REJECT,~,~,Ik] = ...
    K_automeasure_dyn2d_sa(A,ta,dt,Nt,tP,tS,NCOMP,ALGPAR,VIS,mkap);
if ~REJECT.TRUE % if traces not rejected
    kr = mean(vKappa(1:2)); % average horizontal Kappa_r
end