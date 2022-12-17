function [vKappa,mFEX,tStru,SNR,RMS,TRF,Kappa,REJECT,ALGPAR,figh,Ik] = ...
    K_automeasure_dyn2d_sa(A,ta,dt,Nt,tP,tS,Nc,ALGPAR,FIG,mkap)
% Automatically measure Kappa on the given accelerograms 'A', using the 
% method by Anderson & Hough 1984 (fitting a line to freq. vs. log(FAS), 
% where FAS is the Fourier-domain amplitude spectrum of the S-wave signal
% from the acceleration record).
%
% The signal and noise windows and the frequency band for measurement are
% chosen by the algorithm.
%
% This function can take algorithm parameter value input to dynamically
% change the algorithm (optional). Required input are acceleration time
% series data, time vector, sample size, number of samples, onset times
% of P and S waves, number of components to work on, and the optional
% struct 'ALGPAR' which contains algorithm parameter values.
% Last option FIG (optional): plot all traces, spectra in one figure, 
% show 'on', hide 'off', or do not plot anything '' (DEFALT).
%
% Requires Signal-Processing-Toolbox for 'fir1()'.
% 
% OUTPUT:
%   vKappa= (NCOMP,1) accepted Kappa values
%   mFEX  = (NCOMP,2) used f-range for each component
%   tStru = struct('t95',t95,'tNw',tNw,'tSw',tNw); % window times
%   SNR   = (N_BP,NCOMP) signal to noise ratio (f-range,c)
%   RMS   = (N_BP,NCOMP) residual spectrum - kappa_line (f-range,c)
%   TRF   = (N_BP,NCOMP) trade-off values RMS/f_range_length
%   Kappa = (N_BP,NCOMP) one value per f-range and comp
%   REJECT= struct('TRUE',<logical>,'msg',<string>,'rID',<string>)
%   ALGPAR= same as input (or function defaults if not input)
%   figh  = figure handle (empty if no figure)
%   Ik    = (NCOMP,1) indices of accepted Kappa on hor. comp.
%
% 2017-02-16 tsonne: created
% 2017-02-20 tsonne: VERSION 2
%   Added Kappa output.
%   Optimized loop over f-bands some more (precalc. FNS, lFS).
%   Changes in plot: More markers and labels for f_L, f_U.
% 2017-02-22 tsonne: made SSTc and FIRN dynamic variables.
% 2017-07-10 tsonne: VERSION 2b
%   can plot manual kappa input on spectra.
% 2017-07-12 tsonne: wrote output description, added output 'Ik'.
% 2017-07-12 tsonne: VERSION 2c
%   Will plot data even if rejected (produce full overview memo).
% 2017-07-14 tsonne: VERSION 2d
%   'REJECT' struct now contains specific reason ID 'rID' for rejection,
%   which allows later analysis. The possible 'rID's are:
%   ''         = not rejected
%   'tooshort' = (t_end - tS) < TauS2
%   'nonoise1' = could not find noise window start time
%   'nonoise2' = could not find noise window start time
%   'lowsnr'   = have noise window, but all f-bands have too low SNR
% 2019-07-08 tsonne: VERSION 2d_sa
%   standalone script, functions included.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<10, mkap={}; end % mkap is 1x20 cell containing manual values
if nargin<9, FIG=''; end % do not plot anything
% Fill in the algorithm's parameter values:
if nargin<8, ALGPAR=struct(); end
TauS1 = getpar(ALGPAR,'TauS1',0.5); % S-window starts this much BEFORE tS
TauS2 = getpar(ALGPAR,'TauS2',3.0); % min. s-wave window length
TauS3 = getpar(ALGPAR,'TauS3',10.0); % max. s-wave window length
TauN1 = getpar(ALGPAR,'TauN1',0.5); % time before tP to stop noise window
TauN2 = getpar(ALGPAR,'TauN2',3.5); % shortest possible noise window
TauN3 = getpar(ALGPAR,'TauN3',10.0); % max noise window length
TauN4 = getpar(ALGPAR,'TauN4',15.0); % min time after 95% energy for noise
% frequency bands:
BP = getpar(ALGPAR,'BP',[1 12;2 15;4 16;6 18;8 22;(5:10)'*[2 4];...
    [(2:2:14)' 40*ones(7,1)];[(2:2:14)' 45*ones(7,1)]]);
N_BP = getpar(ALGPAR,'N_BP',size(BP,1));
FB_CON = getpar(ALGPAR,'FB_CON',1./sqrt(BP(:,2)-BP(:,1)));
SSTc = getpar(ALGPAR,'SSTc',2.0); % spectral smoothing corner "time"
SNR_LIMIT = getpar(ALGPAR,'SNR_LIMIT',5); % SNR lim for accepting spec
msg_id = getpar(ALGPAR,'msg_id','(ID info N/A)'); % some identifying info 
%                                                   in case of rejection
MAX_LINF = getpar(ALGPAR,'MAX_LINF',50); % linear freq. plot upper limit
YLAB = getpar(ALGPAR,'YLAB',{'X','Y','Z'}); % label for each comp.
% Number of components to work on:
NCOMP = Nc;
% Can be 1, 2, 3. If 3, use geometric mean of hor. comp., not Z-comp.
%
t_end = ta(end);
te1 = t_end*0.002;
% set reasonable value for FIR filter coefficient length, assuming
%  signal T > 3 sec, dt < 0.02 s:
if SSTc > 1.0
    FIRN = 32;
elseif SSTc > 0.3
    FIRN = 64;
elseif SSTc > 0.02
    FIRN = 128;
else
    FIRN = 256;
end
REJECT = struct('TRUE',false,'msg','','rID','');
%
% Start and end times for windows on each component:
t95 = nan(NCOMP,1); % time for energy window (95%)
tNw = nan(NCOMP,2); % time for noise window (spec)
tSw = nan(NCOMP,2); % time for s-wave window (spec)
%
% Kappa measurement for all frequency ranges:
Kappa = zeros(N_BP,NCOMP); % one value per f-range and comp
SNR   = zeros(N_BP,NCOMP); % signal to noise ratio (f-range,c)
RMS   = zeros(N_BP,NCOMP); % residual spectrum - kappa_line (f-range,c)
TRF   = zeros(N_BP,NCOMP); % trade-off values RMS/f_range_length
POLYC = zeros(N_BP,2*NCOMP); % polynomial coefficients
Ik    = zeros(NCOMP,1); % indices of accepted Kappa on hor. comp.
C_SK  = cell(NCOMP,1); % cell for both spectral kappa signatures
C_IS  = cell(NCOMP,1); % indices for both selected f-ranges
% Output for each component:
vKappa = nan(NCOMP,1); % final Kappa value
mFEX   = nan(NCOMP,2); % final frequency limits (f_E, f_X)
tStru  = struct('t95',t95,'tNw',tNw,'tSw',tNw); % window times
%
% check pick and trace quality:
%
% First check: Min. trace length
if (t_end - tS) < TauS2 % if trace too short, skip it
    REJECT.TRUE = true;
    REJECT.msg  = sprintf('Reject: Too short trace: %s',msg_id);
    REJECT.rID  = 'tooshort';
    return
end
% Decide whether to plot and show:
if isempty(FIG)
    FIG = false; % no plots
else
    if islogical(FIG)
        if FIG
            VIS='off'; % plot, but hide figure
        end
    else
        VIS=FIG; % show or hide specified by input
        FIG=true; % do plot
        if ~any(strcmpi(VIS,{'on','yes','y','off','no','n'}))
            % unknown figure option, do not plot anything then
            fprintf('Warning: Unknown figure option: %s\n',VIS);
            FIG=false;
        end
    end
end
figh = []; % figure handle
if FIG
    % Figure, plot parameters:
    FIG_SZ=num2cell([1 1.35]*9); % figure size parameters
    TCOL ='r'; % text color in plots
    T1a = 'Units'; % text option
    T1b = 'normalized'; % text option
    T2a = 'HorizontalAlignment';
    T2b = 'left';
    T2c = 'right';
    T3a = 'FontWeight';
    T3b = 'bold';
    TX  = {T1a,T1b,T2a,T2b,T3a,T3b};
    col = {'color',0.5*[1 1 1]};
    %
    figh(1)=ts_fig(FIG_SZ{:},VIS); % open one figure for all components
end
% VERSION 2b: with manual picks:
if ~isempty(mkap) % check if values given
    mkap(cellfun(@isempty,mkap)) = {nan};
    mk = cell2mat(mkap(1,8:20));
    %   AZ1,FE1,FX1,K1,SNR1,RMS1,AZ2,FE2,FX2,K2,SNR2,RMS2,Azim
    %    8   9  10  11  12   13  14  15  16  17  18   19   20
    %    1   2   3   4   5    6   7   8   9  10  11   12   13
    MAZ = mk([1  7]);
    MFE = mk([2  8]);
    MFX = mk([3  9]);
    MKA = mk([4 10]);
else
    fprintf('%s, both comp: No manual pick info!\n',msg_id);
end
%
for c = 1:NCOMP % go through each comp.
    
    if c == 3 % third  run based on COMP 1 time windows, but will
        % use geometric mean of COMP 1 & 2 spectra later.
        acc = A(:,1);
        maa = max(abs(acc));
    else % COMP 1 & 2
        acc = A(:,c);
        maa = max(abs(acc));
    end
    
    % SGM duration:
    [~,tt95] = eduration(acc,dt);
    t95(c) = tt95(2);
    
    % determine S-wave window:
    % tS       = S onset
    % t95(c)   = t95
    % tSw(c,1) = start
    % tSw(c,2) = end
    % TauS1   = start window this much before tS
    % TauS2   = min length tSw
    % TauS3   = max length tSw
    %
    tSw(c,1) = tS - TauS1;
    if     t95(c) - tS <= TauS2
        tSw(c,2) = tS + TauS2;
    elseif t95(c) - tS >= TauS3
        tSw(c,2) = tS + TauS3;
    else
        tSw(c,2) = t95(c);
    end
    I_SW = (ceil(tSw(c,1)/dt):floor(tSw(c,2)/dt))';
    I_SW(I_SW<1)=[]; % prevent start at neg.time for strange cases
    N_SW = length(I_SW);
    NFFT_S = 2^nextpow2(N_SW);
    tStru.t95(c) = t95(c);
    tStru.tSw(c,:) = tSw(c,:);
    
    % determine Noise window:
    % tP       = P onset
    % t95(c)   = t95
    % tNw(c,1) = start
    % tNw(c,2) = end
    % TauN1   = min time between tNw(c,1) and tP
    % TauN2   = min length tNw
    % TauN3   = max length tNw
    % TauN4   = min time after t95(c) to start tNw(c,1)
    %
    % Noise start time:
    tp1 = (tP - TauN1); % max time from start to signal
    if ~REJECT.TRUE
        [tNw,REJECT] = ...
            K_auto_nst1(c,tp1,tNw,dt,t_end,t95,...
                        TauN2,TauN3,TauN4,REJECT,msg_id);
    end
    % Noise end time:
    if ~REJECT.TRUE
        [tNw,REJECT] = ...
            K_auto_net1(c,tp1,tNw,t_end,t95,TauN2,TauN4,REJECT,msg_id);
        I_NW = (ceil(tNw(c,1)/dt):floor(tNw(c,2)/dt))';
        N_NW = length(I_NW);
        NFFT_N = NFFT_S; % fft interpolation to same freq as S-spectra
        tStru.tNw(c,:) = tNw(c,:);
    end
    
    
    % get FAS of S wave and noise windows:
    %   taper time series window (5%) first.
    % S-wave spectrum:
    if c == 3 % COMP 3 = geomean(COMP1,COMP2)
        tap1 = ts_costaper(A(I_SW,1),0.05,'no');
        tap2 = ts_costaper(A(I_SW,2),0.05,'no');
        [~,F_xx] = ts_fft(tap2,dt,NFFT_S);
        [f,F_S ] = ts_fft(tap1,dt,NFFT_S);
        F_S = sqrt(F_S.*F_xx); % geometric mean S wave spec
    else % COMP 1 & 2
        tap1 = ts_costaper(acc(I_SW),0.05,'no');
        [f,F_S] = ts_fft(tap1,dt,NFFT_S);
    end
    df = f(2)-f(1);
    % Noise spectrum:
    if ~REJECT.TRUE
        if c == 3 % COMP 3 = geomean(COMP1,COMP2)
            tap1 = ts_costaper(A(I_NW,1),0.05,'no');
            tap2 = ts_costaper(A(I_NW,2),0.05,'no');
            [~,F_xx] = ts_fft(tap2,dt,NFFT_N);
            [~,F_N ] = ts_fft(tap1,dt,NFFT_N);
            F_N = sqrt(F_N.*F_xx); % geometric mean noise spec
        else % COMP 1 & 2
            tap2 = ts_costaper(acc(I_NW),0.05,'no');
            [~,F_N] = ts_fft(tap2,dt,NFFT_N);
        end
    end
    
    % smoothen spectra:
    % zero-phase lowpass FIR filtering
%     fs = 1/f(1);
%     fd_h = fdesign.lowpass('Fp,Fst,Ap,Ast',1.6,2.4,1,40,fs);
%     d_d = design(fd_h,'equiripple'); %Lowpass FIR filter
%     F_N = exp(filtfilt(d_d.Numerator,1,log(F_N)));
%     F_S = exp(filtfilt(d_d.Numerator,1,log(F_S)));
    % using fir1():
    coeff = fir1(FIRN,SSTc*2*df);
    F_S = exp(filtfilt(coeff,1,log(F_S))); % smooth S-spectrum
    lFS = log(F_S); % for line fit
    if ~REJECT.TRUE
        F_N = exp(filtfilt(coeff,1,log(F_N))); % smooth Noise spectrum
        %
        FNS = F_N./F_S; % inverse SNR (for harm.mean)
        % FOR EACH FREQ. RANGE, GET KAPPA:
        for bp = 1:N_BP
            % set frequency range: (indices of spectra)
            In = ceil(BP(bp,1)/df):floor(BP(bp,2)/df);
            nf = numel(In); % to calc mean
            % fit line
            %lnF_S = log(F_S(In));
            %pcS = polyfit(f(In),lnF_S,1); % 8-10 times slower
            pcS = fast_polyfit(f(In),lFS(In),1);
            POLYC(bp,2*c-1:2*c) = pcS;
            % get line fit residual:
            %Sval = polyval(pcS,f(In)); % 6 times slower
            Sval = f(In)*pcS(1) + pcS(2);
            Sres = lFS(In) - Sval;
            %RMS(bp,c)  = sqrt(mean(Sres.*Sres)); % direct calc 15% faster
            RMS(bp,c) = sqrt(sum(Sres.*Sres)/nf); %   than using 'mean()'
            % Kappa:
            Kappa(bp,c) = -pcS(1)/pi;
            % Signal-to-noise ratio:
            %snr = F_S(In)./(F_N(In));
            %SNR(bp,c) = harmmean(snr(~isnan(snr))); % (see above)
            SNR(bp,c) = nf/sum(FNS(In));
        end
        %
        % find acceptable f-range for Kappa from RMS,SNR:
        % first set min S/N condition:
        Isnr = find(SNR(:,c) >= SNR_LIMIT);
        if isempty(Isnr)
            REJECT.TRUE = true;
            REJECT.msg  = sprintf('Too low SNR in all f-bands! %s',msg_id);
            REJECT.rID  = 'lowsnr';
        end
        % consider possible rejection, find best option:
        if ~REJECT.TRUE
            % Trade-off RMS vs f-bandwidth: Find optimum
            TRF(:,c) = RMS(:,c) .* FB_CON ;
            [~,Iopt] = min(TRF(Isnr,c));
            Ik(c) = Isnr(Iopt); % idx of accepted Kappa value
            vKappa(c) = Kappa(Ik(c),c); % accepted Kappa
            mFEX(c,:) = BP(Ik(c),:); % used f-range
            %
            % HF decay spectra for accepted Kappa values:
            C_SK{c} = polyval(POLYC(Ik(c),2*c-1:2*c),f);
            % indices of spectra within range:
            C_IS{c} = ceil(BP(Ik(c),1)/df):floor(BP(Ik(c),2)/df);
        end
    end
    
    % IF MANUAL KAPPA VALUES GIVEN, PREPARE FOR PLOT:
    PLOT_MANK = false;
    if c<3 && ~isempty(mkap) % comp 1,2, check for empty cell
        %  ...  also check if f_L is empty
        if ~isnan(MFE(c))
            PLOT_MANK = true;
            MFR = [MFE(c) MFX(c)]; % f-range (manual)
            % set frequency range: (indices of spectral data)
            In = [ceil(MFR(1)/df):floor(MFR(2)/df)];
            nf = numel(In); % to calc mean
            % fit line:
            pcS = fast_polyfit(f(In),lFS(In),1);
            % get line fit residual:
            Sval = f*pcS(1) + pcS(2); % entire freq axis
            Sres = lFS(In) - Sval(In); % only chosen range
            MRMS = sqrt(sum(Sres.*Sres)/nf); % rms line residual
            % Kappa: on current spectrum using given f-band
            MKA2 = -pcS(1)/pi; % can be different from original man.kappa
            % Signal-to-noise ratio:
            if ~REJECT.TRUE
                MSNR = nf/sum(FNS(In));
            end
            % f,A-coordinates for manual f_L, f_U:
            MFLU = [f(In(1)) f(In(end))  ;  exp(Sval([In(1) In(end)]))'];
        else
            fprintf('%s, Comp %d: No manual pick info!\n',msg_id,c);
        end
    end
    
    if FIG %%%%%%%%%% PLOT: %%%%%%%%%%
        if c<3
            ys = -0.35*(c-1);
            % acceleration time series
            subplot('Position',[0.075 0.92+ys 0.91 0.06])
            plot(ta,acc,col{:}); hold on % full acc trace
            if ~REJECT.TRUE
                plot(ta(I_NW),acc(I_NW),'g') % noise window
            end
            plot(ta(I_SW),acc(I_SW),'r') % S wave window
            plot([tP tP],[0 maa],'c',... % P pick indicator
                [tS tS],[0 maa],'c') % S pick indicator
            plot(t95(c)*[1 1],[-maa 0],'b') % 95% energy indic.
            text(tP-te1,0.85*maa,'P',T2a,T2c) % P label
            text(tS-te1,0.85*maa,'S',T2a,T2c) % S label
            text(t95(c)+te1,-0.85*maa,'t_{95}',T2a,T2b) % t95% label
            text(0.02,0.85,['COMP: ' YLAB{c}],TX{[1 2 5 6]}) % COMP ID
            
            axis([0 ta(end) [-maa maa]*1.1])
            ylabel('Acc [cm/s^2]')
            xlabel('Time [s]',TX{[1 2]},'Position',[0.5 -0.3 1])
        end
        if c==3
            ys = -0.62;
        end
        % Plot spectra:
        %
        % Log-Linear spectrum:
        subplot('Position',[0.075 0.68+ys 0.56 0.20])
        if ~REJECT.TRUE
            semilogy(f,F_N,col{:}); hold on % noise spectrum
        end
        semilogy(f,F_S,'k'); hold on        % S wave spectrum
        if ~REJECT.TRUE
            semilogy(f(C_IS{c}),exp(C_SK{c}(C_IS{c})),'r--') % reg.line
            % mark and label lower and upper freq.: (added)
            fA4 = [BP(Ik(c),:);exp(C_SK{c}(C_IS{c}([1 end])))'];
            semilogy(fA4(1,1)*[1 1],fA4(2,1)*[0.5 2],'r') % f_L marker
            semilogy(fA4(1,2)*[1 1],fA4(2,2)*[0.5 2],'r') % f_U marker
            text(fA4(1,1),fA4(2,1)*2.5,'f_L','Color',TCOL)
            text(fA4(1,2)-1,fA4(2,2)*3.0,'f_U','Color',TCOL)
        end
        % mark and label MANUAL f_L,f_U if given:
        if PLOT_MANK
            semilogy(MFLU(1,:),MFLU(2,:),'b--') % reg.line
            semilogy(MFLU(1,1)*[1 1],MFLU(2,1)*[0.5 2],'b') % f_L marker
            semilogy(MFLU(1,2)*[1 1],MFLU(2,2)*[0.5 2],'b') % f_U marker
            text(MFLU(1,1),MFLU(2,1)*2.5,'f_{L,m}','Color','b')
            text(MFLU(1,2)-1,MFLU(2,2)*3.0,'f_{U,m}','Color','b')
            text(0.4,0.86,... % Kappa value, freq-range
                sprintf(['manual \\kappa = %5.3f ',...
                '(%5.3f), f = [%3.1f, %3.1f]'],...
                MKA2,MKA(c),MFR(1),MFR(2)) ,  T1a,T1b,'Color','b')
            text(0.62,0.78,sprintf('man.comp.az.: %3.0f\\circ',MAZ(c)),...
                T1a,T1b,'Color','b')
        end
        % 
        if ~REJECT.TRUE
            text(0.50,0.94,... % Kappa value, freq-range
                sprintf('auto \\kappa = %5.3f, f = [%2.0f, %2.0f]',...
                Kappa(Ik(c),c),BP(Ik(c),1),BP(Ik(c),2)),...
                T1a,T1b,'Color',TCOL)
        end
        text(0.02,0.05,['COMP: ' YLAB{c}],TX{[1 2 5 6]}) % COMP ID
        xlim([0 MAX_LINF])
        ylabel('FAS [cm/s]')
        xlabel('Frequency [Hz]',TX{[1 2]},'Position',[0.5 -0.1 1])
        %
        % Log-Log spectrum:
        subplot('Position',[0.715 0.67+ys 0.265 0.22])
        if ~REJECT.TRUE
            loglog(f,F_N,col{:}); hold on % noise spectrum
        end
        loglog(f,F_S,'k'); hold on        % S wave spectrum
        if ~REJECT.TRUE
            loglog(f,exp(C_SK{c}),'r--')  % decay line
            loglog(fA4(1,1)*[1 1],fA4(2,1)*[0.2 5],'r') % f_L marker
            loglog(fA4(1,2)*[1 1],fA4(2,2)*[0.2 5],'r') % f_U marker
        end
        % mark and label MANUAL f_L,f_U if given:
        if PLOT_MANK
            loglog(f,exp(Sval),'b--') % reg.line
            loglog(MFLU(1,1)*[1 1],MFLU(2,1)*[0.2 5],'b') % f_L marker
            loglog(MFLU(1,2)*[1 1],MFLU(2,2)*[0.2 5],'b') % f_L marker
            % no text labels...
        end
        ylabel('FAS [cm/s]')
        xlabel('Frequency [Hz]',TX{[1 2]},'Position',[0.5 -0.12 1])
        tsSquareLogAxes(gca,F_S,[0.1 100])
    end
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function y = getpar(S,field,default)
%GETPAR get parameter value from a struct
% S         struct
% field     field to be extracted from the S struct
% default   default value if field is not a member of the S struct
%
if isfield(S,field)
  y = getfield(S,field);
elseif nargin>2
  y = default;
else
  error('Need value for struct field: %s',field);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [Tdur,T1T2,I1I2] = eduration(cA,dt,limits)
%
% Calculates the duration of an acceleration time history
% Finds the duration between 5-95% of the energy in the time history
% 
% Optional: give wanted limits as [L1 L2]
%           e.g. limits = [0.02 0.98], to get time from 2% to 98% energy
%
% SYNTAX: [T_90,T1T2] = eduration(cA,dt[, [L1 L2] ]);
% 
% INPUT: cA - Acceleration time history amplitudes
%        dt - sample interval (s)
% 
% OUTPUT: Tdur - duration (s)
%         T1T2 - [t1 t2] relative time confining the duration
%         I1I2 - [I1 I2] indices of duration start and stop on input vector
%
% 2016-01-13 tsonne: created
%
if nargin<3, limits=[0.05 0.95]; end
assert(numel(limits)==2,'Given limits var must be 2el.-vector!');
%
% CALCULATING THE DURATION
cA   = detrend(cA);
cumA = cumsum(cA.^2);
I1   = find(cumA <= limits(1)*cumA(end), 1,'last');
I2   = find(cumA >= limits(2)*cumA(end), 1,'first');
Tdur = (I2-I1)*dt;
I1I2 = [I1 I2];
T1T2 = dt*(I1I2-1); % time starts at zero, not dt
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [tNw,REJECT] = ...
    K_auto_nst1(c,tp1,tNw,dt,t_end,t95,TauN2,TauN3,TauN4,REJECT,msg_id)
% subfunction:
% Noise start time determination
% 20170712tsonne
% 20170714tsonne: added REJECT.rID
% Noise start time:
if     TauN2 <= tp1 && tp1 < TauN3
    tNw(c,1) = dt;
elseif tp1 >= TauN3
    tNw(c,1) = tp1 - TauN3;
elseif tp1 < TauN2 && (t_end-t95(c)-TauN4) >= TauN2 && ...
        (t_end-t95(c)-TauN4) <= TauN3
    tNw(c,1) = t95(c) + TauN4;
elseif tp1 < TauN2 && (t_end-t95(c)-TauN4) > TauN3
    tNw(c,1) = t_end - TauN3;
else
    REJECT.TRUE = true;
    REJECT.msg  = sprintf(...
       'No adequate noise window obtainable! %s',msg_id);
    REJECT.rID  = 'nonoise1';
end
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [tNw,REJECT] = ...
    K_auto_net1(c,tp1,tNw,t_end,t95,TauN2,TauN4,REJECT,msg_id)
% subfunction:
% Noise end time determination
% 20170712tsonne
% 20170714tsonne: added REJECT.rID
% Noise end time:
if     tp1 >= TauN2
    tNw(c,2) = tp1;
elseif tp1 <= TauN2 && (t_end-t95(c)-TauN4) >= TauN2
    tNw(c,2) = t_end;
else
    REJECT.TRUE = true;
    REJECT.msg  = sprintf(...
       'No adequate noise window obtainable! %s',msg_id);
    REJECT.rID  = 'nonoise2';
end
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function y_tap = ts_costaper(y,a,E_yn)
% Cosine tapering of time|data series "y" according to taper ratio "a".
% tsonne 2015-06-06
% 2016-03-07 tsonne:
%     option to avoid automatic energy correction: E_yn = 'no',
%     i.e. do not correct trace energy for loss due to tapering
% 2017-03-02 tsonne:
%     input 'y' can be matrix, each column is treated as signal.
if nargin<3, E_yn=1; end % default: correct energy for taper-loss
TURN = false;
[N,nc] = size(y);
if N==1 % assume single row is one signal
    y=y';
    [N,nc] = size(y);
    TURN = true;
end
t = linspace(0,1,N)';
% cosine taper:
c = [0.5 * (1 - cos(pi*t(t<=a)/a));
     ones(length(t(t>a & t<(1-a))),1);
     0.5*(1-cos(pi*(1-t(t>=(1-a)))/a))] ;% /(1-aa);
if nc>1
    y_tap = bsxfun(@times,y,c);
else
    y_tap = y.*c; % tapered signal
end
% energy matching:
if ischar(E_yn) && any(strcmpi(E_yn,{'no','n'}))
    if TURN, y_tap =y_tap'; end
		return
elseif E_yn==0
    if TURN, y_tap =y_tap'; end
		return
end
E = sum(y.*conj(y),1);
Et = sum(y_tap.*conj(y_tap),1);
% same E as original:
if nc>1
    y_tap = bsxfun(@times,y_tap,sqrt(E/Et));
else
    y_tap = y_tap * sqrt(E/Et);
end
if TURN, y_tap =y_tap'; end
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [f,FAS,FPS] = ts_fft(a,dt,ZERO_PADDING)
% FOURIER ABSOLUTE SPECTRUM (one-sided),  scaled by dt
% optional: phase spectrum.
% if ZERO_PADDING wanted:
%   case 'yes':  automatically use next >= power of 2 as NFFT
%   case number: use  NFFT = ZERO_PADDING
%   else: no zero padding, use input vector length as is
%
% Syntax 
%     [f,FAS,FPS] = ts_fft(a,dt,ZERO_PADDING);
% 
% where f is frequency, FAS is the absolute, and FPS is the phase
% spectra. Input variables are the time history a, timestep dt, and
% the optional zero padding request
% first frequency is df=1/T (After zero-padding)
% tsonne 2015-05-29: created
% tsonne 2015-10-28: better ZERO_PADDING option behaviour
% tsonne 2017-03-02: can input 'a' as matrix now, but use same dt and
%                    output 'f' always as vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TURN = false;
[N,nc] = size(a);
if N==1 % assume single row is one signal
    y=y';
    [N,nc] = size(a);
    TURN = true;
end
if nargin < 3 || isempty(ZERO_PADDING)
    NFFT = N;
else
    if ischar(ZERO_PADDING)
        switch lower(ZERO_PADDING)
            case {'y','yes','pad','next2','nextpow2'}
                NFFT = 2^nextpow2(N);
            case {'n','no'}
                NFFT = N;
            otherwise
                fprintf('# ERROR: Zero padding option value unknown!\n');
                return
        end
    else
        NFFT = ZERO_PADDING;
    end
end
% fft, cut, scale:
Y = fft(a,NFFT);
Y = Y(2:NFFT/2+1,:); % not using the zero frequency value
Y = Y*dt; % scale by sampling interval
% frequency from 1*df to Nyquist frequency:
f = (1:NFFT/2)'/(NFFT*dt);
FAS = abs(Y);
if nargout>2
    FPS = atan2(imag(Y), real(Y));
end
if TURN
    f   = f';
    FAS = FAS';
    if nargout>2
        FPS = FPS';
    end
end
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function p = fast_polyfit(x,y,n)
% FAST_POLYFIT excludes the matrix condition test to save time (tsonne)
% 2017-02-17 tsonne: Mod from polyfit to speed up, also by commenting out
% the warning function calls around 'p = R\(Q'*y);'
%
%POLYFIT Fit polynomial to data.
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 5.17.4.11 $  $Date: 2009/09/03 05:25:22 $
if ~isequal(size(x),size(y))
    error('MATLAB:polyfit:XYSizeMismatch',...
          'X and Y vectors must be the same size.')
end
x = x(:);
y = y(:);
if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end
% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end
% Solve least squares problem.
[Q,R] = qr(V,0);
%ws = warning('off','all'); 
p = R\(Q'*y);    % Same as p = V\y;
%warning(ws);
if size(R,2) > size(R,1)
   warning('MATLAB:polyfit:PolyNotUnique', ...
       'Polynomial is not unique; degree >= number of data points.')
end
if nargout > 1
    r = y - V*p;
    % S is a structure containing three elements: the triangular factor
    % from a QR decomposition of the Vandermonde matrix, the degrees of
    % freedom and the norm of the residuals.
    S.R = R;
    S.df = max(0,length(y) - (n+1));
    S.normr = norm(r);
end
p = p.';          % Polynomial coefficients are row vectors by convention.
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function tsSquareLogAxes(thisax,FAS,XLIM,YMAG)
%
% This function ensures that log-log axes appear square so that functions
% of f^n appear as straight lines with slope n
%
% Syntax: tsSquareLogAxes(thisax,FAS,XLIM,YMAG)
% Input: thisax = axis handle
%        FAS    = Fourier amplitude spectrum
%        XLIM   = (optional) set x-axis limits (Default: 10.^[-1.5 2.5])
%        YMAG   = (optional) y-axis magnitude range, e.g. 4 to show exactly
%                 4 orders of magnitude as in [0.01 100].
%                 OR: YLIM vector, e.g. [0.01 100] to set. 
%                 (Default: same as x-axis, square plot) 
%                 If not same as x-axis, plot area will not
%                 be square, but x- to y-axis unit lengths will be same, 
%                 such that f^n has slope n on plot.
%
% tsonne 2016-09-17
%        2016-11-02 just corrected description above.
%        2017-02-24 update to allow YMAG setting.
%        2018-04-13 allow YLIM vector instead of YMAG
if nargin<3, XLIM = 10.^[-1.5 2.5]; end
m = log10(max(FAS));  % max Y value exp
% upper Y plot limit exp
% if ceil(m/0.5)-m/0.5 < 0.01
%     % if very close to a 0.5 log10 step, go 0.5 up
%     y = ceil(m/0.5)*0.5+0.5;
% else
%     y = ceil(m/0.5)*0.5;
% end
y = ceil(m/0.5)*0.5;
x = log10(XLIM);      % x-axis exp
rx = x(2)-x(1);        % x value range (to be same on y-axis)
if nargin<4
    axis(thisax, [XLIM 10.^[(y-rx) y]])
    set(thisax, 'xtick',10.^(ceil(x(1)):1:floor(x(2))))
    set(thisax, 'ytick',10.^ceil((y-rx-1):1:y))
    axis(thisax, 'square')
else
    if numel(YMAG)==1
        ry = YMAG;
    else
        y = log10(YMAG(2));
        ry = y-log10(YMAG(1));
    end
    axis(thisax, [XLIM 10.^[(y-ry) y]])
    set(thisax, 'xtick',10.^(floor(x(1)):1:ceil(x(2))))
    set(thisax, 'ytick',10.^ceil((y-ry-1):1:y))
    set(thisax, 'PlotBoxAspectRatio', [rx ry ry])
end
end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%