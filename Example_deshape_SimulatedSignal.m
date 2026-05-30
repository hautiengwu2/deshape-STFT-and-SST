clear all ; close all ;


% load necessary toolboxes
% this include some auxillary tools
addpath('./tool') ;


% this is the main folder with the main de-shape code
addpath('./deshape') ;


% change this random seed if you want to see other results
initstate(188882) ;



%% generate simulated signals
Hz = 1000 ;

dd = zeros(1e5+1,1) ;
dd(1:Hz:end) = 1 ;
aa = exp(-([-2e4:8e4]/6e4).^2) ;



% Generate the first intrinsic mode type (IMT) function with nontrivial wave-shape
% function (WSF)
% remove the first harmonic
[h,dh] = hermf(Hz,1,10) ;
h = h ./ max(h) ;

% set special case -- remove the fundamental components
if 0
    xhat =fft(dh) ;
    xhat(1)=0 ; xhat(2)=0 ;
    xhat(end)=0;
    dh = ifft(xhat) ; dh = 10 * dh / norm(dh) ;
    x1 = 3*conv(aa'.*dd, dh,'same') ; % 1Hz
else
    x1 = 3*conv(aa'.*dd, h,'same') ; % 1Hz
end

t = [1:100001]' / Hz ;
IF1 = ones(size(x1)) ;

% Generate the second intrinsic mode type (IMT) function with nontrivial wave-shape
% function (WSF)
% generate instantaneous frequency (IF)
ff = abs(cumsum(randn(size(x1)))) ; IF2 = ff./(max(abs(ff))/2) + pi/2+1.4 ;
IF2 = smooth(IF2, 10000) + t/33 ;
phi = cumsum(IF2) ./ Hz ;

% generate amplitude modulation (AM)
AM2 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 2500) ;
AM2 = 3*AM2 ./ max(AM2) - .5 ;

% this is for wave-shape function (WSF)
gg = mod(phi,1);
[a,b] = findpeaks(gg);
b = [1; b; 2*b(end)-b(end-1)] ;
s2 = zeros(size(phi)) ;
for ii = 1: length(b)-1
    idx = b(ii):b(ii+1) ;
    s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
end
aa = exp(-(abs([-8e4:2e4])/8e4).^(1.1)) ;
x2 = aa' .* s2(1:length(AM2)) ;
x2 = 2.8 * (x2 - mean(x2)) ;


% Generate the third intrinsic mode type (IMT) function with nontrivial wave-shape
% function (WSF)
% generate instantaneous frequency (IF)
ff = abs(cumsum(randn(size(x1)))) ;
IF3 = ff./max(abs(ff)) + 4+0.4 ;
IF3 = smooth(IF3, 10000) ;
IF3(35001:end) = IF3(35001:end) - t(1:65001)/40 ;
phi = cumsum(IF3) ./ Hz ;

% generate amplitude modulation (AM)
AM3 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 6000) ;
AM3 = 3*AM3 ./ max(AM3) - .7 ;

% this is for wave-shape function (WSF)
gg = mod(phi,1);
[a,b] = findpeaks(gg);
s3 = zeros(size(phi)) ;
for ii = 1: length(b)
    idx = max(1,b(ii)-30) : min(b(ii)+30, length(s3)) ;
    s3(idx) = 1 ;
end
x3 = AM3 .* s3 ;
x3(1:35000) = 0 ;
IF3(1:35000) = nan ;


% generate complicated noise
sigma = 0.8 ;
noise = random('T',2,length(x3),1) ;
noise = sigma * noise ;
var(noise)
snrdb = 20 * log10(std(x1+x2+x3)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;


% this is the simulated CLEAN signal
x = x1 + x2 + x3 ;

% this is the simulated NOISY signal
Y = x + noise ;

% downsample to 50Hz.
x = x(1:20:end) ;
Y = Y(1:20:end) ;
t = t(1:20:end) ;
Hz = Hz / 20 ;


% plot the generates simulated signal
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
subplot(411) ;
plot(t,x1(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; 
set(gca,'fontsize', 24) ; axis tight ;
subplot(412) ;
plot(t,x2(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; 
set(gca,'fontsize', 24) ; axis tight ;  
subplot(413) ;
plot(t,x3(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; 
set(gca,'fontsize', 24) ; axis tight ;  
subplot(414) ;
plot(t,x,'k') ; hold on ; set(gca,'fontsize', 18) ; 
set(gca,'fontsize', 24) ; axis tight ; 
xlabel('Time (sec)') ;




figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
subplot(411) ;
plot(t,x1(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([20 80 -inf inf]) ;
set(gca,'fontsize', 24) ; axis([25 65 -inf inf]) ;
subplot(412) ;
plot(t,x2(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([20 80 -inf inf]) ;
set(gca,'fontsize', 24) ; axis([25 65 -inf inf]) ;
subplot(413) ;
plot(t,x3(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([20 80 -inf inf]) ;
set(gca,'fontsize', 24) ; axis([25 65 -inf inf]) ;
subplot(414) ;
plot(t,x,'k') ; hold on ; set(gca,'fontsize', 18) ; axis([20 80 -inf inf]) ;
set(gca,'fontsize', 24) ; axis([25 65 -inf inf]) ;
xlabel('Time (sec)') ;
 



figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*7/10]) ;
subplot(311) ;
plot(t, x, 'k') ; axis([25 65 -3 7]) ;
set(gca,'fontsize', 24) ;
subplot(312) ;
plot(t,noise(1:20:end),'k') ; axis([25 65 -3 7]) ;
set(gca,'fontsize', 24) ;
subplot(313) ;
plot(t,Y,'k') ; axis([25 65 -3 7]) ;
set(gca,'fontsize', 24) ; xlabel('Time (sec)') ;
 




%% setup parameters for de-shape

basicTF.win = 400;

% hop to speed up the calculation
basicTF.hop = 10;

% sampling rate
basicTF.fs = Hz;

% frequency axis resolution
basicTF.fr = 0.02;
basicTF.feat = 'SST11';

advTF.num_tap = 1; 

% window type
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}

% smoothen with a larger alpha
advTF.Smo = 1;

% rejection, can ignore
advTF.Rej = 0;

% threshold
advTF.ths = 1E-9;

% highest/lowest frequency to explore
advTF.HighFreq = 16/Hz ;
advTF.LowFreq = 0.1/Hz ;

% can ignore
advTF.lpc = 0;

% power of soft log 
cepR.g = 0.1;  
cepR.Tc = 0;

% used in fine tuning the algorithm. can ignore
P.num_s = 1;
P.num_c = 1;




%% run de-shape algorithm

% this is clean signal. For noisy signal, replace x by Y
[tfrC, cepsC, ceptic, icepsC, dstfrC, dstfrsqC, tfrsqC, tfrtic] = ...
    deshape(x-mean(x), basicTF, advTF, cepR, P);

% the new sampling periods associated with the hop
DeltaT = basicTF.hop/basicTF.fs;



%% plot all results
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
 
% this is time-varying cepstrum
subplot(2, 2, 1) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), ceptic,...
    abs(cepsC), 0.995); 
axis xy; colormap(1-gray);
xlabel('time (sec)'); ylabel('quef (q)') ; set(gca,'fontsize', 16) ;
title('STCT') ; axis([0 inf 0 2]) ;

subplot(2, 2, 2) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), ceptic, abs(cepsC), 0.995); 
axis xy; colormap(1-gray);
xlabel('time (sec)'); ylabel('quef (q)') ; set(gca,'fontsize', 16) ;
set(gca,'fontsize', 16) ; axis([0 inf 0 2]) ;
hold on; plot(t,1./IF2(1:20:end), 'r', 'linewidth',5) ; 
plot(t,1./IF3(1:20:end), 'b', 'linewidth',5) ; 
plot(t,1./IF1(1:20:end), 'm', 'linewidth',5) ; 

% this is inverse time-varying cepstrum
subplot(2, 2, 3) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, abs(icepsC), 0.995); 
axis xy; colormap(1-gray);
xlabel('time (sec)'); ylabel('freq (Hz)') ; set(gca,'fontsize', 16) ;
title('iSTCT'); axis([0 inf 0 8]) ;

subplot(2, 2, 4) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, abs(icepsC), 0.995); 
axis xy; colormap(1-gray);
xlabel('time (sec)'); ylabel('freq (Hz)') ; set(gca,'fontsize', 16) ;
axis([0 inf 0 8]) ;
hold on; plot(t,IF2(1:20:end), 'r', 'linewidth',5) ; 
plot(t,IF3(1:20:end), 'b', 'linewidth',5) ; 
plot(t,IF1(1:20:end), 'm', 'linewidth',5) ; 





figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;

% this is STFT
subplot(2, 2, 1) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, tfrC, 0.995); 
axis xy; colormap(1-gray);
title('STFT'); xlabel('time (sec)'); ylabel('freq (Hz)') ; 
set(gca,'fontsize', 16) ; axis([0 inf 0 8]) ;

% this is de-shape STFT
subplot(2, 2, 2) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, abs(dstfrC), 0.995); 
axis xy; colormap(1-gray);
xlabel('time (sec)'); ylabel('freq (Hz)') ; set(gca,'fontsize', 16) ;
title('ds-STFT') ; axis([0 inf 0 8]) ;

% this is de-shape SST
subplot(2, 2, 3) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, dstfrsqC, 0.995); 
axis xy; colormap(1-gray);
title('ds-SST'); xlabel('time (sec)'); ylabel('freq (Hz)'); 
axis([0 inf 0 8]) ; set(gca,'fontsize', 16) ;

% this is de-shape SST with ground truth superimposed
subplot(2, 2, 4) ;
imageSQ(0:DeltaT:DeltaT*(size(dstfrsqC,2)-1), tfrtic*basicTF.fs, dstfrsqC, 0.995); 
axis xy; colormap(1-gray);
hold on; plot(t,IF2(1:20:end), 'r', 'linewidth',5) ; 
plot(t,IF3(1:20:end), 'b', 'linewidth',5) ; 
plot(t,IF1(1:20:end), 'b', 'linewidth',5) ; 
title('ds-SST'); xlabel('time (sec)'); ylabel('freq (Hz)'); 
axis([0 inf 0 8]) ; set(gca,'fontsize', 16) ;

