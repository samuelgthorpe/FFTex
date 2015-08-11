%  
%  
% [varargout] = FFTex4(f,N,A,theta,sr,maxf,varargin)
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
% The purpose of this function is to illustrate the basic principles of the
% Discrete Fourier Transform (DFT). 
%
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% maxf: max frequency to display in plot (Hz)
%
%  
% OUTPUTS 
% ------------------------------------------------------------------------| 
% A couple of plots, along with optional outputs:
% varargout{1} ts: time series structure specifying sampling rate, 
%                  time points, and sinusoid amplitude.
% varargout{2} fs: Fourier series structure specifying spectral resolution, 
%                  frequency bin centers, and Fourier coefficients.
%
%  
% NOTES 
% ------------------------------------------------------------------------| 
% Examples of usage:
% FFtex1 by itself plots a 1 Hz sinusoid over a single cycle of 1 second.
% [ts,fs] = FFTex1(1,2,1); returns a 1 Hz sinusoid plotted over 2 cycles
%                          with an amplitude of 2.
% 
%  
% Written 09/11/2012 
% By Sam Thorpe 


function [varargout] = FFTex4(epoch,band,maxf,varargin)


% % PARSE VARARGIN
if nargin<1 || isempty(epoch), epoch = 2; end;
if nargin<2 || isempty(maxf),  maxf = 50; end;
if nargin<3 || isempty(band), band = 11; end;
load('fakeEEG.mat');
hm = 'H65_headmodel_AdultJER.mat';


% % COMPUTE TIME SERIES & SPECTRUM
N = 1;
sr = info.sr;
samps = 1:N*sr;
gce = 1:60;
t = info.time(samps);
y = data(samps,gce,epoch);
yF = fft(y)/length(y);
df = sr/length(yF);
freqs = 0:df:sr-df;
yS = abs(yF);
yT = angle(yF);


% % PARSE VARARGOUT
if nargout 
    ts = structure(sr,t,y);
    varargout{1} = ts;     
end;
if nargout>1  
    fs = structure(df,freqs,yF);
    varargout{2} = fs;     
end;


%                                 PLOTS 
% % ----------------------------------------------------------------------| 
%  


% % PLOT SPECTRUM
figure;
yL = [0 1.1*max(yS(:))];
if length(band)>1
    hpp = patch(band([1 1 2 2]),yL([1 2 2 1]),[0 0 1],'facealpha',0.2);
    hold on;
end;
plotx3(gce,freqs,yS,'linewidth',2);
xlabel('Frequency (Hz)','fontweight','bold','fontsize',16);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title(' Example EEG: Amplitude Spectrum','fontweight','bold','fontsize',18);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');
set(gca,'ylim',yL,'fontweight','bold');


% % PLOT TIME SERIES
figure;
line([0 N],[0 0],'linewidth',3,'linestyle','--','color',[0.65 0.65 0.65]);
hold on;
plotx3(gce,t,y,'linewidth',2); 
xlabel('Time (sec)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title(' Example EEG: Time Series','fontweight','bold','fontsize',18);
grid('on');
set(gca,'fontweight','bold');
set(gca,'xlim',[0 N],'fontweight','bold');
yL = [-max(abs(y(:))) max(abs(y(:)))];
set(gca,'ylim',yL,'fontweight','bold');


% % signal/phase & topo plot
bi = find(freqs>=band(1) & freqs<=band(end));
clr = repmat([0 0 1],length(gce),1);
for q = 1:length(bi)
    
    % % signal/phase 
    genrosex(gce,yT(bi(q),:),yS(bi(q),:),'colors',clr,'yt',0.045,'radians');
    axes('position',[0.15 0.9 0.7 0.1]);
    line([0 1],[0 0],'linewidth',2.5,'color',[0 0 0]);hold on;
    fnt = [num2str(freqs(bi(q))),' Hz Phase Angle by Channel'];
    text(0.15,0.25,fnt,'fontsize',14,'fontweight','bold');
    axis('off');
    
    % % topo plot
    topdat = zeros(1,65); topdat(gce) = yS(bi(q),:);
    topo3_genhm(hm,topdat,gce,[],'label');
    axes('position',[0.15 0.9 0.7 0.1]);
    line([0 1],[0 0],'linewidth',2.5,'color',[0 0 0]);hold on;
    fnt = [num2str(freqs(bi(q))),' Hz Topography'];
    text(0.3,0.25,fnt,'fontsize',14,'fontweight','bold');
    axis('off');
end


%                              SUBFUNCTIONS 
% % ----------------------------------------------------------------------| 
%  


%                                END ALL 
% % ----------------------------------------------------------------------| 
%  
