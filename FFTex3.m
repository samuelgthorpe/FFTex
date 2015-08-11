%  
%  
% [varargout] = FFTex3(f,N,A,theta,sr,maxf,varargin)
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
% The purpose of this function is to illustrate the basic principles of the
% Fast Fourier Transform (FFT) applied to some example EEG data.
%
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% chan: channel index (can range from 1-65)
% epoch: epoch index (can range from 1-40)
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
% FFtex3 by itself plots data from channel 20, epoch 2. 
% 
%  
% Written 09/11/2012 
% By Sam Thorpe 


function [varargout] = FFTex3(chan,epoch,maxf,varargin)


% % PARSE VARARGIN
if nargin<1 || isempty(chan),  chan = 20; end;
if nargin<2 || isempty(epoch), epoch = 2; end;
if nargin<3 || isempty(maxf),  maxf = 50; end;
load('fakeEEG.mat');


% % COMPUTE TIME SERIES & SPECTRUM
N = 1;
sr = info.sr;
samps = 1:N*sr;
t = info.time(samps);
y = data(samps,chan,epoch);
yF = fft(y)/length(y);
df = sr/length(yF);
freqs = 0:df:sr-df;
yS = abs(yF);
yT = angle(yF);


%                                 PLOTS 
% % ----------------------------------------------------------------------| 
%  


% % signal/phase plot
bi = find(freqs<=maxf);
clr = repmat([0 0 1],length(bi),1);
genrosex(freqs(bi),yT(bi),yS(bi),'colors',clr,'unit',' Hz','yt',0.045,'radians');
%genrosex(freqs(bi),yT(bi),yS(bi),'unit',' Hz','yt',0.045,'radians');
axes('position',[0.15 0.9 0.7 0.1]);
line([0 1],[0 0],'linewidth',2.5,'color',[0 0 0]);hold on;
fnt = ' Example EEG: Phase Angles';
text(0.2,0.25,fnt,'fontsize',14,'fontweight','bold');
axis('off');


% % PLOT SPECTRUM
figure;
bar(freqs,yS,'facecolor','b');
xlabel('Frequency (Hz)','fontweight','bold','fontsize',16);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title(' Example EEG: Amplitude Spectrum','fontweight','bold','fontsize',18);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');


% % PLOT TIME SERIES
figure;
line([0 N],[0 0],'linewidth',3,'linestyle','--','color',[0.65 0.65 0.65]);
hold on;
plot(t,y,'linewidth',3,'color','b'); 
xlabel('Time (sec)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title(' Example EEG: Time Series','fontweight','bold','fontsize',18);
grid('on');
set(gca,'fontweight','bold');
set(gca,'xlim',[0 N],'fontweight','bold');
yL = [-max(abs(y)) max(abs(y))];
set(gca,'ylim',yL,'fontweight','bold');


%                            PARSE VARARGOUT
% % ----------------------------------------------------------------------| 
%  


if nargout 
    ts.sr = sr; ts.t = t; ts.y = y;
    varargout{1} = ts;     
end;
if nargout>1  
    fs.df = df; fs.freqs = freqs; fs.yF = yF;
    varargout{2} = fs;     
end;


%                                END ALL 
% % ----------------------------------------------------------------------| 
%  
