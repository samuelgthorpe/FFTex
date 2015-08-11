%  
%  
% [varargout] = FFTex1(f,N,A,theta,sr,maxf,varargin)
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
% The purpose of this function is to illustrate the basic principles of the
% Fast Fourier Transform (FFT). The function takes as input the 
% frequency (f) of a sine wave and returns a plot of a sine wave of this 
% frequency. Additional parameters can also be specified (see below). 
%
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% f: frequency
% N: number of seconds
% A: Amplitude of sinusoid
% theta: phase of sinusoid (radians)
% sr: sampling rate (Hz)
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


function [varargout] = FFTex1(f,N,A,theta,sr,maxf,varargin)


% % parse inputs
if ~nargin  || isempty(f), f = 1; end;
if nargin<2 || isempty(N), N = 1; end;
if nargin<3 || isempty(A), A = 1; end;
if nargin<4 || isempty(theta), theta = 0; end;
if nargin<5 || isempty(sr), sr = 1000; end;
if nargin<6 || isempty(maxf), maxf = 10; end;


% % compute time series and spectrum
t = 0:1/sr:N-1/sr;
y = A*sin(2*pi*f*t + theta);
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
btmp = find(abs(freqs-f)==min(abs(freqs-f)));
clr = repmat([0.5 0.5 0.5],length(bi),1);
clr(btmp(1),:) = [0 0 1];
genrosex(freqs(bi),yT(bi),yS(bi),'colors',clr,'unit',' Hz','yt',0.045,'radians');
axes('position',[0.15 0.9 0.7 0.1]);
line([0 1],[0 0],'linewidth',2.5,'color',[0 0 0]);hold on;
fnt = [num2str(f),' Hz sinusoid: Phase Angles'];
text(0.2,0.25,fnt,'fontsize',14,'fontweight','bold');
axis('off');


% % plot spectrum
figure; hold on;
bar(freqs,yS,'facecolor','b');
xlabel('Frequency (Hz)','fontweight','bold','fontsize',16);
ylabel('Amplitude','fontweight','bold','fontsize',15);
title([num2str(f),' Hz sinusoid: Amplitude Spectrum'],'fontweight','bold','fontsize',17);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');
set(gca,'ylim',[0 A],'fontweight','bold');


% % plot time series
figure;
line([0 N],[0 0],'linewidth',3,'linestyle','--','color',[0.65 0.65 0.65]);
hold on;
ytmp = A*sin(2*pi*f*t);
plot(t,ytmp,'linewidth',2,'linestyle','--','color',[0.65 0.65 0.65]); 
plot(t,y,'linewidth',3,'color','b'); 
xlabel('Time (sec)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title([num2str(f),' Hz sinusoid: Time Series'],'fontweight','bold','fontsize',17);
grid('on');
set(gca,'fontweight','bold');
set(gca,'xlim',[0 N],'fontweight','bold');
set(gca,'ylim',[-A A],'fontweight','bold');


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
