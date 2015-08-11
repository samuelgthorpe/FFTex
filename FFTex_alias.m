%  
%  
% [varargout] = FFTex_alias(f,sr,maxf,varargin) 
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
%  
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% f: frequency of sine wave
% sr: sampling rate
%  
% OUTPUTS 
% ------------------------------------------------------------------------| 
%  
%  
% NOTES 
% ------------------------------------------------------------------------| 
%  
%  
% Written 03/07/2013 
% By Sam Thorpe 


function [varargout] = FFTex_alias(f,sr,theta,maxf,varargin) 


if ~nargin  || isempty(f),  f  = 1;  end;
if nargin<2 || isempty(sr), sr = 50; end;
if nargin<3 || isempty(theta), theta = 0; end;
if nargin<4 || isempty(maxf), maxf = 10; end;
sintag = 1;
costag = 0;
if any(strcmpi(varargin,'cos'))
    sintag = 0;
    costag = 1;
end;


fs = 1000;
N = 1;
t = 0:1/fs:N-1/fs;
freqs = 0:1/N:fs-1/N;
tsr = 0:1/sr:N-1/sr;
df = sr/length(tsr);
srfreqs = 0:df:sr-df;
if sintag
    y = sin(2*pi*f*t + theta);
    ysr = sin(2*pi*f*tsr + theta);
elseif costag
    y = cos(2*pi*f*t + theta);
    ysr = cos(2*pi*f*tsr + theta);
end
yF = fft(y)/length(y);
yFsr = fft(ysr)/length(ysr);
yS = abs(yF);
ySsr = abs(yFsr);
ylabstr = 'Amplitude (\muV)';
yT = angle(yF);
yTsr = angle(yFsr);


% % PARSE VARARGOUT
if nargout 
    ts = structure(sr,tsr,ysr);
    varargout{1} = ts;     
end;
if nargout>1  
    fs = structure(df,srfreqs,yFsr);
    varargout{2} = fs;     
end;


%                                 PLOTS 
% % ----------------------------------------------------------------------| 
%  


% % SOME PLOT PARAMETERS
lw = 3;
stxt = 15;
ltxt = 18;


% % PLOT TIME SERIES
hF1 = figure('position',[1725 511 560 420]);
line([0 N],[0 0],'linewidth',lw,'linestyle','--','color',[0.65 0.65 0.65]);
hold on;
plot(t,y,'linewidth',lw,'color',[0 0 1]);
xlabel('Time (sec)','fontweight','bold','fontsize',stxt);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',stxt);
title([num2str(f),' Hz sinusoid: Time Series'],'fontweight','bold','fontsize',ltxt);
grid('on');
set(gca,'fontweight','bold');
set(gca,'xlim',[0 1],'fontweight','bold');
set(gca,'ylim',[-1 1],'fontweight','bold');


% % PLOT SPECTRUM
hF2 = figure('position',[1725 7 560 420]);
bar(freqs,yS,'facecolor','b'); hold on;
xlabel('Frequency (Hz)','fontweight','bold','fontsize',stxt);
ylabel(ylabstr,'fontweight','bold','fontsize',stxt);
title([num2str(f),' Hz sinusoid: Power Spectrum'],'fontweight','bold','fontsize',ltxt);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');
set(gca,'ylim',[0 1],'fontweight','bold');


% % APPEND DOWNSAMPLED PLOTS
if any(strcmpi(varargin,'pause')), pause; end;
figure(hF1);
hp = plot(tsr,ysr,'--go','linewidth',2); 
set(hp,'markersize',8,'markerfacecolor',[0.5 1 0.5],'markeredgecolor',[0.5 0.5 0.5]);
figure(hF2);
hp = plot(srfreqs,ySsr,'go','linewidth',2);
set(hp,'markersize',8,'markerfacecolor',[0.5 1 0.5],'markeredgecolor',[0.5 0.5 0.5]);


% % signal/phase plot
bi = find(srfreqs<=maxf);
clr = repmat([0 1 0],length(bi),1);
genrosex(srfreqs(bi),yTsr(bi),ySsr(bi),'colors',clr,'unit',' Hz','yt',0.045,'radians');
axes('position',[0.15 0.9 0.7 0.1]);
line([0 1],[0 0],'linewidth',2.5,'color',[0 0 0]);hold on;
fnt = [num2str(f),' Hz sinusoid: Phase Angles'];
text(0.2,0.25,fnt,'fontsize',14,'fontweight','bold');
axis('off');


%                              SUBFUNCTIONS 
% % ----------------------------------------------------------------------| 
%  


%                                END ALL 
% % ----------------------------------------------------------------------| 
%  


