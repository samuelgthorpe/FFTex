%  
%  
% [varargout] = FFTex2(varargin) 
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
% The purpose of this function is to illustrate the basic principles of the
% Fast Fourier Transform (FFT). The function takes as input the 
% number of harmonics desired to approximate a square wave. 
%
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% Nharms: number of harmonics
% 
%  
% OUTPUTS 
% ------------------------------------------------------------------------| 
%  
%  
% NOTES 
% ------------------------------------------------------------------------| 
% FFTex2 by itself only plots the first harmonic approximation (i.e. just a
% sinusoid). For the output of this funtion to be interesting try
% successive approximations:
% for n = 1:5, FFTex2(n) end;
%
%  
% Written 03/08/2013 
% By Sam Thorpe 


function [varargout] = FFTex2(Nharms,varargin) 


% % general parameters
if ~nargin  || isempty(Nharms), Nharms = 1;  end;
f = 1:2:2*Nharms;
N = 1;   
plotphase = 0;
theta = 0; 
maxf = max(10,2*max(f)+1);
A = 1;
fs = 1000;
t = 0:1/fs:N-1/fs;
yA = zeros(length(t),Nharms);
hclr = zeros(Nharms,3);


% % make sum of sinusoids
for Q = 1:Nharms
    yA(:,Q) = 2*(2*A)/(pi*f(Q))*sin(2*pi*f(Q)*t - theta);
    hclr(Q,:) = (Q-1)/Nharms*[0.75 0.75 0.75];
end
y = sum(yA,2);


% % make square wave for comparison
ySQ = zeros(length(t),1);
for Q = 1:N
   ySQ((Q-1)*fs + (1:fs/2)) = 1;
   ySQ((Q-1)*fs + (fs/2+1:fs)) = -1;
end


% % compute spectra
yF = fft(y)/length(y);
ySQF = fft(ySQ)/length(ySQ);
df = fs/length(yF);
freqs = 0:df:fs-df;
yS = 2*abs(yF);
ySQS = 2*abs(ySQF);
yT = angle(yF);


%                                 PLOTS 
% % ----------------------------------------------------------------------| 
%  


% % PLOT ANGLE
if plotphase
    figure('position',[500 511 560 420]);
    bar(freqs,yT.*yS,'facecolor','r');
    xlabel('Frequency (Hz)','fontweight','bold','fontsize',15);
    ylabel('Angle*Amplitude (radians*\muV)','fontweight','bold','fontsize',15);
    title([num2str(f),' Hz sinusoid: Phase Angle'],'fontweight','bold','fontsize',18);
    grid('on');
    set(gca,'xlim',[0 maxf],'fontweight','bold');
    set(gca,'ylim',[-pi pi],'fontweight','bold');
    ytl = {'-pi','-pi/2','0','pi/2','pi'};
    set(gca,'ytick',[-pi -pi/2 0 pi/2 pi],'yticklabel',ytl);
end

  
% % PLOT SPECTRUM
figure('position',[500 511 560 420]);
bar(freqs,yS,'facecolor','b');
xlabel('Frequency (Hz)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title('Sum of odd sine series: Amplitude Spectrum','fontweight','bold','fontsize',18);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');
set(gca,'ylim',1.5*[0 A],'fontweight','bold');


% % PLOT SPECTRUM
figure('position',[500 511 560 420]);
bar(freqs,ySQS,'facecolor','r');
xlabel('Frequency (Hz)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title('Square Wave: Amplitude Spectrum','fontweight','bold','fontsize',18);
grid('on');
set(gca,'xlim',[0 maxf],'fontweight','bold');
set(gca,'ylim',1.5*[0 A],'fontweight','bold');


% % PLOT TIME SERIES
figure('position',[500 511 560 420]);
line([0 N],[0 0],'linewidth',3,'linestyle','--','color',[0.65 0.65 0.65]);
hold on;
for Q = 1:Nharms
    plot(t,yA(:,Q),'linewidth',1,'linestyle','-','color',hclr(Q,:));
end
plot(t,ySQ,'linewidth',2,'color',[1 0 0]); 
plot(t,y,'linewidth',3,'color',[0 0 1]); 
xlabel('Time (sec)','fontweight','bold','fontsize',15);
ylabel('Amplitude (\muV)','fontweight','bold','fontsize',15);
title('Sum of odd sine series: Time Series','fontweight','bold','fontsize',18);
grid('on');
set(gca,'fontweight','bold');
set(gca,'xlim',[0 N],'fontweight','bold');
set(gca,'ylim',1.5*[-A A],'fontweight','bold');


%                            PARSE VARARGOUT
% % ----------------------------------------------------------------------| 
%  

if nargout,   varargout{1} = t;     end;
if nargout>1, varargout{2} = y;     end;
if nargout>2, varargout{3} = freqs; end;
if nargout>3, varargout{4} = yF;    end;

%                                END ALL 
% % ----------------------------------------------------------------------| 
%  


