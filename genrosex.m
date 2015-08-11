%  
%  
% [varargout] = genrosex(index,angles,sigs,varargin) 
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
%  
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% index: list of indices for plotx buttondown functionality, (such as 
%        a channel list), if left empty defaults to [1:length(angles)];
% angles: phase values
% sigs:   scaler (power or snr values to scale radius)
% 
% (optional) 'colors': specify dot colors (optional+1)
% (optional) 'labels': specify text labels (optional+1)
% (optional) 'mpv': plot line of best fit
%  
%  
% OUTPUTS 
% ------------------------------------------------------------------------| 
%  h: various handles
%
%  
% NOTES 
% ------------------------------------------------------------------------| 
% example of usage:
% tmp = rand(1,100).*exp(i*pi/4*(1+randn(1,100)/4));
% inds1 = find(angle(tmp)> pi/4);
% inds2 = find(angle(tmp)<=pi/4);
% colarray = zeros(length(tmp),3);
% colarray(inds1,:) = repmat([0 1 0],length(inds1),1); 
% colarray(inds2,:) = repmat([0 0 1],length(inds2),1); 
% genrosex([],angle(tmp),abs(tmp),'colors',colarray,'mpv')
%
%  
% Written 05/11/2011 
% By Sam Thorpe 


function [varargout] = genrosex(index,angles,sigs,varargin) 


% % parse varargins
if isempty(index), index = [1:length(angles)]; end;
if any(strcmpi(varargin,'fig'))
    hF = varargin{find(strcmpi(varargin,'fig'))+1}; 
else
    hF = figure;
end;
if any(strcmpi(varargin,'ax'))
    ax = varargin{find(strcmpi(varargin,'ax'))+1}; 
else
    ax = gca; 
end;
if any(strcmpi(varargin,'colors'))
    if nargin > 3+find(strcmpi(varargin,'colors'))
        colarray = varargin{find(strcmpi(varargin,'colors'))+1};
    else
        colarray = rand(length(index),3);
    end;
else
    colarray = rand(length(index),3); 
end;
if any(strcmpi(varargin,'unit'))
    ut = varargin{find(strcmpi(varargin,'unit'))+1}; 
else
    ut = []; 
end;
if any(strcmpi(varargin,'yt'))
    yt = varargin{find(strcmpi(varargin,'yt'))+1}; 
else
    yt = 0; 
end;


% % Bounding circles
rad = max(sigs);
allrad = linspace(0,rad,4);
tmp = -pi:pi/100:pi;
for rr = 1:3
    Xt = allrad(rr+1)*cos(tmp);
    Yt = allrad(rr+1)*sin(tmp);
    circ(rr) = plot(ax,Xt,Yt,'linewidth',2,'linestyle',':','color','k');
    ctxt(rr) = text(-allrad(rr+1),0,num2str(round(100*allrad(rr+1))/100));
    set(ctxt(rr),'fontweight','bold');
    hold on;
end;


% % additional lines
lax = line([-rad rad],[0 0],'linewidth',2,'color',[0.3 0.3 0.3],'linestyle',':');
lay = line([0 0],[-rad rad],'linewidth',2,'color',[0.3 0.3 0.3],'linestyle',':');
laxyp = line([-rad rad],[-rad rad],'linewidth',2,'color',[0.3 0.3 0.3],'linestyle',':');
laxyn = line([-rad rad],[ rad -rad],'linewidth',2,'color',[0.3 0.3 0.3],'linestyle',':'); 


% % plot data on top
X = sigs.*cos(angles);
Y = sigs.*sin(angles);
for pin = 1:length(index)
    p(pin) = plot(ax,X(pin),Y(pin),'o');
    set(p(pin),'markersize',10,...
         'markerfacecolor',colarray(pin,:),...
         'markeredgecolor',[0.2 0.2 0.2]); 
         %'markeredgecolor',colarray(pin,:)); 
    hold on;
end;


% % labels?
if any(strcmpi(varargin,'label'))
    if nargin > 3+find(strcmpi(varargin,'label'))
        labels = varargin{find(strcmpi(varargin,'label'))+1}{1};
        thestrings = cell(1,length(labels));
        li = varargin{find(strcmpi(varargin,'label'))+1}{2};
        if isnumeric(labels)
            for Q = 1:length(labels),
                thestrings{Q} = num2str(labels(Q));
            end;
        else
            thestrings = labels;
        end
    else
        labels = index;
        thestrings = cell(1,length(labels));
        li = 1:length(index);
        for Q = 1:length(labels),
            thestrings{Q} = num2str(labels(Q));
        end;
    end;
    txtcolor = [0.2 0.2 0.2];
    h.txtL = text(X(li),Y(li),thestrings);
    set(h.txtL,'fontsize',11,'fontweight','bold','color',txtcolor);
end;


% % mean line?
if any(strcmpi(varargin,'mpv'))
    mpv = circmean(X,Y);
    lmpv = line([cos(mpv) -cos(mpv)],...
                [sin(mpv) -sin(mpv)],...
                'linewidth',2,'color',[1 0 0]);
end


% % set tick marks
set(ax,'fontweight','bold','fontsize',13);
if ~any(strcmpi(varargin,'radians'))
    set(ax,'xtick',[-rad 0 rad],'xticklabel',[225 270 315]);
    set(ax,'ytick',[0 rad],'yticklabel',[180 135]);
else
    %set(ax,'xtick',[-rad 0 rad],'xticklabel',{'-3\pi/4','-\pi/2','-\pi/4'});
    %set(ax,'ytick',[0 rad],'yticklabel',{'\pi,-\pi','3\pi/4'});
    set(ax,'xtick',[-rad 0 rad],'xticklabel',{'-3p/4','-p/2','-p/4'});
    set(ax,'ytick',[0 rad],'yticklabel',{'p,-p','3p/4'});
    set(ax,'fontname','symbol','fontweight','bold','fontsize',14);
end
ax2 = copyobj(ax,gcf);
set(ax2,'xaxislocation','top');
set(ax2,'yaxislocation','right');
if ~any(strcmpi(varargin,'radians'))
    set(ax2,'xtick',[0 rad],'xticklabel',[90 45]);
    set(ax2,'ytick',[0],'yticklabel',[0]);
else
    %set(ax2,'xtick',[0 rad],'xticklabel',{'\pi/2','\pi/4'});
    set(ax2,'xtick',[0 rad],'xticklabel',{'p/2','p/4'});
    set(ax2,'ytick',[0],'yticklabel',[0]);
    set(ax2,'fontname','symbol','fontweight','bold','fontsize',14);
end
set(ax,'xlim',[-rad rad]);
set(ax,'ylim',[-rad rad]);
axis(ax,'square');
set(ax2,'xlim',[-rad rad]);
set(ax2,'ylim',[-rad rad]);
axis(ax2,'square');
yf = 0.05;
tp = get(ax,'position');
set(ax,'position',[tp(1) tp(2)-yt tp(3:end)]);
set(ax2,'position',[tp(1) tp(2)-yt tp(3:end)]);


% % handles
h.fig = hF;
h.ax = ax;
h.plot.p = p;
h.plot.circs = circ;
h.plot.ctxt = ctxt;
h.plot.lax = lax;
h.plot.lay = lay;
h.plot.laxyp = laxyp;
h.plot.laxyn = laxyn;


% plotxxx functionality
axes(ax);
h.txt = text('unit','norm','position',[.05 .9],...
    'string',' ','edgecolor','none','linewidth',2,'fontweight','bold');
set(h.txt,'parent',get(h.plot.p(1),'parent'));
for k = 1:numel(h.plot.p)
    set(h.plot.p(k),'HitTest','on','DisplayName',num2str(index(k)),...
        'ButtonDownFcn',{@lineseriescall,index,k,h.txt,ut});
end    


% % varargouts
if nargout > 0, varargout{1} = h; end


%                             SUBFUNCTIONS 
% % ----------------------------------------------------------------------| 
%  

function [mpv] = circmean(x,y)

if size(x,1)>1, x = x.'; end;
if size(y,1)>1, y = y.'; end;
mps = sum(x.*y)./sum(x.^2);
mpv = atan(mps);


function lineseriescall(src,evt,ind,k,h,ut)

 set(h,'string',[]); pause(0.1);
 set(h,'edgecolor',get(src,'markerfacecolor'),'string',[num2str(ind(k)),ut]);

 
%                               END ALL 
% % ----------------------------------------------------------------------| 
%  





