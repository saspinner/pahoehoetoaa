%%%%% plot Regime Diagram%%%%%%
function JplotRD(N,FF,IJ,fs)
%input: N is the height ratio; FF is the Froude number; IJ is the result;  
% fs: is free-slip or not
%output: Regime Diagram



% set printing options
if fs == 1
str       = {'Regime_FreeSurface'};
elseif fs == 0
    str   = {'Regime_ConstSlip'};
elseif fs == 2
    str  = {'Regime_NoSlip'};
end
figname   = char(strcat(str(1)));
format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';
printfig  = true;


% prepare formating options
HA = {'HorizontalAlignment','left','center','right'};
VA = {'VerticalAlignment','bottom','middle','top'};
UN = {'Units','Normalized','Inches'};
TX = {'Interpreter','Latex'};
TL = {'TickLabelInterpreter','Latex'};
LW = {'LineWidth',1,1.25,1.5,2};
FS = {'FontSize',10,15,18,21,24};
MS = {'MarkerSize',6,8,12};
LS = {'LineStyle','-','--','-.',':'};
% LC = {'Color',color};

% prepare axes/borders dimensions
axh = 3*2;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.2;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2;
fh = axb + 1*axh + 1*avs +           axt;
fw = axl + 1*axw + 1*ahs + 1.5*cbw + axr;

% initialize figure and axes
f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax(1) = axes('Units','Inches','position',[axl         axb         axw axh]);


axes(ax(1)); hold on;
contourf(N,FF,IJ,'edgecolor','none');
colormap([ 33,102,172;253,219,199;239,138,98;178,24,43]./255)
hold on


set(gca,'YScale', 'log')
xlabel('Height ratio, n',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Froude number, F',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);

print(f,format,resl,rend,figname,'-loose');