% This file contains the default settings for Chen
% Defaults are set for figure labels and text properties
% set(groot,'DefaultTextFontName','Times New Roman')
% set(groot,'DefaultAxesFontName','Times New Roman')
% set(groot,'DefaultTextFontName','Helvetica')
% set(groot,'DefaultAxesFontName','Helvetica')
% set(0,'DefaultTextFontWeight','bold')

set(groot,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(groot,'DefaultTextFontWeight','normal')
set(groot,'DefaultTextFontSize',18)
set(groot,'DefaultTextColor','black')

set(groot,'DefaultAxesFontSize',16)
set(groot,'DefaultAxesFontWeight','normal')
set(groot,'DefaultAxesLineWidth',1.0)
set(groot,'DefaultLineLineWidth',1.5)







% % New (perula) colormap
% newcol=[
% 0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% 
% % use 
% % plot([0,1],[0,1],'color',newcol(1,:))
% 
% set(groot,'defaultAxesColorOrder',newcol)

%%% Run the command below to explicitly export graphs as vectors
% exportgraphics(gcf,'figure.pdf','ContentType','vector');