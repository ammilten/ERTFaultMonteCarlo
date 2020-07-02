function [m,n,c] = PlaneElevGlobal(dipdir, dip, outcropX, outcropY, topofile)

% Everything EXCEPT Xlocal, Ylocal are global definitions

% Topo file
% domainfile = '../data/domain.mat';
% topofile = '../data/topo.mat';

load(topofile);

x = (LLx:dx:(LLx+dx*(size(topo,2)-1)))';
y = (LLy:dx:(LLy+dx*(size(topo,1)-1)))';
[X,Y] = meshgrid(x,y);

outcropZ = interp2(X,Y,topo,outcropX,outcropY);

m = -sind(dipdir)*tand(dip);
n = -cosd(dipdir)*tand(dip);
c = sind(dipdir)*tand(dip)*outcropX + cosd(dipdir)*tand(dip)*outcropY + outcropZ;




