function [m,n,c] = PlaneElevLocal(dipdir, dip, outcropX, outcropY, topofile, domainfile)

% Everything EXCEPT Xlocal, Ylocal are global definitions

% Topo file
% domainfile = '../data/domain.mat';
% topofile = '../data/topo.mat';

load(domainfile);
% load(topofile);
% 
% x = (LLx:dx:(LLx+dx*(size(topo,2)-1)))';
% y = (LLy:dx:(LLy+dx*(size(topo,1)-1)))';
% [X,Y] = meshgrid(x,y);
% 
% outcropZ = interp2(X,Y,topo,outcropX,outcropY);
% 
% mg = -sind(dipdir)*tand(dip);
% ng = -cosd(dipdir)*tand(dip);
% cg = sind(dipdir)*tand(dip)*outcropX + cosd(dipdir)*tand(dip)*outcropY + outcropZ;

[mg, ng, cg] = PlaneElevGlobal(dipdir,dip,outcropX,outcropY,topofile);

theta = atan((corners(2,2)-corners(1,2))/(corners(2,1)-corners(1,1)));
m = ng*sin(theta) + mg*cos(theta);
n = ng*cos(theta) + mg*sin(theta);
c = mg*corners(1,1) + ng*corners(1,2) + cg;




