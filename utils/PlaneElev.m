function [a,b,c,d] = PlaneElev(dipdir, dip, outcropX, outcropY, outcropZ)

a = sind(dipdir)*sind(dip);
b = cosd(dipdir)*sind(dip);
c = cosd(dip);
d = sind(dipdir)*sind(dip)*outcropX + cosd(dipdir)*sind(dip)*outcropY + cosd(dip)*outcropZ;




