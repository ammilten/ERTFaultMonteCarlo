grdx = [ 101 200];
grdy = [ 91 500];
grdz = [ 110 1];

GRID = struct('grdx',grdx,'grdy',grdy,'grdz',grdz);

strike = grdx(1):grdx(3):grdx(1) + grdx(3)*(grdx(2)-1);
dip = grdy(1):grdy(3):grdy(1) + grdy(3)*(grdy(2)-1);
