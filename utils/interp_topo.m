function z = interp_topo(topofile,coords)

load(topofile);

x = (LLx:dx:(LLx+dx*(size(topo,2)-1)))';
y = (LLy:dx:(LLy+dx*(size(topo,1)-1)))';
[X,Y] = meshgrid(x,y);

if size(coords,1) == 1
    z = interp2(X,Y,topo,coords(1),coords(2));
else
    z = interp2(X,Y,topo,coords(:,1),coords(:,2));
end