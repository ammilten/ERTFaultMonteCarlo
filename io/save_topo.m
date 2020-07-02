function save_topo(GRID,Z,toponame)


tid = fopen(toponame,'w');

fprintf(tid,'X,Y,Z');
for i = 1:GRID.grdx(2)
    for j = 1:GRID.grdy(2)
        
        x = GRID.grdx(1) + (i-1) * GRID.grdx(3);
        y = GRID.grdy(1) + (j-1) * GRID.grdy(3);
                        
        fprintf(tid,'\n%f,%f,%f',x,y,max(Z(:,j)));
    end
end

fclose(tid);