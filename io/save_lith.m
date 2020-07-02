function save_lith(model,grid,fname)


fid = fopen(fname,'w');

% Size of cube
fprintf(fid,'3\n%d %d %d',size(model,1),size(model,2),size(model,3));
% Spacing in meters
fprintf(fid,'\n%f %f %f',grid.grdx(3),grid.grdy(3),grid.grdz(3));
% Start in meters
fprintf(fid,'\n%f %f %f',grid.grdx(1),grid.grdy(1),grid.grdz(1));

for k = size(model,3):-1:1
    for j = 1:size(model,2)
        
        fprintf(fid,'\n%d',model(1,j,k));
        for i = 2:size(model,1)
            fprintf(fid,' %d',model(i,j,k));
        end
        
    end
end

fclose(fid);
