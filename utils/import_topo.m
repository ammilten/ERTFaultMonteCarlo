F = '/home/ammilten/Documents/LBNL/LiDAR/BareEarth_5m_Arc.asc';

fid = fopen(F);

ln = split(fgetl(fid));
ncols = str2double(ln{2});

ln = split(fgetl(fid));
nrows = str2double(ln{2});
ln = split(fgetl(fid));

LLx = str2double(ln{2});
ln = split(fgetl(fid));
LLy = str2double(ln{2});

ln = split(fgetl(fid));
dx = str2double(ln{2});

ln = split(fgetl(fid));
nodata = str2double(ln{2});

topo = zeros(nrows,ncols);
for i = nrows:-1:1
    ln = split(fgetl(fid));
    lnnum = str2double(ln(2:end));
    lnnum(lnnum==nodata) = nan;
    topo(i,:) = lnnum;
end

fclose(fid);

save('../data/topo.mat','topo','LLx','LLy','dx')

figure;
imagesc(topo)
axis xy

