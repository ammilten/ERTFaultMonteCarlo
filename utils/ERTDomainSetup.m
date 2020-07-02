load ../data/electrodes.mat

UR = [max(electrodes.Easting), max(electrodes.Northing)];
LR = [min(electrodes.Easting), min(electrodes.Northing)];
corners_electrode_global = [LR;UR];

electrodes_local = Global2LocalCoords(corners_electrode_global, table2array(electrodes(:,1:2)));

yshift = 200; %200; %meters
xshift = 100; %100; %meters

domain_poly_local = [...
    min(electrodes_local(:,1))-xshift, mean(electrodes_local(:,2))-yshift;... %LR
    max(electrodes_local(:,1))+xshift, mean(electrodes_local(:,2))-yshift;... %UR
    max(electrodes_local(:,1))+xshift, mean(electrodes_local(:,2))+yshift;... %UL
    min(electrodes_local(:,1))-xshift, mean(electrodes_local(:,2))+yshift; ... %LL
    min(electrodes_local(:,1))-xshift, mean(electrodes_local(:,2))-yshift ... %LR
    ];

domain_poly = Local2GlobalCoords(corners_electrode_global,domain_poly_local);
corners = domain_poly(1:2,:);
dims = [max(electrodes_local(:,1))-min(electrodes_local(:,1))+2*xshift, yshift*2];

figure; 
hold on
plot(domain_poly(:,1),domain_poly(:,2))
plot(electrodes.Easting, electrodes.Northing)
scatter(corners(:,1),corners(:,2))
daspect([1 1 1])
legend({'Model Domain','ERT Line','Transformation Corners'},'Location','eastoutside')

save('../data/domain.mat','domain_poly','corners','dims')





