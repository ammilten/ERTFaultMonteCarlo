function Local = Global2LocalCoords(corners,Global)
% This function uses the domain defined by 'corners' and the query points
% in 'Global' to calculate 'Local' coordinates. The origin is at the LR
% corner. Local X coordinates extend from LR -> UR. Local Y coordinates
% extend from LR -> LL
%
% corners: 2x2 matrix
%     Coordinates of the lower right (origin) and upper right (x-axis) 
%     corners of the local domain. Columns are UTM easting and northing. 
%     Row 1 is LR and row 2 is UR.
%         LRe, LRn
%         URe, URn
% Global: nx2 matrix
%     Matrix containing UTM X and Y coordinates of query points
%


%%
% 
% corners = [...
%     332750, 4310950;...
%     331950, 4307000;...
%     ];

theta = atan((corners(2,2)-corners(1,2))/(corners(2,1)-corners(1,1)));

R = [...
    cos(theta), sin(theta);...
    -sin(theta), cos(theta) ...
    ];

Local = (Global-corners(1,:))*R';

