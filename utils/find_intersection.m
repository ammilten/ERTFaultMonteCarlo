function intersections = find_intersection(a,b,c,d,xmin,xmax,ymin,ymax,X,Z)

y = (ymin:ymax)';

% Intersection of north boundary
Zpx = 1/c * (d - a*X - b*ymax);
dfx = Z - Zpx;
try
    xn = interp1(dfx,X,0);
catch ME
    if size(dfx,1) ~= size(unique(dfx),1)
        dfx = dfx + 0.01 * rand(size(dfx));
        xn = interp1(dfx,X,0);
    else
        disp(ME)
    end
end


% Intersection of south boundary
Zpx = 1/c * (d - a*X - b*ymin);
dfx = Z - Zpx;
try
    xs = interp1(dfx,X,0);
catch ME
    if size(dfx,1) ~= size(unique(dfx),1)
        dfx = dfx + 0.01 * rand(size(dfx));
        xs = interp1(dfx,X,0);
    else
        disp(ME)
    end
end


% Intersection of west boundary
Zpy = 1/c * (d - a*xmin - b*y);
dfy = Z(1)*ones(length(y),1) - Zpy;
try
    yw = interp1(dfy,y,0);
catch ME
    if size(dfy,1) ~= size(unique(dfy),1)
        dfy = dfy + 0.01 * rand(size(dfy));
        yw = interp1(dfy,y,0);
    else
        disp(ME)
    end
end

% Intersection of east boundary
Zpy = 1/c * (d - a*xmax - b*y);
dfy = Z(1)*ones(length(y),1) - Zpy;
try
    ye = interp1(dfy,y,0);
catch ME
    if size(dfy,1) ~= size(unique(dfy),1)
        dfy = dfy + 0.01 * rand(size(dfy));
        ye = interp1(dfy,y,0);
    else
        disp(ME)
    end
end

intersections = [...
    xn, ymax, interp1(X,Z,xn);...
    xs, ymin, interp1(X,Z,xs);...
    xmin, yw, Z(1);...
    xmax, ye, Z(end)...
    ];
intersections(any(isnan(intersections), 2), :) = [];

