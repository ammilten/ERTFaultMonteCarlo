function XYZ = get_topo(domainfile,topofile,varargin)

% domainfile = '../data/domain.mat';
% topofile = '../data/topo.mat';

load(domainfile);
load(topofile);

if ~isempty(varargin)
    dx2 = varargin{1};
else
    dx2 = dx;
end

bbox = polyshape(domain_poly(1:4,1),domain_poly(1:4,2));

x = (LLx:dx:(LLx+dx*(size(topo,2)-1)))';
y = (LLy:dx:(LLy+dx*(size(topo,1)-1)))';

x2 = (LLx:dx2:(LLx+dx*(size(topo,2)-1)))';
y2 = (LLy:dx2:(LLy+dx*(size(topo,1)-1)))';

[X,Y] = meshgrid(x,y);
[X2,Y2] = meshgrid(x2,y2);
[~,I,~] = intersect([X(:),Y(:)],[X2(:),Y2(:)],'rows');

Z = reshape(topo(I),size(X2));
X = X2;
Y = Y2;
x = x2;
y = y2;
clear X2 Y2 x2 y2

inbox = X >= min(domain_poly(:,1)) & ...
    X <= max(domain_poly(:,1)) & ...
    Y >= min(domain_poly(:,2)) & ...
    Y <= max(domain_poly(:,2));

Xc = X(inbox);
Yc = Y(inbox);
Zc = Z(inbox);

in = isinterior(bbox,[Xc,Yc]);

% figure;
% imagesc(x,y,topo)
% axis xy
% hold on
% scatter(Xc(in),Yc(in),'k.')
% hold off

XYZ = [Xc(in), Yc(in), Zc(in)];

