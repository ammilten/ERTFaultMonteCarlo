function [eta, Bv] = ObjectModel_bottom(x,PARAMS)

eta = zeros(2,length(x));
eta(1,:) = (x(end)-x) * PARAMS.St;
eta(2,:) = (x(end)-x) * PARAMS.Sv - PARAMS.theta_a;

Bv = PARAMS.Bv * ones(2,length(x));



