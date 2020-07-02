function [eta, Bv] = ObjectModel(x,PARAMS)

ys = 15000;

eta = zeros(2,length(x));
eta(1,:) = (x(end)-x) * PARAMS.St;
eta(2,:) = x(end)*PARAMS.St - PARAMS.theta_s + PARAMS.Sv*ys - x*PARAMS.Sv;

Bv = PARAMS.Bv * ones(2,length(x));



