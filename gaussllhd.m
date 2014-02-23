function [ L ] = gaussllhd(Y,X,theta,sigma)

% Theta is a kx1 column vector
% Y is a nx1 column vector
% X is a nxk column vector
% Sigma is a kxk matrix 

L = (2*pi)^(-size(theta,1)/2)*det(sigma)^(-.5)*...
    exp(-.5*(Y - X*theta)'/sigma*(Y - X*theta));

end

