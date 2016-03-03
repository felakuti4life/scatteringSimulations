function [ x,y,z ] = obl2cart( zeta,lambda,phi, c )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
        x = c * cosh(zeta) * cos(lambda) * cos(phi); 

    y = c * cosh(zeta) * cos(lambda) * sin(phi);
    z = c * sinh(zeta) * sin(lambda);
end

