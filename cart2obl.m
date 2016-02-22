function [ zeta,lambda, phi ] = cart2obl( x,y,z,c )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    phi = atan(y/x);
    zeta = real(acosh((((x^2 + y^2)^0.5) + z*1i)/c));
    lambda = imag(acosh((((x^2 + y^2)^0.5) + z*1i)/c));
    
end

