function [a,x,y] = mapSphereToEllipsoid(zeta, lambda, c)
    N = 300;
    draw_lambda = linspace(0, 2*pi, N);
    draw_phi = linspace(0,2*pi, N);
    [x,y,z] = obl2cart(zeta,draw_lambda, 0,c);
    plot(x,z);
    hold on
    %lambda = 0.4;
    [inc_x, inc_y, inc_z] = obl2cart(zeta, lambda,0, c);
    scatter(inc_x, inc_z, 'g');
    theta = (atan((c+inc_x)/inc_z) +  atan((c-inc_x)/inc_z))/2
    r = zeta*0.2;
    xTrans = inc_x - (r*sin(theta));
    yTrans = inc_z - (r*cos(theta));
    plotTransposedCircle(r, xTrans, yTrans, N);
    title(sprintf('Y transpose: %f X Transpose: %f', yTrans, xTrans))
    axis equal;
    hold off;
end
                                                              
function [x,y] = plotTransposedCircle(r,xTrans,yTrans,N)
    theta = linspace(0, 2*pi, N);
    x = r * cos(theta)  + xTrans;
    y = r * sin(theta) + yTrans;
    plot(x,y, 'r');
end

function [ zeta,lambda, phi ] = cart2obl( x,y,z,c )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    phi = atan(y/x);
    zeta = real(acosh((((x^2 + y^2)^0.5) + z*1i)/c));
    lambda = imag(acosh((((x^2 + y^2)^0.5) + z*1i)/c));
    
end

function [ x,y,z ] = obl2cart( zeta,lambda,phi, c )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    x = c * cosh(zeta) * cos(lambda) .* cos(phi);
    y = c * cosh(zeta) * cos(lambda) .* sin(phi);
    z = c * sinh(zeta) * sin(lambda);

end