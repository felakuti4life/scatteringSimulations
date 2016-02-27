function [r,xTrans,yTrans] = mapSphereToEllipsoid(zeta, lambda, c,plotit)
    N = 300;
    draw_lambda = linspace(0, 2*pi, N);
    draw_phi = linspace(0,2*pi, N);
    [x,y,z] = obl2cart(zeta,draw_lambda, 0,c);
    if plotit
        plot(x,z);
        hold on
    end
    %lambda = 0.4;
    [inc_x, inc_y, inc_z] = obl2cart(zeta, lambda,0, c);
    if plotit
        scatter(inc_x, inc_z, 'g');
    end
    
    if inc_z > 0
        theta = (atan((c+inc_x)/inc_z) +  atan((inc_x-c)/inc_z))/2;
    else
        theta = (atan((c+inc_x)/inc_z) +  atan((inc_x-c)/inc_z))/2 - pi;
    end
    a = (((inc_x+ c)^2 + (inc_z^2))^0.5 + ((inc_x-c)^2 + inc_z^2)^0.5)/2;
    b = (a^2 - c^2)^0.5;

    %calculate local curvature
    curv = (a*b)/ ((((a^2*(sin(lambda)^2))+(b^2*(cos(lambda))^2))^0.5)^3);
    r = 1/curv;
    
    xTrans = inc_x - (r*sin(theta));
    yTrans = inc_z - (r*cos(theta));
    if plotit
        plotTransposedCircle(r, xTrans, yTrans, N);
        title(sprintf('Y transpose: %f X Transpose: %f', yTrans, xTrans))
        axis equal;
        hold off;
    end
end
                                                              
function [x,y] = plotTransposedCircle(r,xTrans,yTrans,N)
    theta = linspace(0, 2*pi, N);
    x = r * cos(theta)  + xTrans;
    y = r * sin(theta) + yTrans;
    plot(x,y, 'r');
end
