function scatter_pressure = elliptical_scatter( incident, k, a, c )
%ELLIPTICAL_SCATTER Summary of this function goes here
%   Scattering of an ellipse
%   incident: [x y z] vector holding cartesian coordiantes of incident wave
%   k: wave number
%   a: zeta of the ellipse (hyperbolic radius)
%   c: foci point
    [zeta, lambda, phi] = cart2obl(incident(1), incident(3), incident(2), c);
    [r, cir_x, cir_y] = mapSphereToEllipsoid(a, lambda,c, 0);
    if(incident(1) < 0)
        cir_x = -cir_x;
    end
    scatter_pressure = spherical_scatter([incident(1) + cir_x, incident(2) + cir_y, incident(3)], k, r);
    
    LAM = linspace(0,2*pi, 40);
    theta = linspace(0, 2*pi, 40);
    x = r * cos(theta)  + cir_x;
    y = r * sin(theta) + cir_y;
    plot(x,y, 'r');
    hold on
    [san_x,san_y,san_z] = obl2cart(zeta,lambda,phi,c);
    scatter(san_x,san_z,'k');
    %scatter(incident(1),incident(2));
    title(sprintf('zeta: %0.2f lambda: %0.2f', zeta, lambda));
    [hitp_x, hitp_y, hitp_z] = obl2cart(a, lambda, c, 0);
    plot([incident(1),hitp_x], [incident(2), hitp_z], 'k');
    [el_x, el_y, el_z] = obl2cart(a, LAM, 0,c);
    plot(el_x,el_z, 'm--');
    axis equal
    pause(0.1);
    %hold off
end

