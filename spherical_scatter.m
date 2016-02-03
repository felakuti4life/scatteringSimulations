function scatter_pressure = spherical_scatter( incident, k, a)
%SPHERICAL_SCATTER scattering function for a plane wave hitting a sphere
%   Cylinder is assumed to be located at 0.
%   Input args:
%   incident: [x y z] vector holding cartesian coordinates of incident wave
%   a: radius of shpere
    r = (incident(1)^2 + incident(2)^2)^0.5;
    phi = atan(incident(2)/incident(1));
    ka = k*a;
    n = ceil(2*(ka+1));
    M = [1:n-1].';

    J = [-besselj(1,ka); besselj(M-1,ka)-besselj(M+1,ka)];
    Y = [ bessely(1,ka); bessely(M+1,ka)-bessely(M-1,ka)];

    Gamma = atan2(J,Y);
    epsilon = [1; 2*ones(n-1,1)];
    cosOfe = cos([0:n-1].'*phi);
    bess = besselj([0:n-1].',k*r)+1i*bessely([0:n-1].',k*r);
    ith = 1i.^([0:n-1].'+1);
    scatter_pressure = -sum(epsilon.*cos(Gamma).*exp(-1i*Gamma).*cosOfe.*bess.*ith);
    
end