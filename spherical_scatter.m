function scatter_pressure = spherical_scatter( incident, k, a)
%SPHERICAL_SCATTER scattering function for a plane wave hitting a sphere
%   Sphere is assumed to be located at 0.
%   Input args:
%   incident: [x y z] vector holding cartesian coordinates of incident wave
%   a: radius of shpere
    r = (incident(1)^2 + incident(2)^2)^0.5;
    phi = atan(incident(2)/incident(1));
    theta = acos(incident(3)/r);
    
    
    ka = k*a;
    n = ceil(2*(ka+1));
    M = [1:n-1].';
    %C = M.*(M+1);
    terms = 2*M+1;

    for m = M'
        l = legendre(m,-incident(1)/r);
        pn(m,:)=l(1,:);
    end
    
    hank = sphh(M, k*r);
    bess = sphj(M, k*a);
    smooth_factor = bess./sphh(M,k*a);
    ith = 1i.^(M.'+1);

    scatter_pressure = -sum(terms.*ith'.*smooth_factor.*hank.*pn);
end

%spherical bessel of the first kind
function j=sphj(n,x)
    j=sqrt(pi/(2*x))*besselj(n+0.5,x);
end

%spherical hard bessel bessel
function jp=sphjp(n,x)
    jp=n/x.*sphj(n,x)-sphj(n+1,x);
end

%spherical hankel
function h=sphh(n,x)
    h=sqrt(pi/(2*x))*besselh(n+0.5,x);
end

%spherical hard hankel
function hp=sphhp(n,x)
    hp=n/x*sphh(n,x)-sphh(n+1,x);
end


