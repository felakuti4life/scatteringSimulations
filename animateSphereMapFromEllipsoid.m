close all; clear all;
N = 600;
lambda = linspace(0, 2*pi, N);
for i = 1:N
    mapSphereToEllipsoid(0.3, lambda(i), 0.2);
    pause(0.03)
    M(i) = getframe(gcf);
end
movie2avi(M, 'wrong_angle_draw');