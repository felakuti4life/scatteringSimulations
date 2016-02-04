%PLANE WAVE HITTING A SPHERE
close all; clear all;
%% CONFIG
a = 1; %cylinder radius
c = 3; %speed of sound (uniform)
K = 32; %number of harmonics we're modelin'

%Actual space
x_dist = 3;
y_dist = 3;
x_res = 0.05;
y_res = 0.05;


%Actual time
fs = 20;
dur = 6;


%% allocation
disp('allocating matrices...')
x = -x_dist:x_res:x_dist;
y = -y_dist:y_res:y_dist;
[X,Y] = meshgrid(x,y);
scatter_wave = zeros(length(y), length(x), K);
inc_wave = zeros(length(y), length(x), K);

%% Get our model for all wave numbers
disp('generating scattering model...')
for iter = 1:K
    disp(sprintf('generating harmonic %i of %i...', iter, K));
    k = 2*pi*iter/(length(x)*x_res);
    for l = 1:length(y)
        for m = 1:length(x)
            if((x(m)^2+y(l)^2)^0.5 < a)
                scatter_wave(l,m,iter) = 0;
                inc_wave(l,m,iter) = 0;
            else
                scatter_wave(l,m,iter) = spherical_scatter([x(m), y(l),0], k, a);
                inc_wave(l,m,iter) = exp(1i*k*X(l,m));
            end
        end
    end
end

%% make our movie

inc_cm = colormap('winter');
scat_cm = colormap('autumn');
t = 0:1/fs:dur*fs;
inc_map = zeros(length(y),length(x));
scat_map = zeros(length(y),length(x));
disp('generating animation...')
for n = 1:dur*fs
    disp(sprintf('generating frame %i of %i (t = %f)...', n, dur*fs, t(n)));
    inc_map = zeros(length(y),length(x));
    scat_map = zeros(length(y),length(x));
    %sum our waves
    for iter = 1:K
        k = iter*2*pi/(length(x)*x_res);
        inc_map = inc_map + real(inc_wave(:,:,iter).*exp(-1i*k*c*t(n)));
        scat_map = scat_map + real(scatter_wave(:,:,iter).*exp(-1i*k*c*t(n)));
    end
    %close all;
    %figure;
    
    colormap(1 - gray);
    %plot incidents
    
    mesh(X,Y, inc_map);
    hold on
    %plot scatters
    
    mesh(X,Y, scat_map);
    %imagesc(scat_map, [-50, 50]);
    
    %plot a circle for cylinder reference
    plot(a*sin(-pi:0.01:pi), a*cos(-pi:0.01:pi), 'k');
    
    axis([-x_dist, x_dist, -y_dist, y_dist, -x_dist*2, x_dist*4]);
    caxis([-8 20]);
    view(-12,84);
    title(sprintf('Spherical Scattering of Plane Wave at %f seconds', t(n)));
    hold off
    M(n) = getframe(gcf);
    pause(1);
end
disp('bouncing to avi...');
movie2avi(M,'waveHitsSphere.avi');