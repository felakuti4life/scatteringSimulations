%PLANE WAVE HITTING AN ELLIPSE
close all; clear all;
%% CONFIG
a = 1.6; %ellipse hyperbolic radius
ep_c = 0.8; %ellipse scale factor
c = 3; %speed of sound (uniform)
K = 16; %number of harmonics we're modelin'

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

%% Get our model for all wave numbers on a single 2d plane
disp('generating scattering model...')
for iter = 1:K
    disp(sprintf('generating harmonic %i of %i...', iter, K));
    k = 2*pi*iter/(length(x)*x_res);
    for l = 1:length(y)
        for m = 1:length(x)
            [zeta,lambda, phi] = cart2obl(x(m), 0, y(l), ep_c);
            if(zeta < a)
                scatter_wave(l,m,iter) = 0;
                inc_wave(l,m,iter) = 0;
            else
                scatter_wave(l,m,iter) =  elliptical_scatter([x(m), y(l),0], k, a, ep_c);
                inc_wave(l,m,iter) = exp(1i*k*X(l,m)+ 1i * pi);
            end
        end
    end
end
%% make our movie

inc_cm = colormap('winter');
scat_cm = colormap('autumn');
t = 0:1/fs:dur*fs;
%inc_map = zeros(length(y),length(x), length(z));
%scat_map = zeros(length(y),length(x), length(z));
disp('generating animation...')
for n = 1:dur*fs
    disp(sprintf('generating frame %i of %i (t = %f)...', n, dur*fs, t(n)));
    %inc_map = zeros(length(y),length(x),length(z));
    %scat_map = zeros(length(y),length(x),length(z));
    inc_map = zeros(length(y), length(x));
    scat_map = zeros(length(y), length(x));
    %sum our waves
    for iter = 1:K
        k = iter*2*pi/(length(x)*x_res);
        inc_map = inc_map + real(inc_wave(:,:,iter).*exp(-1i*k*c*t(n)));
        scat_map = scat_map + real(scatter_wave(:,:,iter).*exp(-1i*k*c*t(n)));
    end
    %close all;
    %figure;
    
    %we get a lot of small numbers so lets just do this
    %inc_map(inc_map<(1/(2^16))) = 0;
    %scat_map(scat_map<(1/(2^16))) = 0;
    
    colormap(hsv);
    alphamap('rampup');
    %plot incidents
    %mesh(X,Y,inc_map);
    
    %scatter3(X(:),Y(:),Z(:), inc_map(:)+0.0000001,inc_map(:)+0.0000000001);
    %hold on
    
    %plot scatters
    %mesh(X,Y, scat_map*3);
    %scatter3(X(:),Y(:), Z(:), scat_map(:) + 0.000001, scat_map(:)+0.0000000001);
    hold on
    imagesc(scat_map, [-20, 20]);%, [-50, 50]);
    
    %imagesc(inc_map, [-20, 20]);
    hold off
    %axis equal;
    %plot a circle for cylinder reference
    %plot(a*sin(-pi:0.01:pi), a*cos(-pi:0.01:pi), 'k');
    %[sphx, sphy, sphz] = sphere(40);
    %surf(sphx*a, sphy*a, sphz*a);
    %axis([-x_dist, x_dist, -y_dist, y_dist, -x_dist*2, x_dist*4]);
    %caxis([-8 20]);
    %view(-12,84);
    title(sprintf('Elliptical Scattering of Plane Wave at %f seconds: a = %f', t(n), a));
    %hold off
    M(n) = getframe(gcf);
    pause(1);
end
disp('bouncing to avi...');
movie2avi(M,'waveHitsEllipse-a1dot6.avi');