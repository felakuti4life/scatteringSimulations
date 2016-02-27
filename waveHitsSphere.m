%PLANE WAVE HITTING A SPHERE
close all; clear all;
%% CONFIG
a = 1; %sphere radius
c = 3; %speed of sound (uniform)
K = 32; %number of harmonics we're modelin'

%Actual space
x_dist = 3;
y_dist = 3;
z_dist = 3;

x_res = 0.05;
y_res = 0.05;
z_res = 0.05;

%Actual time
fs = 20;
dur = 6;


%% allocation
disp('allocating matrices...')
x = -x_dist:x_res:x_dist;
y = -y_dist:y_res:y_dist;
z = -z_dist:z_res:z_dist;

[X,Y,Z] = meshgrid(x,y,z);
scatter_wave = zeros(length(y), length(x), length(z), K);
inc_wave = zeros(length(y), length(x), length(z), K);

%% Get our model for all wave numbers on a single 2d plane
disp('generating scattering model...')
zero_idx = find(y==0);
for iter = 1:K
    disp(sprintf('generating harmonic %i of %i...', iter, K));
    k = 2*pi*iter/(length(x)*x_res);
    for l = 1:length(y)
        for m = 1:length(x)            
            if((x(m)^2+y(l)^2)^0.5 < a)
                scatter_wave(l,m,zero_idx,iter) = 0;
                inc_wave(l,m,zero_idx,iter) = 0;
            else
                scatter_wave(l,m,zero_idx,iter) =  spherical_scatter([x(m), y(l),0], k, a);
                inc_wave(l,m,zero_idx,iter) = exp(1i*k*X(l,m,zero_idx)+ 1i * pi);
            end
        end
    end
end
inc_map = zeros(length(y), length(x));
%% rotate 2D map to get 3D model
for iter = 1:K
    disp(sprintf('rotaing harmonic %i of %i...', iter, K));
    for l = zero_idx:length(y)
        for n = zero_idx+1:length(z)
            y_idx = floor(length(y)/2 + ((y(l)^2 + z(n)^2)^0.5)/y_res);
            if y_idx > length(y)
                y_idx = length(y);
            end
            for m = zero_idx:length(x)
                s = scatter_wave(y_idx, m, zero_idx, iter);
                i = inc_wave(y_idx,m,zero_idx,iter);
                scatter_wave(l,m,n, iter) = s;
                scatter_wave(zero_idx - l,m,zero_idx - n, iter) =  s;
                scatter_wave(zero_idx - l,m,n, iter) =  s;
                scatter_wave(l,m,zero_idx - n, iter) = s;
                    
                inc_wave(l,m,n,iter) = i;
                inc_wave(zero_idx - l,m,n,iter) = i;
                inc_wave(zero_idx - l,m,zero_idx - n,iter) = i;
                inc_wave(l,m,zero_idx - n,iter) = i;
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
    inc_map = zeros(length(y), length(x),length(z));
    scat_map = zeros(length(y), length(x),length(z));
    %sum our waves
    for iter = 1:K
        k = iter*2*pi/(length(x)*x_res);
        inc_map = inc_map + real(inc_wave(:,:,:,iter).*exp(-1i*k*c*t(n)));
        scat_map = scat_map + real(scatter_wave(:,:,:,iter).*exp(-1i*k*c*t(n)));
    end
    disp([max(scat_map(:)), min(scat_map(:))]);
    %close all;
    %figure;
    
    colormap(hsv);
    alphamap('rampup');
    %plot incidents
    %mesh(X,Y,Z,inc_map);
    %scatter3(X(:),Y(:),Z(:), inc_map(:)+0.0000001,inc_map(:)+0.0000000001);
    %hold on
    %plot scatters
    
    %mesh(X,Y,Z, scat_map);
    %scatter3(X(:),Y(:), Z(:), scat_map(:) + 0.000001, scat_map(:)+0.0000000001);
    imagesc(scat_map(:,:,zero_idx), [-20, 20]);%, [-50, 50]);
    axis equal;
    %plot a circle for cylinder reference
    %plot(a*sin(-pi:0.01:pi), a*cos(-pi:0.01:pi), 'k');
    %[sphx, sphy, sphz] = sphere(40);
    %surf(sphx*a, sphy*a, sphz*a);
    %axis([-x_dist, x_dist, -y_dist, y_dist, -x_dist*2, x_dist*4]);
    %caxis([-8 20]);
    %view(-12,84);
    title(sprintf('Spherical Scattering of Plane Wave at %f seconds', t(n)));
    %hold off
    M(n) = getframe(gcf);
    pause(1);
end
disp('bouncing to avi...');
movie2avi(M,'waveHitsSphere.avi');