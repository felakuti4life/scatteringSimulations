%CYLINDER SCATTER STUFF
%% params
a = 0.3;
x_dist = 1;
y_dist = 1;
x_res = 0.01;
y_res = 0.01;

c = 348;
freq_min = 5;
freq_max = 600;
freq_res = 5;
%% get scatters
x = -x_dist:x_res:x_dist;
y = -y_dist:y_res:y_dist;
sc = zeros(length(x), length(y));
frm = 1;
for f = freq_min:freq_res:freq_max
    for n = 1:length(x)
        for m = 1:length(y)
            sc(n,m) = cylinder_scatter([x(n), y(m)], (2*pi*f)/c,a);
        end
    end
    hold on
    imagesc(abs(sc));
    M(frm) = getframe(gcf);
    frm = frm +1;
    str = sprintf('cylindrical scattering at %f Hz', f);
    title(str);
    hold off
end

for i = 1:5
    zoom on
    M(frm) = getframe(gcf);
    frm = frm+1;
end

movie2avi(M,'cylindricalScatterCloseUp.avi');

    M(frm) = getframe(gcf);