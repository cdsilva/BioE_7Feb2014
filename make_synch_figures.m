function make_synch_figures

clear all
close all

fmt = '-dpng';
res = '-r200';

%%
A = imread('PUsig2.jpg');
n = size(A, 1);
A = A(:, 1:n, :);
A(:, floor(0.8*n):n, :) = 255;
A = circshift(A, [0, -17, 0]);
A = rgb2gray(A);

figure; imshow(A)

%%
nimages = 16;
plot1 = 4;
plot2 = 4;
rot_images = zeros([size(A) nimages], 'uint8');

L = 20;
deg_rot = linspace(0, 360, L+1);
deg_rot = deg_rot(1:end-1);

x = (1:n)-mean(1:n);
y = x;
[X, Y] = meshgrid(x, y);
R2 = X.^2 + Y.^2;
circle_window_idx = find(R2 <= (n/2)^2);
circle_window_idx2 = setdiff((1:n^2)', circle_window_idx);

r = 2.5;
r1 = 3;
delta = 0.5;

%% rotate images
%figure;
for i=1:nimages
    %subplot(plot1, plot2, i)
    Arot = imcomplement(imrotate(imcomplement(A), randsample(deg_rot, 1), 'crop'));
    %imshow(Arot)
    rot_images(:,:,i) = Arot;
end

%% add noise
noisy_images = zeros(size(rot_images), 'uint8');
sigma = 0;
%figure;
for i=1:nimages
    %subplot(plot1, plot2, i)
    Anoisy = imnoise(rot_images(:,:,i), 'gaussian', 0, sigma);
    Anoisy(circle_window_idx2) = 255;
    %imshow(Anoisy)
    noisy_images(:,:,i) = Anoisy;
end

figure;
make_plot(noisy_images, 0, r, r1, delta, false, false, zeros(nimages, 1))
print('drosophila_pics/PU_clean',fmt,res)


%% template

template = imnoise(A, 'gaussian', 0, sigma);
template(circle_window_idx2) = 255;
%figure;
%imshow(template)

opt_angles = zeros(nimages, 1);
%figure;
for i=1:nimages
    min_dist = inf;
    min_idx = 0;
    for j=1:L
        image_diff = double(imcomplement(imrotate(imcomplement(noisy_images(:,:,i)), deg_rot(j), 'crop'))) - double(template);
        temp_err = sum(image_diff(circle_window_idx).^2);
        
        if temp_err < min_dist
            min_dist = temp_err;
            min_idx = j;
        end
    end
    %subplot(plot1, plot2, i)
    %Arot = imcomplement(imrotate(imcomplement(rot_images(:,:,i)), deg_rot(min_idx), 'crop'));
    %imshow(Arot)
    opt_angles(i) = deg_rot(min_idx);
end

figure;
make_plot(rot_images, template, r, r1, delta, false, true, zeros(nimages, 1))
print('drosophila_pics/PU_template1_clean',fmt,res)

figure;
make_plot(rot_images, template, r, r1, delta, false, true, opt_angles)
print('drosophila_pics/PU_template2_clean',fmt,res)


%% add noise
noisy_images = zeros(size(rot_images), 'uint8');
sigma = 5;
%figure;
for i=1:nimages
    %subplot(plot1, plot2, i)
    Anoisy = imnoise(rot_images(:,:,i), 'gaussian', 0, sigma);
    Anoisy(circle_window_idx2) = 255;
    %imshow(Anoisy)
    noisy_images(:,:,i) = Anoisy;
end

figure;
make_plot(noisy_images, 0, r, r1, delta, false, false, zeros(nimages, 1))
print('drosophila_pics/PU_noisy',fmt,res)


%% template

template_clean = A;
template = imnoise(A, 'gaussian', 0, sigma);
template(circle_window_idx2) = 255;
%figure;
%imshow(template)

figure;
for i=1:nimages
    min_dist = inf;
    min_idx = 0;
    for j=1:L
        image_diff = double(imcomplement(imrotate(imcomplement(noisy_images(:,:,i)), deg_rot(j), 'crop'))) - double(template);
        temp_err = sum(image_diff(circle_window_idx).^2);
        
        if temp_err < min_dist
            min_dist = temp_err;
            min_idx = j;
        end
    end
    subplot(plot1, plot2, i)
    Arot = imcomplement(imrotate(imcomplement(rot_images(:,:,i)), deg_rot(min_idx), 'crop'));
    Arot(circle_window_idx2) = 255;
    imshow(Arot)
    opt_angles(i) = deg_rot(min_idx);
end

figure;
make_plot(noisy_images, template, r, r1, delta, false, true, zeros(nimages, 1))
print('drosophila_pics/PU_template1_noisy',fmt,res)

figure;
make_plot(noisy_images, template, r, r1, delta, false, true, opt_angles)
print('drosophila_pics/PU_template2_noisy',fmt,res)

figure;
make_plot(rot_images, template_clean, r, r1, delta, false, true, opt_angles)
print('drosophila_pics/PU_template3_noisy',fmt,res)


%% angular synchronization

H = zeros(2*nimages);
for i1=1:nimages
    for i2=1:i1-1
        min_dist = inf;
        min_idx = 0;
        for j=1:L
            image_diff = double(imcomplement(imrotate(imcomplement(noisy_images(:,:,i1)), deg_rot(j), 'crop'))) - double(noisy_images(:,:,i2));
            temp_err = sum(image_diff(circle_window_idx).^2);
            
            if temp_err < min_dist
                min_dist = temp_err;
                min_idx = j;
            end
        end
        H(2*i1-1:2*i1, 2*i2-1:2*i2) = [cos(deg_rot(min_idx)*(pi/180)) -sin(deg_rot(min_idx)*(pi/180)); sin(deg_rot(min_idx)*(pi/180)) cos(deg_rot(min_idx)*(pi/180))];
    end
end
H = H + H';

R_opt = ang_synch(H, 2);

opt_angles = zeros(nimages, 1);
for i=1:nimages
    opt_angles(i) = atan2(R_opt(2*i, 1), R_opt(2*i-1, 1)) * (180/pi);
end

% figure;
% for i=1:nimages
%     subplot(plot1, plot2, i)
%     Arot = imcomplement(imrotate(imcomplement(rot_images(:,:,i)), opt_angles(i), 'crop'));
%     imshow(Arot)
% end

figure;
make_plot(noisy_images, 0, r, r1, delta, true, false, zeros(nimages, 1))
print('drosophila_pics/PU_angsynch1',fmt,res)


figure;
make_plot(noisy_images, 0, r, r1, delta, true, false, opt_angles)
print('drosophila_pics/PU_angsynch2',fmt,res)

figure;
make_plot(rot_images, 0, r, r1, delta, true, false, opt_angles)
print('drosophila_pics/PU_angsynch3',fmt,res)

function make_plot(image_set, template_image, r, r1, delta, draw_kn, draw_template, opt_angles)

% set(0, 'DefaultFigurePaperSize',[12 12]);
% set(0,'DefaultLineMarkerSize',1)
% set(0,'defaultlinelinewidth', 0.1)
% set(0,'DefaultAxesBox', 'off');
% set(0, 'defaultaxesvisible', 'off')
% set(0, 'defaultfigurecolor', 'w');
% set(0, 'defaultfigurepaperposition', [0 0 1 1])


%set(gcf, 'papersize', [12 12])

nimages = size(image_set, 3);

%set(gcf, 'papersize', [8 8])
%set(gcf, 'paperposition', [0 0 1 1])

plot(r*cos(2*pi*(1:nimages)/nimages),r*sin(2*pi*(1:nimages)/nimages), '.w', 'markersize', 1)
hold on
for i=1:nimages
    x = r1 * cos(2*pi*i/nimages);
    y = r1 * sin(2*pi*i/nimages);
    Arot = imcomplement(imrotate(imcomplement(image_set(:,:,i)), opt_angles(i), 'crop'));
    image([x-delta x+delta], [y-delta y+delta],Arot)
    colormap(gray(256))
    hold on
end
if draw_kn
    for i=1:nimages
        for j=1:i-1
            plot(r * cos(2*pi*[i j]/nimages), r * sin(2*pi*[i j]/nimages), 'linewidth', 0.1);
        end
    end
end
if draw_template
    image([-delta delta], [-delta delta],template_image)
    for i=1:nimages
            plot([r delta] * cos(2*pi*i/nimages), [r delta] * sin(2*pi*i/nimages), 'linewidth', 0.1);
    end
end
axis([-r1-delta r1+delta -r1-delta r1+delta])
set(gca, 'box', 'off')
set(gca, 'YDir', 'reverse');
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca, 'visible', 'off')
set(gcf, 'color', 'w');

ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
set(gca,'units','inches');
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
