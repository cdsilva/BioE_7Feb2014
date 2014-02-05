clear all
close all

%% set parameters for saving figures
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontSize',20)

res = '-r300';
fmt = '-djpeg';

dpERK_data = '../../Stas_data/membrane_lengths/oct16.mat';
dpERK_image_dir = '../../Stas_data/membrane_pictures/membrane2/dpERK_staining';
dpERK_membrane_dir = '../../Stas_data/membrane_pictures/membrane2';
data_dir = '../../Stas_data/image_analysis_paper';

%% load data
load(dpERK_data);
dpERK = dpERK_raw;

[m, n] = size(dpERK);

% scramble data alignments
dpERK_unaligned = zeros(size(dpERK));
rng(12345);
rand_offsets = zeros(m,1);
for i=1:m
    rand_offsets(i) = randi(n);
    dpERK_unaligned(i,:) = circshift(dpERK(i,:),[0 rand_offsets(i)]);
end

%% plot images

image_idx = 44;

im1 = imread(sprintf('%s/emb%02d.tif',dpERK_image_dir, image_idx));

im1_dpERK = im1;
im1_dpERK(:,:,1) = im1_dpERK(:,:,2);
im1_dpERK(:,:,2) = 0;
im1_dpERK(:,:,3) = 0;

im1_DI = im1;
im1_DI(:,:,1) = im1_DI(:,:,3);
im1_DI(:,:,2) = 0;

im1_membrane = imread(sprintf('%s/emb%02d.tif',dpERK_membrane_dir, image_idx));
im1_membrane(:,:,2) = imadjust(im1_membrane(:,:,2));

figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(im1_dpERK)
set(gca,'position',[0 0 1 1],'units','normalized')
print('drosophila_dpERK', fmt, res)

figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(im1_DI)
set(gca,'position',[0 0 1 1],'units','normalized')
print('drosophila_DI', fmt, res)

figure;
set(gcf, 'paperposition',[0 0 8 8])
imshow(im1_membrane)
set(gca,'position',[0 0 1 1],'units','normalized')
print('drosophila_membrane', fmt, res)


%% plot circle and line profiles

r = 1;
theta = linspace(0, 2*pi, n+1);
theta = theta(1:end-1);
theta = theta - pi/4;

x = r * cos(theta);
y = r * sin(theta);

figure;
set(gcf, 'paperposition',[0 0 8 8])
scatter(x, y, 2000, dpERK(image_idx,:), '.')
hold on
delta=0.1;
plot([x(1)-delta x(1)+delta], [y(1)+delta y(1)-delta], 'linewidth',6, 'color','k')
text(x(1)+delta+0.05, y(1)-delta-0.05, 'Open circle')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
set(gca, 'visible', 'off')
set(gca,'position',[0 0 1 1],'units','normalized')
print('circle_profile', fmt, res)

figure;
set(gcf, 'paperposition',[0 0 8 8])
imagesc(repmat(dpERK(image_idx,:), 20, 1))
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
set(gca, 'visible', 'off')
set(gca,'position',[0 0 1 1],'units','normalized')
print('line_profile', fmt, res)


%% plot scrambled data

% plot data
figure;
imagesc(dpERK)
ylabel('data point (unordered)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_unordered',fmt,res)

%% plot data ordered by membrane length

% plot data
figure;
imagesc(dpERK_sort)
ylabel('data point (ordered using membrane thickness)','fontsize',16)
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_ordered_membrane',fmt,res)

%% PCA

load(sprintf('%s/pca_figures_figures.mat', data_dir));

% data ordered by PCA
figure;
imagesc(dpERK(I,:))
ylabel('data point (ordered using PCA)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_ordered_PCA',fmt,res)

% correlation
figure;
plot(L(:,1),coeff(:,1),'.')
xlabel('membrane thickness')
ylabel('$\langle x_i, \psi_1 \rangle$','interpreter','latex')
print('PCA_time_corr',fmt,res)
fprintf('PCA Spearman coeff: %2.4f \n', corr(L(:,1),coeff(:,1), 'type','spearman'));

% plot first two PCA coefficients
figure;
scatter(coeff(:,1),coeff(:,2),200,'b','.')
xlabel('$\langle x_i, \psi_1 \rangle$','interpreter','latex')
ylabel('$\langle x_i, \psi_2 \rangle$','interpreter','latex')
axis equal
print('coeff_12_new', fmt, res)

%% DMAPS

load(sprintf('%s/dmaps_figures.mat', data_dir));

% plot first two PCA coefficients
figure;
scatter(coeff(I,1),coeff(I,2),200,1:m,'.')
xlabel('$\langle x_i, \psi_1 \rangle$','interpreter','latex')
ylabel('$\langle x_i, \psi_2 \rangle$','interpreter','latex')
axis equal
print('coeff_12_colored', fmt, res)


figure;
imagesc(dpERK(I,:))
ylabel('data point (ordered using DMAPS)')
xlabel('position')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_ordered_DMAPS',fmt,res)

% correlation
figure;
plot(L(:,1),V_dmaps(:,2),'.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
print('DMAPS_time_corr',fmt,res)

figure;
rankcorr = zeros(m, 2);
rankcorr(L(:,2),1) = 1:m;
rankcorr(I,2) = 1:m;
plot(rankcorr(:,1),rankcorr(:,2),'.')
% hold on
% p = polyfit(rankcorr(:,1),rankcorr(:,2),1);
% plot(rankcorr(:,1),p(1) .* rankcorr(:,1) + p(2))
xlabel('ranking from membrane thickness')
ylabel('ranking from DMAPS')
text(0.65*m, 0.3*m, sprintf('correlation = %2.2f', corr(rankcorr(:,1),rankcorr(:,2))), 'fontsize', 14)
print('DMAPS_rank_corr',fmt,res)

fprintf('DMAPS Spearman coeff: %2.4f \n', corr(L(:,1), V_dmaps(:,2), 'type','spearman'));
sprintf('%d, ', I)


%% alignment

load(sprintf('%s/1d_alignment_figures.mat', data_dir));

% plot scrambled data
figure;
imagesc(dpERK_unaligned)
xlabel('position (unaligned)')
ylabel('data point (unordered)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_unaligned_unordered',fmt, res)


% plot alignment with angular synchronization
figure;
imagesc(dpERK_aligned)
xlabel('position (aligned using angular synchronization)')
ylabel('data point (unordered)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_aligned_unordered',fmt, res)

% plot dmaps ordering
figure;
imagesc(dpERK_aligned(I,:))
ylabel('data point (ordered using DMAPS)')
xlabel('position (aligned using angular synchronization)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_ordered_angsynch',fmt, res)

% plot correlation
figure;
plot(L(:,1), V_dmaps(:,2), '.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
print('angsynch_time_corr',fmt, res)
fprintf('Angular synchronization Spearman coeff: %2.4f \n', corr(L(:,1),V_dmaps(:,2), 'type','spearman'));
sprintf('%d, ', I)

% VDM

% plot aligned data
figure;
imagesc(data2)
xlabel('position (aligned using VDM)')
ylabel('scrambled data')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_aligned_vdm',fmt, res)

% plot aligned and ordered data
figure;
imagesc(data2(idx,:))
ylabel('data point (ordered using VDM)')
xlabel('position (aligned using VDM)')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
print('data_ordered_vdm',fmt, res)


% plot correlation
figure;
plot(L(:,1), embed_coord(:, coord_idx), '.')
xlabel('membrane thickness')
ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
print('VDM_time_corr',fmt, res)

fprintf('VDM (1-D) Spearman coeff: %2.4f \n', corr(L(:,1),embed_coord(:, coord_idx), 'type','spearman'));
sprintf('%d, ', idx)


%% raw dpERK images-- synchronization

load(sprintf('%s/2d_alignment_figures.mat', data_dir));
im_save_idx = [1 20 30 40 52];

figure;
plot(L(:,1),embed_coord(:,coord_idx),'.')
xlabel('membrane thickness')
ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
print('vdm_2d_time_corr',fmt, res)
clf

figure;
plot(L(:,1),embed_coord(:,coord_idx),'.')
hold on
plot(L(I(im_save_idx),1), embed_coord(I(im_save_idx),coord_idx), 'or')
xlabel('membrane thickness')
ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
print('vdm_2d_time_corr2',fmt, res)

figure;
for i=1:length(im_save_idx)
    plot(L(:,1),embed_coord(:,coord_idx),'.')
    hold on
    plot(L(I(im_save_idx(i)),1), embed_coord(I(im_save_idx(i)),coord_idx), 'or')
    xlabel('membrane thickness')
    ylabel(sprintf('$\\langle \\phi_%d, \\phi_%d \\rangle$', embed_idx(1, coord_idx), embed_idx(2, coord_idx)),'interpreter','latex')
    print(sprintf('vdm_2d_time_corr_%d',i),fmt, res)
    clf
end

fprintf('VDM (2-D) Spearman coeff: %2.4f \n', corr(L(:,1),embed_coord(:,coord_idx), 'type','spearman'));
sprintf('%d, ', I)


figure;
for i=1:length(im_save_idx)
    set(gcf, 'paperposition',[0 0 8 8])
    imshow(uint8(image_set_aligned_vdm(:,:,I(im_save_idx(i)))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    % make red colormap
    cm_green = gray;
    cm_green(:,2) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('dpERK_vdm_%d',i),fmt,res)
    clf
end

figure;
for i=1:m
    set(gcf, 'paperposition',[0 0 1 1])
    imshow(uint8(image_set_aligned_vdm(:,:,i)), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    % make red colormap
    cm_green = gray;
    cm_green(:,2) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('dpERK_aligned_%d',i),fmt,res)
    clf
end

figure;
for i=1:m
    set(gcf, 'paperposition',[0 0 1 1])
    imshow(uint8(image_set(:,:,i)), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    % make red colormap
    cm_green = gray;
    cm_green(:,2) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('dpERK_unaligned_%d',i),fmt,res)
    clf
end

% figure;
% for i=1:m
%     set(gcf, 'paperposition',[0 0 1 1])
%     imshow(uint8(image_set_aligned_vdm(:,:,I(i))), 'InitialMagnification', 'fit')
%     set(gca,'position',[0 0 1 1],'units','normalized')
%     % make red colormap
%     cm_green = gray;
%     cm_green(:,2) = 0;
%     cm_green(:,3) = 0;
%     colormap(cm_green)
%     axis off
%     print(sprintf('dpERK_vdm_all_%d',i),fmt,res)
%     clf
% end
% 
% 
% figure;
% for i=1:m
%     set(gcf, 'paperposition',[0 0 1 1])
%     imshow(uint8(image_set(:,:,I(i))), 'InitialMagnification', 'fit')
%     set(gca,'position',[0 0 1 1],'units','normalized')
%     % make red colormap
%     cm_green = gray;
%     cm_green(:,2) = 0;
%     cm_green(:,3) = 0;
%     colormap(cm_green)
%     axis off
%     print(sprintf('dpERK_vdm_unaligned_%d',i),fmt,res)
%     clf
% end

% figure;
% for i=1:m
%     set(gcf, 'paperposition',[0 0 8 8])
%     imshow(uint8(image_set_aligned_vdm(:,:,I(i))), 'InitialMagnification', 'fit')
%     set(gca,'position',[0 0 1 1],'units','normalized')
%     % make red colormap
%     cm_green = gray;
%     cm_green(:,2) = 0;
%     cm_green(:,3) = 0;
%     colormap(cm_green)
%     axis off
%     print(sprintf('aligned_%d',i),fmt,res)
%     clf
% end

%% raw dpERK images--scattering

load(sprintf('%s/2d_scattering_figures.mat',data_dir));
im_save_idx = [1 20 29 40 51];

figure;
plot(L(:,1), V(:,2),'.')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
print('DMAPS_scat_time_corr',fmt,res)

figure;
plot(L(:,1), V(:,2),'.')
hold on
plot(L(I(im_save_idx),1), V(I(im_save_idx), 2), 'or')
xlabel('membrane thickness')
ylabel('$\phi_2$','interpreter','latex')
print('DMAPS_scat_time_corr2',fmt,res)

figure;
for i=1:length(im_save_idx)
    plot(L(:,1), V(:,2),'.')
    hold on
    plot(L(I(im_save_idx(i)),1), V(I(im_save_idx(i)), 2), 'or')
    xlabel('membrane thickness')
    ylabel('$\phi_2$','interpreter','latex')
    print(sprintf('DMAPS_scat_time_corr_%d',i),fmt,res)
    clf
end
fprintf('Scattering transform Spearman coeff: %2.4f \n', corr(L(:,1),V(:,2), 'type','spearman'));
sprintf('%d, ', I)


figure;
for i=1:length(im_save_idx)
    set(gcf, 'paperposition',[0 0 8 8])
    imshow(uint8(image_set(:,:,I(im_save_idx(i)))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    % make green colormap
    cm_green = gray;
    cm_green(:,2) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('dpERK_scat_%d',i),fmt,res)
    clf
end

% figure;
% for i=1:m
%     set(gcf, 'paperposition',[0 0 1 1])
%     imshow(uint8(image_set(:,:,I(i))), 'InitialMagnification', 'fit')
%     set(gca,'position',[0 0 1 1],'units','normalized')
%     % make green colormap
%     cm_green = gray;
%     cm_green(:,2) = 0;
%     cm_green(:,3) = 0;
%     colormap(cm_green)
%     axis off
%     print(sprintf('dpERK_scat_all_%d',i),fmt,res)
%     clf
% end

return
%% membrane images--scattering

load(sprintf('%s/2d_membrane_scattering_figures.mat', data_dir));
im_save_idx = [10 26 40 46 52];

figure;
plot(L(:,1), V(:,idx),'.')
xlabel('membrane thickness')
ylabel(sprintf('$\\phi_%d$',idx),'interpreter','latex')
print('DMAPS_membrane_scat_time_corr',fmt,res)

figure;
plot(L(:,1), V(:,idx),'.')
hold on
plot(L(I(im_save_idx),1), V(I(im_save_idx), idx), 'or')
xlabel('membrane thickness')
ylabel(sprintf('$\\phi_%d$',idx),'interpreter','latex')
print('DMAPS_membrane_scat_time_corr2',fmt,res)


figure;
for i=1:length(im_save_idx)
    plot(L(:,1), V(:,idx),'.')
    hold on
    plot(L(I(im_save_idx(i)),1), V(I(im_save_idx(i)), idx), 'or')
    xlabel('membrane thickness')
    ylabel(sprintf('$\\phi_%d$',idx),'interpreter','latex')
    print(sprintf('DMAPS_membrane_scat_time_corr_%d',i),fmt,res)
    clf
end

fprintf('Scattering transform (membrane) Spearman coeff: %2.4f \n', corr(L(:,1), V(:,idx), 'type','spearman'));

figure;
for i=1:length(im_save_idx)
    set(gcf, 'paperposition',[0 0 8 8])
    imshow(logical(image_set_membrane(:,:,I(im_save_idx(i)))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    axis off
    print(sprintf('membrane_scat_%d',i),fmt,res)
    clf
    
    imshow(uint8(image_set_membrane_raw(:,:,I(im_save_idx(i)))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    % make green colormap
    cm_green = gray;
    cm_green(:,1) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('membrane_scat_raw_%d',i),fmt,res)
    clf
end

figure;
for i=1:m
    set(gcf, 'paperposition',[0 0 1 1])
    imshow(logical(image_set_membrane(:,:,I(i))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    axis off
    print(sprintf('membrane_scat_all_%d',i),fmt,res)
    clf
    
    imshow(uint8(image_set_membrane_raw(:,:,I(i))), 'InitialMagnification', 'fit')
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    % make green colormap
    cm_green = gray;
    cm_green(:,1) = 0;
    cm_green(:,3) = 0;
    colormap(cm_green)
    axis off
    print(sprintf('membrane_scat_raw_all_%d',i),fmt,res)
    clf
end

