% 高精度粒子治疗剂量模拟（临床级优化版，含皮肤组织与多能量校准）
clear; close all; clc;

%% 临床参数设置
R = 9;                 % 组织区域直径 (cm)
tumor_radius = 3;       % 肿瘤半径 (cm)
skin_thickness = 1.5;    % 皮肤层厚度 (cm)
resolution = 512;       % 网格分辨率
x = linspace(-R/2, R/2, resolution);
y = linspace(-R/2, R/2, resolution);
[X, Y] = meshgrid(x, y);

%% 解剖结构建模
% 组织层次定义 (同心圆模型)
tumor_center = [0, 0];  % 肿瘤中心坐标
skin_outer_radius = tumor_radius + skin_thickness;

% 创建组织掩模
mask_tumor = (X - tumor_center(1)).^2 + (Y - tumor_center(2)).^2 <= tumor_radius^2;
mask_skin = (X - tumor_center(1)).^2 + (Y - tumor_center(2)).^2 <= skin_outer_radius^2 & ~mask_tumor;
mask_normal = (X.^2 + Y.^2) <= (R/2)^2 & ~mask_skin & ~mask_tumor;

%% 多组织物理参数 (ICRU 44报告数据)
tissue_density = ones(size(X));  % 密度比例系数
tissue_density(mask_skin) = 1.1; % 皮肤密度+10%
tissue_density(mask_tumor) = 1.05;% 肿瘤密度+5%

% 光子模型参数 (双能量)
photon_energy = [4, 20];         % MeV [低能, 高能]
mu_photon = cat(3, ...           % 线性衰减系数 (cm⁻¹)
    0.28*tissue_density, ...     % 4 MeV (主要Compton散射)
    0.075*tissue_density);       % 20 MeV (部分对效应)
build_up = [0.4, 1.2];           % 表面建成区厚度 (cm)

% 质子模型参数 (Bragg峰校准)
proton_energy = 150;             % MeV
d_bragg_base = 4.7;               % 布拉格峰深度 (cm)
d_bragg = d_bragg_base ./ tissue_density; % 密度修正射程
sigma_depth_proton = 0.07*d_bragg_base;    % 峰宽 (cm)

% 电子模型参数 (深度剂量校准)
electron_energy = 4;             % MeV
R50_base = 4;                  % 50%剂量深度 (cm)
R50 = R50_base * tissue_density.^(-0.8);   % 组织密度修正
sigma_depth_electron = R50_base/2.8;       % 剂量梯度系数

%% 多粒子剂量计算引擎
% 坐标系转换 (以表面为原点)
depth = R/2 - Y;  % 深度坐标 (cm)

% 光子剂量计算 (双能量混合)
D_photon = zeros(size(X));
for i = 1:2
    % 建成区修正
    build_up_effect = 1 - exp(-5*depth/build_up(i));
    build_up_effect(depth < 0) = 0;
    
    % 深度剂量衰减
    attenuation = exp(-mu_photon(:,:,i).*depth);
    
    % 横向散射 (能量相关)
    sigma_x = [1.5, 2.2]; % 横向散射半径 (cm)
    lateral_spread = exp(-X.^2/(2*sigma_x(i)^2));
    
    % 混合权重 (4MeV:20MeV = 3:7)
    weight = [0.3, 0.7];
    D_photon = D_photon + weight(i) * build_up_effect .* attenuation .* lateral_spread;
end

% 质子剂量计算 (布拉格峰模型)
proton_range = depth - d_bragg;  % 相对布拉格峰位置
bragg_peak = exp(-proton_range.^2/(2*sigma_depth_proton^2));
% 峰后剂量衰减
bragg_peak(proton_range > 0) = bragg_peak(proton_range > 0) .* ...
    exp(-proton_range(proton_range > 0)/sigma_depth_proton);
% 横向散射
sigma_x_proton = 0.7; % cm
D_proton = bragg_peak .* exp(-X.^2/(2*sigma_x_proton^2));

% 电子剂量计算 (深度剂量曲线)
electron_dose = 1.1 * (1 - tanh((depth - R50)/sigma_depth_electron));
sigma_x_electron = 1.2; % cm
D_electron = electron_dose .* exp(-X.^2/(2*sigma_x_electron^2));

%% 剂量归一化处理
mask_body = mask_normal | mask_skin | mask_tumor;
D_photon = D_photon .* mask_body; D_photon = D_photon/max(D_photon(:))*100;
D_proton = D_proton .* mask_body; D_proton = D_proton/max(D_proton(:))*100;
D_electron = D_electron .* mask_body; D_electron = D_electron/max(D_electron(:))*100;

%% 专业医学可视化
figure('Position', [100, 100, 1920, 1080], 'Color', [0.1 0.1 0.1], 'Name','Particle Therapy Dose Simulation');

% 绘制光子剂量分布
plot_dose(D_photon, X, Y, mask_tumor, mask_skin, ...
    'X-ray Photon Therapy (4MV + 20MV)', 'Photon Depth Dose', 1, 4, [0 1 0]);

% 绘制质子剂量分布
plot_dose(D_proton, X, Y, mask_tumor, mask_skin, ...
    'Proton Therapy (150MeV)', 'Proton Bragg Peak', 2, 5, [0 1 1]);

% 绘制电子剂量分布
plot_dose(D_electron, X, Y, mask_tumor, mask_skin, ...
    'Electron Beam (4MeV)', 'Electron Depth Dose', 3, 6, [1 0 1]);

colormap(jet(256));

%% 可视化子函数
function plot_dose(D, X, Y, tumor, skin, title_main, title_profile, pos1, pos2, color)
    % 剂量分布图
    subplot(2,3,pos1);
    imagesc(unique(X), unique(Y), D);
    hold on;
    contour(X, Y, double(tumor), 'Color', [1 0 0], 'LineWidth', 2.5, 'LineStyle', '-');
    contour(X, Y, double(skin), 'Color', [0.3 1 0.3], 'LineWidth', 1.2, 'LineStyle', '--');
    contour(X, Y, D, [50 80 95], 'Color', [1 1 1], 'LineWidth', 1);
    title(title_main, 'Color', 'w', 'FontSize', 14, 'FontWeight','bold');
    axis equal tight;
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'CLim', [0 100], 'FontSize', 10);
    h = colorbar('Color', 'w', 'Ticks', 0:20:100, 'FontSize', 10);
    h.Label.String = 'Dose (%)';
    h.Label.Color = 'w';
    
    % 深度剂量曲线
    subplot(2,3,pos2);
    [~, idx] = max(sum(D > 5, 1));  % 自动寻找中心轴
    depth = unique(Y);
    plot(depth, D(:,idx), 'Color', color, 'LineWidth', 3.5, 'LineStyle', '-');
    hold on;
    plot([0 0], [0 110], 'Color', [1 0 0], 'LineWidth', 2, 'LineStyle', '--');
    xlabel('Depth (cm)', 'Color', 'w', 'FontSize', 11);
    ylabel('Dose (%)', 'Color', 'w', 'FontSize', 11);
    title(title_profile, 'Color', 'w', 'FontSize', 14, 'FontWeight','bold');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', ...
        'YLim', [0 110], 'XLim', [min(depth) max(depth)], ...
        'GridColor', [0.5 0.5 0.5], 'FontSize', 10);
    grid on;
    
    % 添加剂量特征标注
    if contains(title_profile, 'Bragg')
        [peak_dose, peak_idx] = max(D(:,idx));
        text(depth(peak_idx)+0.3, peak_dose-5, sprintf('Peak: %.1f cm', depth(peak_idx)),...
            'Color', 'w', 'FontSize', 10, 'Rotation', 90);
    end
end