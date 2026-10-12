%% =========================================================================
%  谱空间中的二维浅水方程（Shallow Water Equations, SWE）求解器
% =========================================================================
%  包含 6.4.4 节地转适应两大实验：
%    p.exp_case = 1: 实验 1（初始高度场扰动，式 6.107）-> 高度向风场适应
%    p.exp_case = 2: 实验 2（初始涡度场扰动，式 6.108）-> 高度向风场适应
%  数值方法：
%    - 涡度-散度-扰动高度 (zeta, delta, h') 谱方法表象
%    - 半隐式蛙跳时间推进 (Semi-Implicit Leap-Frog, 式 6.89 - 6.92)
%    - Asselin 滤波器 (gamma = 0.06, 式 6.100)
%    - 2/3 截断反走样滤波 (Orszag 2/3 Rule, 6.2.2 节)
%    - 径向双曲正切海绵阻尼层 (Sponge Damping Layer, 式 6.105)
% =========================================================================

clc; close all;
if ~exist('p', 'var'), p = struct(); end

if ~isfield(p, 'exp_case'), p.exp_case = 1; end % 1: 初始高度扰动; 2: 初始涡度扰动 (表 6)

p.N1 = 128;               % x 方向网格点数 (表 6)
p.N2 = 128;               % y 方向网格点数 (表 6)
p.Lx = 2 * pi;            % x 方向域长度 (表 6)
p.Ly = 2 * pi;            % y 方向域长度 (表 6)

p.f = 1.0;                % 科氏参数 Ro^-1 (表 6)
p.H = 1.0;                % 背景水深 (表 6)
p.g_star = 1.0;           % 约化重力 Fr^-2 (表 6)
p.nu = 0.0;               % 粘性系数 Re^-1 = 0 (表 6)

p.Lv = 0.2;               % 扰源半径 (表 6)
p.alpha_sp = 40.0;        % 特征阻尼率 (表 6)
p.r_sp = 2.5;             % 阻尼区半径 (表 6)
p.d_sp = 1.5;             % 阻尼区过渡厚度 (表 6)

p.gamma = 0.06;           % Asselin 滤波系数 (表 6)
p.dt = 0.0014;            % 时间步长 (表 6)
p.T = 25.0;               % 总模拟时间 (表 6)
if ~exist('p', 'var') || ~isfield(p, 'Step') || isempty(p.Step)
    p.Step = round(p.T / p.dt);
end

p.Ah = -0.3;              % 初始高度扰动幅值 (表 6)
p.Av = 0.2;               % 初始涡度扰动幅值 (表 6)

p.output_fre = 200;       % 绘图刷新频率 (时步)

%% 2. 空间网格与物理坐标 (双周期域 [-pi, pi))
xpos = linspace(-p.Lx/2, p.Lx/2, p.N1 + 1); xpos = xpos(1:p.N1);
ypos = linspace(-p.Ly/2, p.Ly/2, p.N2 + 1); ypos = ypos(1:p.N2);
[X, Y] = ndgrid(xpos, ypos);
R = sqrt(X.^2 + Y.^2);    % 距中心距离

alpha = zeros(p.N1, p.N2);
mask_sp = (R > p.r_sp);
alpha(mask_sp) = p.alpha_sp * tanh((R(mask_sp) - p.r_sp) / p.d_sp);

%% 4. 谱空间波数与算子构建
k_array = zeros(p.N1, p.N2);
kx_array = zeros(p.N1, p.N2);
ky_array = zeros(p.N1, p.N2);
for k1 = 1:p.N1
    for k2 = 1:p.N2
        if (k1 <= p.N1/2)  kx = k1 - 1; else  kx = k1 - 1 - p.N1; end
        if (k2 <= p.N2/2)  ky = k2 - 1; else  ky = k2 - 1 - p.N2; end
        kx_array(k1, k2) = 2*pi/p.Lx * kx;
        ky_array(k1, k2) = 2*pi/p.Ly * ky;
        k_array(k1, k2) = ((2*pi/p.Lx * kx)^2 + (2*pi/p.Ly * ky)^2)^0.5;
    end
end
k_sq = k_array.^2;

% 逆拉普拉斯算子 1/k^2 (k=0 处设为 0)
inv_ksq = zeros(p.N1, p.N2);
inv_ksq(k_sq > 0) = 1.0 ./ k_sq(k_sq > 0);

% 2/3 截断反走样算子 (讲义 6.2.2 节)
dealias = ones(p.N1, p.N2);
k_max_x = max(abs(kx_array(:)));
k_max_y = max(abs(ky_array(:)));
dealias(abs(kx_array) > (2/3)*k_max_x) = 0;
dealias(abs(ky_array) > (2/3)*k_max_y) = 0;

zeta0  = zeros(p.N1, p.N2);
delta0 = zeros(p.N1, p.N2);
h0     = zeros(p.N1, p.N2);

if p.exp_case == 1
    % 实验 1: 仅初始高度场扰动 (式 6.107)
    h0 = p.Ah * exp(-(X.^2 + Y.^2) / p.Lv^2) - pi * p.Ah * (p.Lv^2) / (p.Lx * p.Ly);
    zeta0 = zeros(p.N1, p.N2);
    delta0 = zeros(p.N1, p.N2);
    exp_title = 'Exp 1: Initial Height Perturbation (Ah = -0.3)';
elseif p.exp_case == 2
    % 实验 2: 仅初始风场(涡度)扰动 (式 6.108)
    zeta0 = p.Av * exp(-(X.^2 + Y.^2) / p.Lv^2) - pi * p.Av * (p.Lv^2) / (p.Lx * p.Ly);
    h0 = zeros(p.N1, p.N2);
    delta0 = zeros(p.N1, p.N2);
    exp_title = 'Exp 2: Initial Vorticity Perturbation (Av = 0.2)';
end

% 转换初始场到谱空间 (除以 N1*N2)
z_h_n = fft2(zeta0)  / (p.N1 * p.N2);
d_h_n = fft2(delta0) / (p.N1 * p.N2);
h_h_n = fft2(h0)     / (p.N1 * p.N2);

% 蛙跳格式需要 n-1 层变量
z_h_nm1 = z_h_n;
d_h_nm1 = d_h_n;
h_h_nm1 = h_h_n;

%% 6. 半隐式矩阵分母预计算 (式 6.92 - 6.98)
% 隐式重力波行列式分母: 1 + g* * H * k^2 * dt^2
denom = 1.0 + p.g_star * p.H * k_sq * (p.dt^2);

%% 7. 绘图窗口初始化
figure('Color', 'w', 'Position', [100, 100, 1050, 480]);

% 预画 Rossby 变形半径圆 (半径 = 1.0) 和阻尼区圆 (半径 = 2.5) 的角度参数
theta_c = linspace(0, 2*pi, 200);
Lr = sqrt(p.g_star * p.H) / p.f; % Rossby 变形半径 = 1.0
x_Lr = Lr * cos(theta_c); y_Lr = Lr * sin(theta_c);
x_sp = p.r_sp * cos(theta_c); y_sp = p.r_sp * sin(theta_c);

fprintf('====================================================\n');
fprintf('谱空间浅水模式启动: %s\n', exp_title);
fprintf('网格: %dx%d, 总步数: %d, 时间步长 dt = %.4f, 总时间 T = %.1f\n', ...
        p.N1, p.N2, p.Step, p.dt, p.T);
fprintf('====================================================\n');

tic;

%% 8. 主时间推进循环 (半隐式蛙跳 + Asselin 滤波)
for step = 1:p.Step
    t_now = step * p.dt;

    % 计算显式非重力波右端项 P_zeta, P_delta, P_h (式 6.81 - 6.83)
    [Pz_hat, Pd_hat, Ph_hat] = calc_explicit_P(z_h_n, d_h_n, h_h_n, ...
        kx_array, ky_array, inv_ksq, dealias, alpha, p.f, p.N1, p.N2);

    if step == 1
        % 第 1 步使用前向欧拉启动 (dt)
        z_h_np1 = z_h_n + p.dt * Pz_hat;
        
        R_d = d_h_n + p.dt * Pd_hat + p.g_star * k_sq * p.dt .* h_h_n;
        R_h = h_h_n + p.dt * Ph_hat - p.H * p.dt .* d_h_n;
        
        d_h_np1 = (R_d + p.g_star * k_sq * p.dt .* R_h) ./ denom;
        h_h_np1 = (R_h - p.H * p.dt .* R_d) ./ denom;
    else
        % 第 2 步及以后采用半隐式蛙跳步进 (2*dt, 式 6.89 - 6.92)
        z_h_np1 = z_h_nm1 + 2.0 * p.dt * Pz_hat;

        R_d = d_h_nm1 + 2.0 * p.dt * Pd_hat + p.g_star * k_sq * p.dt .* h_h_nm1;
        R_h = h_h_nm1 + 2.0 * p.dt * Ph_hat - p.H * p.dt .* d_h_nm1;

        d_h_np1 = (R_d + p.g_star * k_sq * p.dt .* R_h) ./ denom;
        h_h_np1 = (R_h - p.H * p.dt .* R_d) ./ denom;

        % Asselin 滤波器平滑 n 层计算模态 (式 6.100)
        z_h_n = z_h_n + p.gamma * (z_h_nm1 - 2.0 * z_h_n + z_h_np1);
        d_h_n = d_h_n + p.gamma * (d_h_nm1 - 2.0 * d_h_n + d_h_np1);
        h_h_n = h_h_n + p.gamma * (h_h_nm1 - 2.0 * h_h_n + h_h_np1);
    end

    % 更新时间层: n-1 <- n, n <- n+1
    z_h_nm1 = z_h_n;
    d_h_nm1 = d_h_n;
    h_h_nm1 = h_h_n;

    z_h_n = z_h_np1;
    d_h_n = d_h_np1;
    h_h_n = h_h_np1;

    % 定期绘图输出 (图 6.13: 涡度 zeta/f 与总高度 h/H)
    if mod(step, p.output_fre) == 0 || step == p.Step
        zeta_phys = real(ifft2(z_h_n) * (p.N1 * p.N2));
        h_phys    = real(ifft2(h_h_n) * (p.N1 * p.N2));
        total_h   = p.H + h_phys; % 总水深 h = H + h'

        % 左图: 相对涡度 zeta / f
        subplot(1, 2, 1);
        pcolor(xpos, ypos, (zeta_phys / p.f)'); shading interp; axis square;
        xlabel('x'); ylabel('y'); colormap(jet); colorbar;
        hold on;
        plot(x_Lr, y_Lr, 'k:', 'LineWidth', 1.8);  % 点线标记 Rossby 变形半径
        plot(x_sp, y_sp, 'k--', 'LineWidth', 1.5); % 虚线标记阻尼层边界
        hold off;
        title(sprintf('Relative Vorticity \\zeta / f at t = %.2f', t_now), 'FontSize', 12);

        % 右图: 总高度 h / H
        subplot(1, 2, 2);
        pcolor(xpos, ypos, (total_h / p.H)'); shading interp; axis square;
        xlabel('x'); ylabel('y'); colormap(jet); colorbar;
        hold on;
        plot(x_Lr, y_Lr, 'k:', 'LineWidth', 1.8);
        plot(x_sp, y_sp, 'k--', 'LineWidth', 1.5);
        hold off;
        title(sprintf('Total Height h / H at t = %.2f', t_now), 'FontSize', 12);

        sgtitle(sprintf('%s (t = %.2f / %.1f)', exp_title, t_now, p.T), 'FontSize', 14);
        drawnow;
    end
end

elapsed_time = toc;
fprintf('\n浅水模式运行完成！耗时: %.2f 秒\n', elapsed_time);

%% =========================================================================
%  显式右端项 P 计算子函数 (式 6.81 - 6.83)
% =========================================================================
function [Pz_hat, Pd_hat, Ph_hat] = calc_explicit_P(zh, dh, hh, ...
    kx, ky, inv_ksq, dealias, alpha, f, N1, N2)

    % 1. Helmholtz 分解求速度 u, v (谱空间代数表达式)
    %    u_hat = -i*kx/k^2 * delta_hat + i*ky/k^2 * zeta_hat
    %    v_hat = -i*ky/k^2 * delta_hat - i*kx/k^2 * zeta_hat
    u_hat = (-1j * kx .* inv_ksq .* dh + 1j * ky .* inv_ksq .* zh) .* dealias;
    v_hat = (-1j * ky .* inv_ksq .* dh - 1j * kx .* inv_ksq .* zh) .* dealias;

    % 2. 谱微分求空间梯度
    zx_hat = (1j * kx .* zh) .* dealias;
    zy_hat = (1j * ky .* zh) .* dealias;
    dx_hat = (1j * kx .* dh) .* dealias;
    dy_hat = (1j * ky .* dh) .* dealias;
    hx_hat = (1j * kx .* hh) .* dealias;
    hy_hat = (1j * ky .* hh) .* dealias;

    % 3. 变换至物理空间
    u  = real(ifft2(u_hat) * (N1 * N2));
    v  = real(ifft2(v_hat) * (N1 * N2));
    z  = real(ifft2(zh .* dealias) * (N1 * N2));
    d  = real(ifft2(dh .* dealias) * (N1 * N2));
    h  = real(ifft2(hh .* dealias) * (N1 * N2));

    zx = real(ifft2(zx_hat) * (N1 * N2));
    zy = real(ifft2(zy_hat) * (N1 * N2));
    dx = real(ifft2(dx_hat) * (N1 * N2));
    dy = real(ifft2(dy_hat) * (N1 * N2));
    hx = real(ifft2(hx_hat) * (N1 * N2));
    hy = real(ifft2(hy_hat) * (N1 * N2));

    % 4. 物理空间非线性与强迫项 (严格对应式 6.81 - 6.83 与 6.104)
    %    P_zeta  = -u*zx - v*zy - (f + z)*d - alpha*z
    %    P_delta = -u*dx - v*dy + (f + z)*z - alpha*d
    %    P_h     = -h*d  - u*hx - v*hy      - alpha*h
    P_z = -u .* zx - v .* zy - (f + z) .* d - alpha .* z;
    P_d = -u .* dx - v .* dy + (f + z) .* z - alpha .* d;
    P_h = -h .* d  - u .* hx - v .* hy      - alpha .* h;

    % 5. 变回谱空间并施加 2/3 截断反走样
    Pz_hat = (fft2(P_z) / (N1 * N2)) .* dealias;
    Pd_hat = (fft2(P_d) / (N1 * N2)) .* dealias;
    Ph_hat = (fft2(P_h) / (N1 * N2)) .* dealias;
end
