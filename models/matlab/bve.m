clc; clear all; close all;

%% 1. 基本参数与全局结构体设定
p.model_case = 1;
p.Lx = 2 * pi;
p.Ly = 2 * pi;

switch p.model_case
    case 1 % 衰减涡格
        p.givenname = 'case1 vortex grid';
        p.N1 = 64;  p.N2 = 64;
        p.Re = 1.0;
        p.T = 6.0;
        p.dt = 3.75e-4;
        p.beta_star = 0.0;
        p.output_fre = 200;

    case 2 % 涡旋相互作用
        p.givenname = 'case2 vortex merger';
        p.N1 = 64;  p.N2 = 64;
        p.Re = 50.0;
        p.T = 6.0;
        p.dt = 8.33e-5;
        p.beta_star = 0.0;
        p.output_fre = 500;
        p.Av = 40.0; p.Lv = 0.4;
        p.xv1 = -1.0; p.yv1 = 0.0;
        p.xv2 =  1.0; p.yv2 = 0.0;

    case 3 % 剪切不稳定
        p.givenname = 'case3 shear instability';
        p.N1 = 64;  p.N2 = 64;
        p.Re = 100.0;
        p.T = 6.0;
        p.dt = 8.33e-5;
        p.beta_star = 0.0;
        p.output_fre = 500;
        p.Ab = 40.0; p.Lb = 0.1; p.xi0 = 0.1;

    case 4 % 涡旋 beta 漂移
        p.givenname = 'case4 beta drift';
        p.N1 = 64;  p.N2 = 64;
        p.Re = 50.0;
        p.T = 6.0;
        p.dt = 8.33e-5;
        p.beta_star = 1.0;
        p.output_fre = 500;
        p.Av = 40.0; p.Lv = 0.4;
        p.xv1 = 0.0; p.yv1 = 0.0;

    case 5 % 二维衰减湍流
        p.givenname = 'case5 2d turbulence';
        p.N1 = 256; p.N2 = 256;
        p.Re = 4.0e4;
        p.T = 300.0;
        p.dt = 0.05;
        p.beta_star = 0.0;
        p.output_fre = 50;
        p.kc = 2 * pi / 0.5;
        p.dk = 2 * pi / 2.0;
end

p.nu = 1.0 / p.Re;
p.Step = round(p.T / p.dt);

%% 2. 空间网格坐标声明
xpos = linspace(-p.Lx/2, p.Lx/2, p.N1 + 1);
xpos = xpos(1:p.N1);
ypos = linspace(-p.Ly/2, p.Ly/2, p.N2 + 1);
ypos = ypos(1:p.N2);
[X, Y] = ndgrid(xpos, ypos);

%% 3. 谱空间波数坐标声明
k_array = zeros(p.N1, p.N2);
kx_array = zeros(p.N1, p.N2);
ky_array = zeros(p.N1, p.N2);
for k1 = 1:p.N1
    for k2 = 1:p.N2
        if (k1 <= p.N1/2)
            kx = k1 - 1;
        else
            kx = k1 - 1 - p.N1;
        end
        if (k2 <= p.N2/2)
            ky = k2 - 1;
        else
            ky = k2 - 1 - p.N2;
        end
        kx_array(k1, k2) = 2*pi/p.Lx * kx;
        ky_array(k1, k2) = 2*pi/p.Ly * ky;
        k_array(k1, k2) = ((2*pi/p.Lx * kx)^2 + (2*pi/p.Ly * ky)^2)^0.5;
    end
end

inv_ksq = zeros(p.N1, p.N2);
inv_ksq(k_array > 0) = -1.0 ./ (k_array(k_array > 0).^2);

[p.dealias] = dealias_make(kx_array, ky_array);

diff_factor = exp(-p.nu * (k_array.^2) * p.dt);

%% 4. 初始条件设定
Iz0 = zeros(p.N1, p.N2);
switch p.model_case
    case 1 % 衰减涡格
        Iz0 = -2.0 * cos(X) .* cos(Y);

    case 2 % 涡旋相互作用
        Iz0 = p.Av * exp(-((X - p.xv1).^2 + (Y - p.yv1).^2) / p.Lv^2) + ...
              p.Av * exp(-((X - p.xv2).^2 + (Y - p.yv2).^2) / p.Lv^2) - ...
              2 * pi * p.Av * (p.Lv^2) / (p.Lx * p.Ly);

    case 3 % 剪切不稳定
        Iz0 = p.Ab * exp(-Y.^2 / p.Lb^2) + p.xi0 * randn(p.N1, p.N2);

    case 4 % 涡旋 beta 漂移
        Iz0 = p.Av * exp(-((X - p.xv1).^2 + (Y - p.yv1).^2) / p.Lv^2);

    case 5 % 二维衰减湍流
        k_array_shift = fftshift(k_array);

        Iz0_hat_amp = exp( - (k_array_shift - p.kc).^2 / p.dk^2 );
        Iz0_hat_shift = zeros(p.N1, p.N2);
        for i = 2: p.N1
            for j = 2: (p.N2/2+1)
                Iz0_hat_shift(i, j) = Iz0_hat_amp(i, j) * exp(1j*2*pi*rand);
            end
        end 

        j = p.N2/2+1;
        for i = (p.N1/2+2):p.N1
            Iz0_hat_shift(i,j) = conj(Iz0_hat_shift(p.N1+2-i,j));
        end
        Iz0_hat_shift(p.N1/2+1,p.N2/2+1) = 0; 

        for i = 2:p.N1
            for j = (p.N2/2+2):p.N2
                if ((p.N1+2-i>=1 && p.N1+2-i<=p.N1) && (p.N2+2-j>=1 && p.N2+2-j<=p.N2))
                    Iz0_hat_shift(i,j) = conj(Iz0_hat_shift(p.N1+2-i,p.N2+2-j));
                end
            end
        end
        Iz0_hat = fftshift(Iz0_hat_shift);

        Iz0 = real(ifft2(Iz0_hat) * (p.N1 * p.N2));
        Iz0 = Iz0 / std(Iz0(:)) * 0.25;
end

Iz0_hat = fft2(Iz0) / (p.N1 * p.N2);
Iz1_hat = Iz0_hat;

%% 5. 记录与可视化初始化
if p.model_case == 1
    KE_record = zeros(p.Step + 1, 1);
    KE_exact  = zeros(p.Step + 1, 1);
    t_record  = zeros(p.Step + 1, 1);
    
    psi0_hat = Iz0_hat .* inv_ksq;
    u0_hat = -psi0_hat .* (1j * ky_array);
    v0_hat =  psi0_hat .* (1j * kx_array);
    KE0 = sum(sum(0.5 * (conj(u0_hat).*u0_hat + conj(v0_hat).*v0_hat)));
    KE_record(1) = real(KE0);
    KE_exact(1)  = real(KE0);
    t_record(1)  = 0;
end

figure(1);
pcolor(xpos, ypos, Iz0'); shading interp; axis square; colorbar;
title(sprintf('%s (Case %d): Initial Vorticity \\zeta_0', p.givenname, p.model_case));

tic;

%% 6. 计算推进主循环
for t = 1:p.Step
    g1 = calc_g(Iz1_hat, inv_ksq, kx_array, ky_array, p.dealias, p.N1, p.N2, p.beta_star);
    g2 = calc_g(Iz1_hat + 0.5 * p.dt * g1, inv_ksq, kx_array, ky_array, p.dealias, p.N1, p.N2, p.beta_star);
    g3 = calc_g(Iz1_hat + 0.5 * p.dt * g2, inv_ksq, kx_array, ky_array, p.dealias, p.N1, p.N2, p.beta_star);
    g4 = calc_g(Iz1_hat + p.dt * g3, inv_ksq, kx_array, ky_array, p.dealias, p.N1, p.N2, p.beta_star);

    Iz2_hat = Iz1_hat + (p.dt / 6.0) * (g1 + 2*g2 + 2*g3 + g4);

    Iz2_hat = Iz2_hat .* diff_factor;

    if p.model_case == 1
        t_now = t * p.dt;
        psi_hat = Iz2_hat .* inv_ksq;
        u_hat = -psi_hat .* (1j * ky_array);
        v_hat =  psi_hat .* (1j * kx_array);
        KE_record(t + 1) = real(sum(sum(0.5 * (conj(u_hat).*u_hat + conj(v_hat).*v_hat))));
        KE_exact(t + 1)  = real(KE0 * exp(-4.0 * p.nu * t_now));
        t_record(t + 1)  = t_now;
    end

    Iz1_hat = Iz2_hat;

    if (mod(t, p.output_fre) == 0 || t == p.Step)
        Iz = real(ifft2(Iz1_hat) * (p.N1 * p.N2));
        t_now = t * p.dt;

        figure(2);
        if p.model_case == 1
            subplot(1, 2, 1);
            pcolor(xpos, ypos, Iz'); shading interp; axis square;
            xlabel('x'); ylabel('y'); colormap(jet); colorbar;
            title(sprintf('Vorticity \\zeta at t = %.3f', t_now));

            subplot(1, 2, 2);
            plot(t_record(1:t+1), KE_record(1:t+1), 'b-', 'LineWidth', 2); hold on;
            plot(t_record(1:t+1), KE_exact(1:t+1), 'r--', 'LineWidth', 1.8); hold off;
            xlabel('Time t'); ylabel('Kinetic Energy');
            title('KE Decay Verification');
            legend('Numerical', 'Exact e^{-4\nu t}');
            grid on; axis square;

        elseif p.model_case == 5
            subplot(1, 2, 1);
            pcolor(xpos, ypos, Iz'); shading interp; axis square;
            xlabel('x'); ylabel('y'); colormap(jet); colorbar;
            title(sprintf('2D Turbulence \\zeta at t = %.1f', t_now));

            subplot(1, 2, 2);
            [k_bins, Ek] = compute_spectrum(Iz1_hat, inv_ksq, kx_array, ky_array, k_array, p.N1);
            loglog(k_bins, Ek, 'b-', 'LineWidth', 2); hold on;
            ref_k = k_bins(k_bins >= 4 & k_bins <= 20);
            if ~isempty(ref_k)
                ref_Ek = Ek(find(k_bins >= 4, 1)) * (ref_k / ref_k(1)).^(-3);
                loglog(ref_k, ref_Ek, 'r--', 'LineWidth', 1.8);
                legend('E(k)', 'k^{-3} theoretical slope', 'Location', 'southwest');
            end
            hold off;
            xlabel('Wavenumber k'); ylabel('Energy E(k)');
            title('Kinetic Energy Spectrum'); grid on; axis square;

        else
            pcolor(xpos, ypos, Iz'); shading interp; axis square;
            xlabel('x'); ylabel('y'); colormap(jet); colorbar;
            title(sprintf('%s (Case %d) at t = %.3f', p.givenname, p.model_case, t_now));
        end
        drawnow;
    end
end

toc;

if p.model_case == 1
    rel_err = max(abs(KE_record - KE_exact) ./ KE_exact);
    fprintf('\n最大相对误差 = %.4e\n', rel_err);
    if rel_err < 1e-3
        fprintf('代码编写正确\n');
    end
end
%%
%  子函数 1: 反走样算子生成
function [dealias] = dealias_make(kx_array, ky_array)
    k_max_x = max(abs(kx_array(:)));
    k_max_y = max(abs(ky_array(:)));
    dealias = ones(size(kx_array));
    dealias(abs(kx_array) > (2/3)*k_max_x) = 0;
    dealias(abs(ky_array) > (2/3)*k_max_y) = 0;
end
%%
%  子函数 2: 右端平流项计算
function g = calc_g(Iz_hat, inv_ksq, kx_array, ky_array, dealias, N1, N2, beta_star)
    psi_hat = Iz_hat .* inv_ksq;

    u_hat = -psi_hat .* (1j * ky_array);
    v_hat =  psi_hat .* (1j * kx_array);

    zx_hat = Iz_hat .* (1j * kx_array);
    zy_hat = Iz_hat .* (1j * ky_array);

    u  = real(ifft2(dealias .* u_hat)  * (N1 * N2));
    v  = real(ifft2(dealias .* v_hat)  * (N1 * N2));
    zx = real(ifft2(dealias .* zx_hat) * (N1 * N2));
    zy = real(ifft2(dealias .* zy_hat) * (N1 * N2));

    J = u .* zx + v .* zy;

    J_hat = (fft2(J) / (N1 * N2)) .* dealias;

    if beta_star ~= 0
        g = -J_hat - beta_star * v_hat;
    else
        g = -J_hat;
    end
end
%%
%  子函数 3: 一维各向同性能谱计算
function [k_bins, Ek] = compute_spectrum(Iz_hat, inv_ksq, kx_array, ky_array, k_array, N1)
    psi_hat = Iz_hat .* inv_ksq;
    u_hat = -psi_hat .* (1j * ky_array);
    v_hat =  psi_hat .* (1j * kx_array);
    
    E_2D = 0.5 * (conj(u_hat).*u_hat + conj(v_hat).*v_hat);
    
    max_k = floor(N1 / 3);
    k_bins = (1:max_k)';
    Ek = zeros(max_k, 1);
    dk = 1.0;
    
    for i = 1:max_k
        k_val = k_bins(i);
        mask = (k_array >= (k_val - dk/2)) & (k_array < (k_val + dk/2));
        Ek(i) = sum(E_2D(mask));
    end
end

