# =========================================================================
#  二维正压涡度方程（BVE / 2D N-S）傅里叶伪谱法通用求解器 (Python 3 版本)
# =========================================================================

import numpy as np
import time

def run_simulation(model_case=1):
    Lx = 2.0 * np.pi
    Ly = 2.0 * np.pi
    
    if model_case == 1:
        # 表 1: 衰减涡格
        case_name = 'Case 1: Taylor-Green Vortex Grid'
        N1, N2 = 64, 64
        Re = 1.0
        T = 6.0
        dt = 3.75e-4
        beta_star = 0.0
        plot_interval = 200

    elif model_case == 2:
        # 表 2: 涡旋相互作用
        case_name = 'Case 2: Double Vortex Merger'
        N1, N2 = 64, 64
        Re = 50.0
        T = 6.0
        dt = 8.33e-5
        beta_star = 0.0
        Av, Lv = 40.0, 0.4
        xv1, yv1 = -1.0, 0.0
        xv2, yv2 =  1.0, 0.0
        plot_interval = 500

    elif model_case == 3:
        # 表 3: 剪切不稳定
        case_name = 'Case 3: Shear Instability'
        N1, N2 = 64, 64
        Re = 100.0
        T = 6.0
        dt = 8.33e-5
        beta_star = 0.0
        Ab, Lb = 40.0, 0.1
        xi0 = 0.1
        plot_interval = 500

    elif model_case == 4:
        # 表 4: 涡旋 beta 漂移
        case_name = 'Case 4: Vortex Beta Drift'
        N1, N2 = 64, 64
        Re = 50.0
        T = 6.0
        dt = 8.33e-5
        beta_star = 1.0
        Av, Lv = 40.0, 0.4
        xv1, yv1 = 0.0, 0.0
        plot_interval = 500

    elif model_case == 5:
        # 表 5: 二维衰减湍流
        case_name = 'Case 5: 2D Turbulence Decay'
        N1, N2 = 256, 256
        Re = 4.0e4
        T = 300.0
        dt = 0.05
        beta_star = 0.0
        kc = 2.0 * np.pi / 0.5
        dk = 2.0 * np.pi / 2.0
        plot_interval = 50
    else:
        raise ValueError('Invalid model_case!')

    nu = 1.0 / Re
    total_steps = int(round(T / dt))

    # 空间网格
    x = np.linspace(-Lx/2, Lx/2, N1, endpoint=False)
    y = np.linspace(-Ly/2, Ly/2, N2, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # 谱空间波数
    kx = np.fft.fftfreq(N1, d=Lx / (2.0 * np.pi * N1))
    ky = np.fft.fftfreq(N2, d=Ly / (2.0 * np.pi * N2))
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    K_sq = KX**2 + KY**2
    K_array = np.sqrt(K_sq)

    inv_Ksq = np.zeros_like(K_sq)
    inv_Ksq[K_sq > 0] = -1.0 / K_sq[K_sq > 0]

    # Orszag 2/3 截断反走样掩膜
    kmax_x = np.max(np.abs(KX))
    kmax_y = np.max(np.abs(KY))
    dealias = np.ones((N1, N2), dtype=bool)
    dealias[np.abs(KX) > (2.0/3.0) * kmax_x] = False
    dealias[np.abs(KY) > (2.0/3.0) * kmax_y] = False

    # 扩散项解析精确因子
    diff_factor = np.exp(-nu * K_sq * dt)

    # 严格按照讲义设定初始条件
    if model_case == 1:
        # 讲义式 6.56
        zeta0 = -2.0 * np.cos(X) * np.cos(Y)

    elif model_case == 2:
        # 讲义式 6.60
        r1_sq = (X - xv1)**2 + (Y - yv1)**2
        r2_sq = (X - xv2)**2 + (Y - yv2)**2
        zeta0 = (Av * np.exp(-r1_sq / Lv**2) +
                 Av * np.exp(-r2_sq / Lv**2) -
                 2.0 * np.pi * Av * (Lv**2) / (Lx * Ly))

    elif model_case == 3:
        # 讲义式 6.61
        zeta0 = Ab * np.exp(-Y**2 / Lb**2) + xi0 * np.random.randn(N1, N2)

    elif model_case == 4:
        # 讲义式 6.63
        r_sq = (X - xv1)**2 + (Y - yv1)**2
        zeta0 = Av * np.exp(-r_sq / Lv**2)

    elif model_case == 5:
        k_shift = np.fft.fftshift(K_array)
        amp = np.exp(-(k_shift - kc)**2 / (dk**2))
        shift_spec = np.zeros((N1, N2), dtype=complex)
        
        # 下半平面随机相位
        for i in range(1, N1):
            for j in range(1, N2 // 2 + 1):
                phase = np.exp(1j * 2.0 * np.pi * np.random.rand())
                shift_spec[i, j] = amp[i, j] * phase

        # ky = 0 交界区
        j_zero = N2 // 2
        for i in range(N1 // 2 + 1, N1):
            shift_spec[i, j_zero] = np.conj(shift_spec[N1 - i, j_zero])
        shift_spec[N1 // 2, j_zero] = 0.0

        # 上半平面共轭镜像
        for i in range(1, N1):
            for j in range(N2 // 2 + 1, N2):
                ii = N1 - i
                jj = N2 - j
                if 0 <= ii < N1 and 0 <= jj < N2:
                    shift_spec[i, j] = np.conj(shift_spec[ii, jj])

        zeta_hat = np.fft.fftshift(shift_spec)
        zeta0 = np.real(np.fft.ifft2(zeta_hat))
        # 讲义表 5: std(zeta) = 0.25
        zeta0 = zeta0 / np.std(zeta0) * 0.25

    zeta_hat = np.fft.fft2(zeta0)

    def calc_rhs(z_hat):
        psi_h = z_hat * inv_Ksq
        u_h  = -1j * KY * psi_h * dealias
        v_h  =  1j * KX * psi_h * dealias
        zx_h =  1j * KX * z_hat * dealias
        zy_h =  1j * KY * z_hat * dealias

        u  = np.real(np.fft.ifft2(u_h))
        v  = np.real(np.fft.ifft2(v_h))
        zx = np.real(np.fft.ifft2(zx_h))
        zy = np.real(np.fft.ifft2(zy_h))

        adv = u * zx + v * zy
        adv_h = np.fft.fft2(adv) * dealias

        if beta_star != 0.0:
            v_lin_h = 1j * KX * psi_h
            return -adv_h - beta_star * v_lin_h
        else:
            return -adv_h

    if model_case == 1:
        psi_h = zeta_hat * inv_Ksq
        u = np.real(np.fft.ifft2(-1j * KY * psi_h))
        v = np.real(np.fft.ifft2( 1j * KX * psi_h))
        KE0 = 0.5 * np.mean(u**2 + v**2)
        ke_times = [0.0]
        ke_num = [KE0]
        ke_exact = [KE0]

    print('=' * 60)
    print(f'开始模拟: {case_name}')
    print(f'网格: {N1}x{N2}, Re: {Re}, 总步数: {total_steps}, dt: {dt:.2e}')
    print('=' * 60)

    t_start = time.time()

    for step in range(1, total_steps + 1):
        g1 = calc_rhs(zeta_hat)
        g2 = calc_rhs(zeta_hat + 0.5 * dt * g1)
        g3 = calc_rhs(zeta_hat + 0.5 * dt * g2)
        g4 = calc_rhs(zeta_hat + dt * g3)

        zeta_hat = zeta_hat + (dt / 6.0) * (g1 + 2.0 * g2 + 2.0 * g3 + g4)
        zeta_hat = zeta_hat * diff_factor

        if model_case == 1:
            t_now = step * dt
            psi_h = zeta_hat * inv_Ksq
            u = np.real(np.fft.ifft2(-1j * KY * psi_h))
            v = np.real(np.fft.ifft2( 1j * KX * psi_h))
            ke_times.append(t_now)
            ke_num.append(0.5 * np.mean(u**2 + v**2))
            ke_exact.append(KE0 * np.exp(-4.0 * nu * t_now))

        if step % (plot_interval * 4) == 0 or step == total_steps:
            pct = step / total_steps * 100.0
            print(f'进度: {pct:5.1f}% | 步数: {step}/{total_steps} | 当前时间: {step*dt:.3f}')

    t_end = time.time()
    print(f'计算完成！耗时: {t_end - t_start:.2f} 秒')

    return {'x': x, 'y': y, 'vorticity': np.real(np.fft.ifft2(zeta_hat)),
            'time': T, 'Re': Re, 'beta': beta_star,
            **({'energy_times': np.array(ke_times), 'energy': np.array(ke_num),
                'exact_energy': np.array(ke_exact)} if model_case == 1 else {})}

if __name__ == '__main__':
    from pathlib import Path
    import argparse
    parser = argparse.ArgumentParser(description='Run a barotropic vorticity benchmark')
    parser.add_argument('--case', type=int, choices=range(1,6), default=1)
    parser.add_argument('--output', type=Path, default=Path('outputs/bve_result.npz'))
    args = parser.parse_args()
    result = run_simulation(args.case)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **result)
    print('Saved numerical result:', args.output)
