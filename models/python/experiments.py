import os, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

from pathlib import Path
out_dir = str(Path(__file__).resolve().parents[2] / 'outputs')
os.makedirs(out_dir, exist_ok=True)

# Common Grid for 64x64
N1, N2 = 64, 64
Lx, Ly = 2 * np.pi, 2 * np.pi
x = np.linspace(-Lx/2, Lx/2, N1, endpoint=False)
y = np.linspace(-Ly/2, Ly/2, N2, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')

kx = np.fft.fftfreq(N1, d=Lx/(2*np.pi*N1))
ky = np.fft.fftfreq(N2, d=Ly/(2*np.pi*N2))
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K_sq = KX**2 + KY**2
inv_Ksq = np.zeros_like(K_sq)
inv_Ksq[K_sq > 0] = -1.0 / K_sq[K_sq > 0]

kmax_x = np.max(np.abs(KX))
kmax_y = np.max(np.abs(KY))
dealias = np.ones((N1, N2), dtype=bool)
dealias[np.abs(KX) > (2.0/3.0)*kmax_x] = False
dealias[np.abs(KY) > (2.0/3.0)*kmax_y] = False

def solve_bve(zeta0, Re, T, dt, beta_star=0.0, save_times=[0.0, 1.5, 3.0, 6.0]):
    nu = 1.0 / Re
    steps = int(round(T / dt))
    zeta_hat = np.fft.fft2(zeta0)
    diff_factor = np.exp(-nu * K_sq * dt)
    
    def rhs(z_hat):
        psi_h = z_hat * inv_Ksq
        u_h = -1j * KY * psi_h * dealias
        v_h = 1j * KX * psi_h * dealias
        zx_h = 1j * KX * z_hat * dealias
        zy_h = 1j * KY * z_hat * dealias
        u = np.real(np.fft.ifft2(u_h))
        v = np.real(np.fft.ifft2(v_h))
        zx = np.real(np.fft.ifft2(zx_h))
        zy = np.real(np.fft.ifft2(zy_h))
        adv = u * zx + v * zy
        adv_h = np.fft.fft2(adv) * dealias
        if beta_star != 0.0:
            v_lin_h = 1j * KX * psi_h
            return -adv_h - beta_star * v_lin_h
        return -adv_h

    snapshots = {}
    if 0.0 in save_times:
        snapshots[0.0] = np.real(np.fft.ifft2(zeta_hat))
        
    for s in range(1, steps + 1):
        t_now = s * dt
        k1 = rhs(zeta_hat)
        k2 = rhs(zeta_hat + 0.5 * dt * k1)
        k3 = rhs(zeta_hat + 0.5 * dt * k2)
        k4 = rhs(zeta_hat + dt * k3)
        zeta_hat = zeta_hat + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        zeta_hat = zeta_hat * diff_factor
        
        for st in save_times:
            if abs(t_now - st) < 0.51 * dt and st not in snapshots:
                snapshots[st] = np.real(np.fft.ifft2(zeta_hat))
                
    return snapshots

times2 = [0.0, 1.5, 3.0, 6.0]

# --- 3. Case 3: Shear Instability ---
print("Generating Case 3...")
Ab, Lb, xi0 = 40.0, 0.1, 0.1
zeta0_c3 = Ab * np.exp(-Y**2 / Lb**2) + xi0 * np.random.randn(N1, N2)
snaps3 = solve_bve(zeta0_c3, Re=100.0, T=6.0, dt=1.0e-4, save_times=times2)

fig, axs = plt.subplots(1, 4, figsize=(16, 3.8))
for idx, t_val in enumerate(times2):
    c = axs[idx].contourf(X, Y, snaps3[t_val], levels=35, cmap='jet')
    fig.colorbar(c, ax=axs[idx], fraction=0.046, pad=0.04)
    axs[idx].set_title(f't = {t_val}', fontsize=12)
    axs[idx].set_xlabel('x'); axs[idx].set_ylabel('y'); axs[idx].set_aspect('equal')
plt.suptitle('测试3：剪切不稳定 (Shear Instability) 涡度场演变', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig3_case3_shear.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Fig 3 saved!")

# --- 4. Case 4: Beta Drift ---
print("Generating Case 4...")
Av, Lv = 40.0, 0.4
zeta0_c4 = Av * np.exp(-(X**2 + Y**2)/Lv**2)
snaps4 = solve_bve(zeta0_c4, Re=50.0, T=6.0, dt=1.0e-4, beta_star=1.0, save_times=times2)

fig, axs = plt.subplots(1, 4, figsize=(16, 3.8))
for idx, t_val in enumerate(times2):
    c = axs[idx].contourf(X, Y, snaps4[t_val], levels=35, cmap='jet')
    fig.colorbar(c, ax=axs[idx], fraction=0.046, pad=0.04)
    axs[idx].set_title(f't = {t_val}', fontsize=12)
    axs[idx].set_xlabel('x'); axs[idx].set_ylabel('y'); axs[idx].set_aspect('equal')
plt.suptitle('测试4：涡旋漂移 (Beta Drift) 涡度场演变及罗斯贝波频散', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig4_case4_beta_drift.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Fig 4 saved!")

# --- 5. Case 5: 2D Turbulence ---
print("Generating Case 5...")
N_turb = 128
x_t = np.linspace(-Lx/2, Lx/2, N_turb, endpoint=False)
y_t = np.linspace(-Ly/2, Ly/2, N_turb, endpoint=False)
Xt, Yt = np.meshgrid(x_t, y_t, indexing='ij')
kxt = np.fft.fftfreq(N_turb, d=Lx/(2*np.pi*N_turb))
kyt = np.fft.fftfreq(N_turb, d=Ly/(2*np.pi*N_turb))
KXt, KYt = np.meshgrid(kxt, kyt, indexing='ij')
Ksq_t = KXt**2 + KYt**2
K_arr_t = np.sqrt(Ksq_t)
inv_Ksq_t = np.zeros_like(Ksq_t)
inv_Ksq_t[Ksq_t > 0] = -1.0 / Ksq_t[Ksq_t > 0]

kmax_xt = np.max(np.abs(KXt))
kmax_yt = np.max(np.abs(KYt))
dealias_t = np.ones((N_turb, N_turb), dtype=bool)
dealias_t[np.abs(KXt) > (2.0/3.0)*kmax_xt] = False
dealias_t[np.abs(KYt) > (2.0/3.0)*kmax_yt] = False

# Conjugate symmetric spectrum
kc_t = 2 * np.pi / 0.5
dk_t = 2 * np.pi / 2.0
k_shift_t = np.fft.fftshift(K_arr_t)
amp_t = np.exp(-(k_shift_t - kc_t)**2 / dk_t**2)
shift_spec_t = np.zeros((N_turb, N_turb), dtype=complex)
for i in range(1, N_turb):
    for j in range(1, N_turb//2 + 1):
        shift_spec_t[i, j] = amp_t[i, j] * np.exp(1j * 2 * np.pi * np.random.rand())
j_zero = N_turb // 2
for i in range(N_turb//2 + 1, N_turb):
    shift_spec_t[i, j_zero] = np.conj(shift_spec_t[N_turb - i, j_zero])
shift_spec_t[N_turb//2, j_zero] = 0.0
for i in range(1, N_turb):
    for j in range(N_turb//2 + 1, N_turb):
        ii = N_turb - i; jj = N_turb - j
        if 0 <= ii < N_turb and 0 <= jj < N_turb:
            shift_spec_t[i, j] = np.conj(shift_spec_t[ii, jj])
zh_turb = np.fft.fftshift(shift_spec_t)
z0_turb = np.real(np.fft.ifft2(zh_turb))
z0_turb = z0_turb / np.std(z0_turb) * 0.25
zh_turb = np.fft.fft2(z0_turb)

nu_turb = 1.0 / 4.0e4
dt_turb = 0.02
diff_factor_t = np.exp(-nu_turb * Ksq_t * dt_turb)

def rhs_turb(z_h):
    psi_h = z_h * inv_Ksq_t
    u_h = -1j * KYt * psi_h * dealias_t
    v_h =  1j * KXt * psi_h * dealias_t
    zx_h = 1j * KXt * z_h * dealias_t
    zy_h = 1j * KYt * z_h * dealias_t
    u = np.real(np.fft.ifft2(u_h))
    v = np.real(np.fft.ifft2(v_h))
    zx = np.real(np.fft.ifft2(zx_h))
    zy = np.real(np.fft.ifft2(zy_h))
    adv = u * zx + v * zy
    return -np.fft.fft2(adv) * dealias_t

for s in range(1500):
    k1 = rhs_turb(zh_turb)
    k2 = rhs_turb(zh_turb + 0.5 * dt_turb * k1)
    k3 = rhs_turb(zh_turb + 0.5 * dt_turb * k2)
    k4 = rhs_turb(zh_turb + dt_turb * k3)
    zh_turb = zh_turb + (dt_turb / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    zh_turb = zh_turb * diff_factor_t

zt_final = np.real(np.fft.ifft2(zh_turb))

psi_ht = zh_turb * inv_Ksq_t
uh_t = -1j * KYt * psi_ht
vh_t =  1j * KXt * psi_ht
E_2D = 0.5 * (np.abs(uh_t)**2 + np.abs(vh_t)**2)
max_k = int(N_turb / 3)
k_bins = np.arange(1, max_k + 1)
Ek = np.zeros(max_k)
for i, kval in enumerate(k_bins):
    mask = (K_arr_t >= (kval - 0.5)) & (K_arr_t < (kval + 0.5))
    Ek[i] = np.sum(E_2D[mask])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
c = ax1.contourf(Xt, Yt, zt_final, levels=40, cmap='jet')
fig.colorbar(c, ax=ax1)
ax1.set_title('(a) 二维衰减湍流涡度场 (t = 30)', fontsize=12)
ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_aspect('equal')

ax2.loglog(k_bins, Ek, 'b-', lw=2, label='数值模拟能谱 E(k)')
ref_k = k_bins[(k_bins >= 4) & (k_bins <= 20)]
if len(ref_k) > 0:
    ref_Ek = Ek[np.where(k_bins >= 4)[0][0]] * (ref_k / ref_k[0])**(-3.0)
    ax2.loglog(ref_k, ref_Ek, 'r--', lw=1.8, label=r'理论拟涡能级联律 $k^{-3}$')
ax2.set_xlabel('波数 k', fontsize=11); ax2.set_ylabel('能量谱密度 E(k)', fontsize=11)
ax2.set_title(r'(b) 动能一维谱与 $k^{-3}$ 标度律对比', fontsize=12)
ax2.legend(fontsize=10); ax2.grid(True, which='both', linestyle=':')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig5_case5_turbulence.png'), dpi=200)
plt.close()
print("Fig 5 saved!")

# --- 6. New Exploration: Asymmetric Vortex Merger / Fujiwhara Effect ---
print("Generating Case 6...")
zeta0_c6 = (60.0 * np.exp(-((X - 0.5)**2 + Y**2)/0.5**2) +
            25.0 * np.exp(-((X + 0.8)**2 + Y**2)/0.3**2))
zeta0_c6 -= np.mean(zeta0_c6)
snaps6 = solve_bve(zeta0_c6, Re=80.0, T=6.0, dt=1.0e-4, save_times=times2)

fig, axs = plt.subplots(1, 4, figsize=(16, 3.8))
for idx, t_val in enumerate(times2):
    c = axs[idx].contourf(X, Y, snaps6[t_val], levels=35, cmap='jet')
    fig.colorbar(c, ax=axs[idx], fraction=0.046, pad=0.04)
    axs[idx].set_title(f't = {t_val}', fontsize=12)
    axs[idx].set_xlabel('x'); axs[idx].set_ylabel('y'); axs[idx].set_aspect('equal')
plt.suptitle('拓展探索：非对称双涡旋相互作用（藤原效应与强涡撕裂弱涡）', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig6_case_new_exploration.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Fig 6 saved!")

# --- 7. SWE Geostrophic Adjustment (Exp 1 vs Exp 2) ---
print("Generating Case 7...")
def run_swe(exp_case):
    f_p, H_p, g_p, Lv_p = 1.0, 1.0, 1.0, 0.2
    alpha_sp_p, r_sp_p, d_sp_p = 40.0, 2.5, 1.5
    gamma_p, dt_p, T_p = 0.06, 0.002, 25.0
    steps_p = int(round(T_p / dt_p))
    
    alpha_p = np.zeros((N1, N2))
    R_p = np.sqrt(X**2 + Y**2)
    mask_sp = R_p > r_sp_p
    alpha_p[mask_sp] = alpha_sp_p * np.tanh((R_p[mask_sp] - r_sp_p) / d_sp_p)
    
    if exp_case == 1:
        Ah_p = -0.3
        h0_p = Ah_p * np.exp(-(X**2 + Y**2)/Lv_p**2) - np.pi*Ah_p*Lv_p**2/(Lx*Ly)
        z0_p = np.zeros((N1, N2))
    else:
        Av_p = 0.2
        z0_p = Av_p * np.exp(-(X**2 + Y**2)/Lv_p**2) - np.pi*Av_p*Lv_p**2/(Lx*Ly)
        h0_p = np.zeros((N1, N2))
        
    d0_p = np.zeros((N1, N2))
    zh_n = np.fft.fft2(z0_p); dh_n = np.fft.fft2(d0_p); hh_n = np.fft.fft2(h0_p)
    zh_nm1 = zh_n.copy(); dh_nm1 = dh_n.copy(); hh_nm1 = hh_n.copy()
    
    inv_Ksq_swe = np.zeros_like(K_sq)
    inv_Ksq_swe[K_sq > 0] = 1.0 / K_sq[K_sq > 0]
    denom = 1.0 + g_p * H_p * K_sq * (dt_p**2)
    
    def compute_P_swe(zh, dh, hh):
        uh = (-1j * KX * inv_Ksq_swe * dh + 1j * KY * inv_Ksq_swe * zh) * dealias
        vh = (-1j * KY * inv_Ksq_swe * dh - 1j * KX * inv_Ksq_swe * zh) * dealias
        zx_h = 1j * KX * zh * dealias; zy_h = 1j * KY * zh * dealias
        dx_h = 1j * KX * dh * dealias; dy_h = 1j * KY * dh * dealias
        hx_h = 1j * KX * hh * dealias; hy_h = 1j * KY * hh * dealias
        u = np.real(np.fft.ifft2(uh)); v = np.real(np.fft.ifft2(vh))
        z = np.real(np.fft.ifft2(zh * dealias)); d = np.real(np.fft.ifft2(dh * dealias)); h = np.real(np.fft.ifft2(hh * dealias))
        zx = np.real(np.fft.ifft2(zx_h)); zy = np.real(np.fft.ifft2(zy_h))
        dx = np.real(np.fft.ifft2(dx_h)); dy = np.real(np.fft.ifft2(dy_h))
        hx = np.real(np.fft.ifft2(hx_h)); hy = np.real(np.fft.ifft2(hy_h))
        Pz = -u*zx - v*zy - (f_p + z)*d - alpha_p*z
        Pd = -u*dx - v*dy + (f_p + z)*z - alpha_p*d
        Ph = -h*d - u*hx - v*hy - alpha_p*h
        return np.fft.fft2(Pz)*dealias, np.fft.fft2(Pd)*dealias, np.fft.fft2(Ph)*dealias

    snap_t0 = {'h': h0_p, 'zeta': z0_p}
    for s in range(1, steps_p + 1):
        Pz, Pd, Ph = compute_P_swe(zh_n, dh_n, hh_n)
        if s == 1:
            zh_np1 = zh_n + dt_p * Pz
            Rd = dh_n + dt_p * Pd + g_p * K_sq * dt_p * hh_n
            Rh = hh_n + dt_p * Ph - H_p * dt_p * dh_n
            dh_np1 = (Rd + g_p * K_sq * dt_p * Rh) / denom
            hh_np1 = (Rh - H_p * dt_p * Rd) / denom
        else:
            zh_np1 = zh_nm1 + 2.0 * dt_p * Pz
            Rd = dh_nm1 + 2.0 * dt_p * Pd + 2.0 * g_p * K_sq * dt_p * hh_n
            Rh = hh_nm1 + 2.0 * dt_p * Ph - 2.0 * H_p * dt_p * dh_n
            dh_np1 = (Rd + 2.0 * g_p * K_sq * dt_p * Rh) / denom
            hh_np1 = (Rh - 2.0 * H_p * dt_p * Rd) / denom
            
            # Asselin filter
            zh_n = zh_n + gamma_p * (zh_np1 - 2.0*zh_n + zh_nm1)
            dh_n = dh_n + gamma_p * (dh_np1 - 2.0*dh_n + dh_nm1)
            hh_n = hh_n + gamma_p * (hh_np1 - 2.0*hh_n + hh_nm1)
            
        zh_nm1, dh_nm1, hh_nm1 = zh_n, dh_n, hh_n
        zh_n, dh_n, hh_n = zh_np1, dh_np1, hh_np1

    snap_tend = {'h': np.real(np.fft.ifft2(hh_n)), 'zeta': np.real(np.fft.ifft2(zh_n))}
    return snap_t0, snap_tend

snap1_0, snap1_end = run_swe(1)
snap2_0, snap2_end = run_swe(2)

fig, axs = plt.subplots(2, 2, figsize=(10, 8.5))
# Exp 1: Height initial and final
c1 = axs[0, 0].contourf(X, Y, snap1_0['h'], levels=30, cmap='RdBu_r')
fig.colorbar(c1, ax=axs[0, 0])
axs[0, 0].set_title('实验1 (初始高度扰动): 初始高度场 h (t = 0)', fontsize=11)
axs[0, 0].set_xlabel('x'); axs[0, 0].set_ylabel('y'); axs[0, 0].set_aspect('equal')

c2 = axs[0, 1].contourf(X, Y, snap1_end['h'], levels=30, cmap='RdBu_r')
fig.colorbar(c2, ax=axs[0, 1])
axs[0, 1].set_title('实验1 (小尺度重力波频散辐射): 最终高度场 h (t = 25)', fontsize=11)
axs[0, 1].set_xlabel('x'); axs[0, 1].set_ylabel('y'); axs[0, 1].set_aspect('equal')

# Exp 2: Vorticity initial and final height
c3 = axs[1, 0].contourf(X, Y, snap2_0['zeta'], levels=30, cmap='jet')
fig.colorbar(c3, ax=axs[1, 0])
axs[1, 0].set_title('实验2 (初始涡度扰动): 初始涡度场 $\\zeta$ (t = 0)', fontsize=11)
axs[1, 0].set_xlabel('x'); axs[1, 0].set_ylabel('y'); axs[1, 0].set_aspect('equal')

c4 = axs[1, 1].contourf(X, Y, snap2_end['h'], levels=30, cmap='RdBu_r')
fig.colorbar(c4, ax=axs[1, 1])
axs[1, 1].set_title('实验2 (地转适应形成稳定低压): 适应后高度场 h (t = 25)', fontsize=11)
axs[1, 1].set_xlabel('x'); axs[1, 1].set_ylabel('y'); axs[1, 1].set_aspect('equal')

plt.suptitle('选做拓展：旋转浅水模式地转适应过程模拟 (Exp 1 重力波频散 vs Exp 2 地转涡旋建立)', fontsize=13, y=0.99)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig7_swe_geostrophic_adj.png'), dpi=200)
plt.close()
print("Fig 7 saved! All figures completed successfully.")
