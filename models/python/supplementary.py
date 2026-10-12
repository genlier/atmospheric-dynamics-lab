import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

from pathlib import Path
out_dir = str(Path(__file__).resolve().parents[2] / 'outputs')
os.makedirs(out_dir, exist_ok=True)

# -------------------------------------------------------------
# 1. Two-layer Quasi-Geostrophic (QG) Model: Baroclinic Instability
# -------------------------------------------------------------
print("Simulating Two-layer QG Model (Baroclinic Instability)...")
N1, N2 = 64, 64
Lx, Ly = 2 * np.pi, 2 * np.pi
x = np.linspace(-Lx/2, Lx/2, N1, endpoint=False)
y = np.linspace(-Ly/2, Ly/2, N2, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')

kx = np.fft.fftfreq(N1, d=Lx/(2*np.pi*N1))
ky = np.fft.fftfreq(N2, d=Ly/(2*np.pi*N2))
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K_sq = KX**2 + KY**2

# Orszag 2/3 dealiasing mask
kmax_x = np.max(np.abs(KX))
kmax_y = np.max(np.abs(KY))
dealias = (np.abs(KX) <= (2.0/3.0)*kmax_x) & (np.abs(KY) <= (2.0/3.0)*kmax_y)

# Parameters for Two-layer Phillips QG model
# q1 = \nabla^2 psi1 + F (psi2 - psi1) + beta y
# q2 = \nabla^2 psi2 + F (psi1 - psi2) + beta y
# psi+ = psi1 + psi2 -> q+ = \nabla^2 psi+
# psi- = psi1 - psi2 -> q- = (\nabla^2 - 2F) psi-
F = 4.0        # F = f0^2 / (g' H)
beta_qg = 1.0  # Planetary vorticity gradient
nu_qg = 1e-4   # Viscosity
r_drag = 0.05  # Ekman bottom drag on lower layer

inv_Ksq_plus = np.zeros_like(K_sq)
inv_Ksq_plus[K_sq > 0] = -1.0 / K_sq[K_sq > 0]

inv_Ksq_minus = -1.0 / (K_sq + 2.0 * F)

# Base state: Zonal shear flow U1 = 1.0, U2 = -1.0
# Perturbation on PV
q1_pert = 0.2 * np.sin(2 * X) * np.sin(Y) + 0.05 * np.random.randn(N1, N2)
q2_pert = -0.2 * np.sin(2 * X) * np.sin(Y) + 0.05 * np.random.randn(N1, N2)

# Initial PV fields
q1 = -F * ( -2.0 * Y ) + q1_pert  # background + perturbation
q2 =  F * ( -2.0 * Y ) + q2_pert

q1_hat = np.fft.fft2(q1_pert)
q2_hat = np.fft.fft2(q2_pert)

def invert_qg(q1_h, q2_h):
    qp_h = q1_h + q2_h
    qm_h = q1_h - q2_h
    psip_h = qp_h * inv_Ksq_plus
    psim_h = qm_h * inv_Ksq_minus
    psi1_h = 0.5 * (psip_h + psim_h)
    psi2_h = 0.5 * (psip_h - psim_h)
    return psi1_h, psi2_h

def rhs_qg(q1_h, q2_h):
    psi1_h, psi2_h = invert_qg(q1_h, q2_h)
    
    # Velocities: U1_mean = 0.5, U2_mean = -0.5
    u1_h = -1j * KY * psi1_h * dealias
    v1_h =  1j * KX * psi1_h * dealias
    u2_h = -1j * KY * psi2_h * dealias
    v2_h =  1j * KX * psi2_h * dealias
    
    q1x_h = 1j * KX * q1_h * dealias
    q1y_h = 1j * KY * q1_h * dealias
    q2x_h = 1j * KX * q2_h * dealias
    q2y_h = 1j * KY * q2_h * dealias
    
    u1 = np.real(np.fft.ifft2(u1_h)) + 0.8
    v1 = np.real(np.fft.ifft2(v1_h))
    u2 = np.real(np.fft.ifft2(u2_h)) - 0.8
    v2 = np.real(np.fft.ifft2(v2_h))
    
    q1x = np.real(np.fft.ifft2(q1x_h))
    q1y = np.real(np.fft.ifft2(q1y_h)) + (beta_qg + 2*F*0.8)
    q2x = np.real(np.fft.ifft2(q2x_h))
    q2y = np.real(np.fft.ifft2(q2y_h)) + (beta_qg - 2*F*0.8)
    
    adv1 = u1 * q1x + v1 * q1y
    adv2 = u2 * q2x + v2 * q2y
    
    adv1_h = np.fft.fft2(adv1) * dealias
    adv2_h = np.fft.fft2(adv2) * dealias
    
    # Bottom drag on layer 2
    drag2_h = r_drag * (-K_sq) * psi2_h
    
    return -adv1_h, -adv2_h - drag2_h

# Time stepping RK4
dt_qg = 0.005
steps_qg = 2400  # T = 12.0
diff_factor_qg = np.exp(-nu_qg * K_sq * dt_qg)

for s in range(1, steps_qg + 1):
    k1_1, k1_2 = rhs_qg(q1_hat, q2_hat)
    k2_1, k2_2 = rhs_qg(q1_hat + 0.5*dt_qg*k1_1, q2_hat + 0.5*dt_qg*k1_2)
    k3_1, k3_2 = rhs_qg(q1_hat + 0.5*dt_qg*k2_1, q2_hat + 0.5*dt_qg*k2_2)
    k4_1, k4_2 = rhs_qg(q1_hat + dt_qg*k3_1, q2_hat + dt_qg*k3_2)
    
    q1_hat += (dt_qg / 6.0) * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)
    q2_hat += (dt_qg / 6.0) * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)
    
    q1_hat *= diff_factor_qg
    q2_hat *= diff_factor_qg

psi1_h, psi2_h = invert_qg(q1_hat, q2_hat)
psi1 = np.real(np.fft.ifft2(psi1_h))
psi2 = np.real(np.fft.ifft2(psi2_h))
q1_final = np.real(np.fft.ifft2(q1_hat))
temp_interface = psi1 - psi2  # Interfacial temperature/thickness perturbation

fig, axs = plt.subplots(1, 3, figsize=(14, 4))
c1 = axs[0].contourf(X, Y, psi1, levels=30, cmap='RdBu_r')
fig.colorbar(c1, ax=axs[0])
axs[0].set_title('(a) 上层流函数 $\\psi_1$ (斜压波卷紧)', fontsize=11)
axs[0].set_xlabel('x'); axs[0].set_ylabel('y'); axs[0].set_aspect('equal')

c2 = axs[1].contourf(X, Y, temp_interface, levels=30, cmap='coolwarm')
fig.colorbar(c2, ax=axs[1])
axs[1].set_title('(b) 界面热力层厚度 $\\psi_1 - \\psi_2$ (热成风锋面)', fontsize=11)
axs[1].set_xlabel('x'); axs[1].set_ylabel('y'); axs[1].set_aspect('equal')

c3 = axs[2].contourf(X, Y, q1_final, levels=30, cmap='jet')
fig.colorbar(c3, ax=axs[2])
axs[2].set_title('(c) 上层位涡扰动 $q_1\'$ (中纬度气旋/反气旋)', fontsize=11)
axs[2].set_xlabel('x'); axs[2].set_ylabel('y'); axs[2].set_aspect('equal')

plt.suptitle('选做拓展：两层准地转模式斜压不稳定数值模拟 (Phillips 模型)', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig8_qg_baroclinic.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Fig 8 saved!")

# -------------------------------------------------------------
# 2. Thermal Convection (2D Boussinesq Thermal Plume / Rising Bubble)
# -------------------------------------------------------------
print("Simulating Thermal Convection (2D Boussinesq Plume)...")
# Equations:
# d(zeta)/dt + J(psi, zeta) = nu \nabla^2 zeta + g*alpha * d(theta)/dx
# d(theta)/dt + J(psi, theta) = kappa \nabla^2 theta
# \nabla^2 psi = zeta

nu_tc = 1e-3
kappa_tc = 1e-3
B_tc = 2.0   # g * alpha buoyancy parameter

# Initial condition: warm bubble at center-bottom (0, -1.2)
r_bubble = 0.4
theta0 = 2.0 * np.exp(-(X**2 + (Y + 1.2)**2) / r_bubble**2)
zeta0 = np.zeros((N1, N2))

zeta_h = np.fft.fft2(zeta0)
theta_h = np.fft.fft2(theta0)

inv_Ksq_tc = np.zeros_like(K_sq)
inv_Ksq_tc[K_sq > 0] = -1.0 / K_sq[K_sq > 0]

diff_zeta = np.exp(-nu_tc * K_sq * 0.002)
diff_theta = np.exp(-kappa_tc * K_sq * 0.002)
dt_tc = 0.002
steps_tc = 1500 # T = 3.0

def rhs_tc(zh, th):
    psi_h = zh * inv_Ksq_tc
    uh = -1j * KY * psi_h * dealias
    vh =  1j * KX * psi_h * dealias
    zx_h = 1j * KX * zh * dealias
    zy_h = 1j * KY * zh * dealias
    tx_h = 1j * KX * th * dealias
    ty_h = 1j * KY * th * dealias
    
    u = np.real(np.fft.ifft2(uh))
    v = np.real(np.fft.ifft2(vh))
    zx = np.real(np.fft.ifft2(zx_h))
    zy = np.real(np.fft.ifft2(zy_h))
    tx = np.real(np.fft.ifft2(tx_h))
    ty = np.real(np.fft.ifft2(ty_h))
    
    adv_z = u * zx + v * zy
    adv_t = u * tx + v * ty
    
    # Buoyancy baroclinic generation torque: B * d(theta)/dx
    buoy_h = B_tc * (1j * KX * th) * dealias
    
    return -np.fft.fft2(adv_z)*dealias + buoy_h, -np.fft.fft2(adv_t)*dealias

snapshots_th = [theta0.copy()]
snapshots_zt = [zeta0.copy()]
snap_times_tc = [0.0, 1.0, 2.0, 3.0]
current_time = 0.0

for s in range(1, steps_tc + 1):
    k1_z, k1_t = rhs_tc(zeta_h, theta_h)
    k2_z, k2_t = rhs_tc(zeta_h + 0.5*dt_tc*k1_z, theta_h + 0.5*dt_tc*k1_t)
    k3_z, k3_t = rhs_tc(zeta_h + 0.5*dt_tc*k2_z, theta_h + 0.5*dt_tc*k2_t)
    k4_z, k4_t = rhs_tc(zeta_h + dt_tc*k3_z, theta_h + dt_tc*k3_t)
    
    zeta_h += (dt_tc / 6.0) * (k1_z + 2*k2_z + 2*k3_z + k4_z)
    theta_h += (dt_tc / 6.0) * (k1_t + 2*k2_t + 2*k3_t + k4_t)
    
    zeta_h *= diff_zeta
    theta_h *= diff_theta
    current_time += dt_tc
    
    for st in snap_times_tc[1:]:
        if abs(current_time - st) < 0.51 * dt_tc and len(snapshots_th) < len(snap_times_tc):
            snapshots_th.append(np.real(np.fft.ifft2(theta_h)))
            snapshots_zt.append(np.real(np.fft.ifft2(zeta_h)))

fig, axs = plt.subplots(2, 4, figsize=(16, 7.5))
for i, st in enumerate(snap_times_tc):
    c1 = axs[0, i].contourf(X, Y, snapshots_th[i], levels=35, cmap='hot_r')
    fig.colorbar(c1, ax=axs[0, i], fraction=0.046, pad=0.04)
    axs[0, i].set_title(f'温度扰动 $\\theta$ (t={st})', fontsize=11)
    axs[0, i].set_xlabel('x'); axs[0, i].set_ylabel('y'); axs[0, i].set_aspect('equal')
    
    c2 = axs[1, i].contourf(X, Y, snapshots_zt[i], levels=35, cmap='bwr')
    fig.colorbar(c2, ax=axs[1, i], fraction=0.046, pad=0.04)
    axs[1, i].set_title(f'涡度场 $\\zeta$ (t={st})', fontsize=11)
    axs[1, i].set_xlabel('x'); axs[1, i].set_ylabel('y'); axs[1, i].set_aspect('equal')

plt.suptitle('选做拓展：二维Boussinesq热对流与热泡上升对流卷模拟 (浮力斜压扭矩生成对称涡偶极子)', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'fig9_thermal_convection.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Fig 9 saved! Both bonus figures completed.")
