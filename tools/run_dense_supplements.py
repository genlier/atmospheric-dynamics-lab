"""Record dense checkpoints while executing the existing numerical solvers."""
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'dense';OUT.mkdir(parents=True,exist_ok=True)
remaining=(ROOT/'models/python/experiments.py').read_text(encoding='utf-8-sig')
bonus=(ROOT/'models/python/supplementary.py').read_text(encoding='utf-8-sig')
records={}; env={'np':np}

def keep(step,total):
    return step in targets[total]
targets={total:set(np.rint(np.linspace(0,total,601)).astype(int)) for total in [60000,12500,2400,1500]}

def record(name,step,dt,fields):
    assert all(np.isfinite(a).all() for a in fields)
    records.setdefault(name,[]).append((step*dt,np.stack(fields)))

env.update(keep=keep,record=record)
def run(source):exec(compile(source,'<numerical checkpoint run>','exec'),env)
def replace(source,a,b):
    assert source.count(a)==1,a
    return source.replace(a,b)

prefix=remaining[remaining.index('# Common Grid'):remaining.index('times2 =')]
prefix=replace(prefix,'    snapshots = {}','    record("asymmetric",0,dt,[zeta0])\n    snapshots = {}')
prefix=replace(prefix,'        zeta_hat = zeta_hat * diff_factor','        zeta_hat = zeta_hat * diff_factor\n        if keep(s,steps): record("asymmetric",s,dt,[np.real(np.fft.ifft2(zeta_hat))])')
run(prefix)
run('times2=[0.0,1.5,3.0,6.0]\n')
start=remaining.index('zeta0_c6 =');stop=remaining.index('fig, axs =',start)
print('Running asymmetric vortex: 60000 steps.',flush=True);run(remaining[start:stop])

start=remaining.index('def run_swe(');stop=remaining.index('fig, axs =',start)
source=remaining[start:stop]
source=replace(source,"    snap_t0 = {'h': h0_p, 'zeta': z0_p}","    snap_t0 = {'h': h0_p, 'zeta': z0_p}\n    record('swe'+str(exp_case),0,dt_p,[z0_p,1+h0_p])")
source=replace(source,'        zh_n, dh_n, hh_n = zh_np1, dh_np1, hh_np1','        zh_n, dh_n, hh_n = zh_np1, dh_np1, hh_np1\n        if keep(s,steps_p): record("swe"+str(exp_case),s,dt_p,[np.real(np.fft.ifft2(zh_n)),1+np.real(np.fft.ifft2(hh_n))])')
print('Running both shallow-water experiments: 12500 steps each.',flush=True);run(source)

start=bonus.index('N1, N2 =');stop=bonus.index('fig, axs =',start)
source=bonus[start:stop]
def qg_record(step):
    p1,p2=env['invert_qg'](env['q1_hat'],env['q2_hat'])
    record('qg',step,env['dt_qg'],[np.fft.ifft2(p1).real,np.fft.ifft2(p1-p2).real,np.fft.ifft2(env['q1_hat']).real])
env['qg_record']=qg_record
source=replace(source,'for s in range(1, steps_qg + 1):','qg_record(0)\nfor s in range(1, steps_qg + 1):')
source=replace(source,'    q2_hat *= diff_factor_qg','    q2_hat *= diff_factor_qg\n    if keep(s,steps_qg): qg_record(s)')
print('Running two-layer QG: 2400 steps.',flush=True);run(source)

start=bonus.index('nu_tc =');stop=bonus.index('fig, axs =',start)
source=bonus[start:stop]
source=replace(source,'for s in range(1, steps_tc + 1):','record("thermal",0,dt_tc,[theta0,zeta0])\nfor s in range(1, steps_tc + 1):')
source=replace(source,'    theta_h *= diff_theta','    theta_h *= diff_theta\n    if keep(s,steps_tc): record("thermal",s,dt_tc,[np.fft.ifft2(theta_h).real,np.fft.ifft2(zeta_h).real])')
print('Running thermal: 1500 steps.',flush=True);run(source)

for name,frames in records.items():
    times=np.array([f[0] for f in frames]);values=np.stack([f[1] for f in frames])
    assert len(times)==601 and np.all(np.diff(times)>0)
    np.savez_compressed(OUT/(name+'.npz'),times=times,values=values)
print('Recorded 601 actual frames for each supplementary run.',flush=True)
