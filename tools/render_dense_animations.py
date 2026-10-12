"""Save actual numerical checkpoints as PNG frames and encode them in order."""
from pathlib import Path
import argparse,json,subprocess,shutil
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'dense'
PUBLIC=ROOT/'docs'/'animations';PUBLIC.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],'font.size':12})
NAMES=['Decaying vortex grid','Vortex interaction','Shear instability','Vortex drift','Two-dimensional turbulence']
FIELDS={'asymmetric':['Relative vorticity'],'shallow-water':['Exp. 1: vorticity','Exp. 1: total height','Exp. 2: vorticity','Exp. 2: total height'],'qg':['Upper-layer streamfunction','Streamfunction difference','Upper-layer PV perturbation'],'thermal':['Temperature perturbation','Relative vorticity']}

def load(label):
    if label.startswith('case'):
        a=loadmat(OUT/(label+'.mat'),simplify_cells=True);p=a['p'];values=a['animation_snapshots'].transpose(2,0,1)[:,None,:,:]
        return a['animation_times'],values,{'grid':int(p['N1']),'dt':float(p['dt']),'integrationSteps':int(p['Step']),'Re':float(p['Re']),'beta':float(p['beta_star']),'name':NAMES[int(label[-1])-1],'fields':['Relative vorticity']},a
    if label=='shallow-water':
        a=np.load(OUT/'swe1.npz');b=np.load(OUT/'swe2.npz');assert np.array_equal(a['times'],b['times']);times=a['times'];values=np.concatenate([a['values'],b['values']],axis=1)
    else:
        a=np.load(OUT/(label+'.npz'));times=a['times'];values=a['values']
    dt={'asymmetric':.0001,'shallow-water':.002,'qg':.005,'thermal':.002}[label]
    return times,values,{'grid':64,'dt':dt,'integrationSteps':int(round(times[-1]/dt)),'fields':FIELDS[label]},None

def render(label,ffmpeg):
    times,values,meta,source=load(label);assert np.isfinite(values).all();count=len(times);assert count==601
    meta.update(frames=count,fps=30,firstTime=float(times[0]),lastTime=float(times[-1]),variants=[])
    n=values.shape[-1];x=np.linspace(-np.pi,np.pi,n,endpoint=False)
    for which,field in enumerate(meta['fields']):
        folder=OUT/'frames'/f'{label}-{which}';folder.mkdir(parents=True,exist_ok=True)
        diagnostic=label in ['case1','case5']
        fig,axes=plt.subplots(1,2 if diagnostic else 1,figsize=(10,5) if diagnostic else (6,6),dpi=100)
        ax=axes[0] if diagnostic else axes
        image=ax.imshow(values[0,which].T,origin='lower',extent=[x[0],x[-1],x[0],x[-1]],cmap='jet',interpolation='bilinear',aspect='equal')
        ax.set_xlabel('x');ax.set_ylabel('z' if label=='thermal' else 'y');title=ax.set_title('');fig.colorbar(image,ax=ax,fraction=.046,pad=.04)
        if diagnostic:
            other=axes[1];other.set_box_aspect(1);other.grid(True)
            if label=='case1':
                energies=np.sum(source['animation_spectra'],axis=0)
                exact=energies[0]*np.exp(-4*times/meta['Re'])
                other.plot(times,exact,'r--',label='Exact',lw=1.8);line,=other.plot([],[],'b-',label='Numerical',lw=2)
                other.set_xlim(0,times[-1]);other.set_ylim(0,.26);other.set_xlabel('Time t');other.set_ylabel('Kinetic energy');other.set_title('Energy decay');other.legend()
            else:
                other.set_xscale('log');other.set_yscale('log');line,=other.plot([],[],'b-',lw=2,label='E(k)');ref,=other.plot([],[],'r--',lw=1.8,label='k$^{-3}$ reference');other.set_xlim(1,source['animation_k'][-1]);other.set_ylim(1e-12,1e-2);other.set_xlabel('Wavenumber k');other.set_ylabel('Energy E(k)');other.set_title('Kinetic energy spectrum');other.legend(loc='lower left')
        fig.tight_layout()
        for index,time in enumerate(times):
            image.set_data(values[index,which].T);lo=float(values[index,which].min());hi=float(values[index,which].max());image.set_clim(lo if lo<hi else lo-1,hi if lo<hi else hi+1)
            title.set_text(f'{field}\nt = {time:.5g} · frame {index+1}/{count}')
            if label=='case1':line.set_data(times[:index+1],energies[:index+1])
            if label=='case5':
                k=source['animation_k'];spectrum=source['animation_spectra'][:,index];line.set_data(k,spectrum);m=(k>=4)&(k<=20);ref.set_data(k[m],spectrum[np.flatnonzero(k>=4)[0]]*(k[m]/k[m][0])**-3)
            fig.canvas.draw();Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3]).save(folder/f'frame_{index:04d}.png')
        plt.close(fig)
        name=f'{label}-{which}.mp4';video=PUBLIC/name
        subprocess.run([ffmpeg,'-y','-loglevel','error','-framerate','30','-i',str(folder/'frame_%04d.png'),'-frames:v',str(count),'-c:v','libx264','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart','-map_metadata','-1',str(video)],check=True)
        meta['variants'].append({'field':field,'video':name})
        print(f'Encoded {label}, field {which}: {count} actual PNG frames.',flush=True)
    (PUBLIC/(label+'.json')).write_text(json.dumps(meta,indent=2),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('labels',nargs='+');parser.add_argument('--ffmpeg',default=shutil.which('ffmpeg') or 'ffmpeg');args=parser.parse_args()
    for label in args.labels:render(label,args.ffmpeg)
