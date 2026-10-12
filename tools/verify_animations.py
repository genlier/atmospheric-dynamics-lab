"""Check frame counts, timing and media metadata for the published trajectories."""
from pathlib import Path
import json,subprocess,shutil,hashlib
import numpy as np
from render_dense_animations import load

ROOT=Path(__file__).resolve().parents[1];PUBLIC=ROOT/'docs'/'animations'
probe=shutil.which('ffprobe') or 'ffprobe'
labels=[f'case{i}' for i in range(1,6)]+['asymmetric','shallow-water','qg','thermal']
manifest=[]
for label in labels:
    times,values,meta,_=load(label)
    stored=json.loads((PUBLIC/(label+'.json')).read_text())
    assert len(times)==stored['frames']==601
    assert times[0]==0 and abs(times[-1]-stored['integrationSteps']*stored['dt'])<1e-10
    assert np.isfinite(values).all() and values.shape[-1]==stored['grid']
    for which,variant in enumerate(stored['variants']):
        path=PUBLIC/variant['video'];folder=ROOT/'outputs'/'dense'/'frames'/f'{label}-{which}'
        assert len(list(folder.glob('frame_*.png')))==601
        information=json.loads(subprocess.check_output([probe,'-v','error','-show_streams','-show_format','-of','json',str(path)]))
        stream=next(s for s in information['streams'] if s['codec_type']=='video')
        assert int(stream['nb_frames'])==601 and stream['avg_frame_rate']=='30/1'
        assert abs(float(information['format']['duration'])-601/30)<.001
        tags=information['format'].get('tags',{})
        assert set(tags)<= {'major_brand','minor_version','compatible_brands','encoder'}
        manifest.append({'case':label,'field':which,'file':variant['video'],'frames':601,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
times,fields,meta,_=load('case1');n=meta['grid'];x=np.linspace(-np.pi,np.pi,n,endpoint=False);X,Y=np.meshgrid(x,x,indexing='ij')
exact=-2*np.cos(X)*np.cos(Y)*np.exp(-2*times[-1]/meta['Re'])
assert np.max(np.abs(fields[-1,0]-exact))<1e-12
(PUBLIC/'manifest.json').write_text(json.dumps({'sourceSHA256':hashlib.sha256((ROOT/'models/matlab/bve.m').read_bytes()).hexdigest(),'videos':manifest},indent=2),encoding='utf-8')
print(f'Verified {len(manifest)} videos: 601 actual PNG frames each, full integration duration, finite fields, media metadata and analytic vortex decay.')
