"""Refresh the small initial-field previews from the dense numerical runs."""
from pathlib import Path
import json
from render_dense_animations import load

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'docs'/'data.js'
text=path.read_text(encoding='utf-8-sig').strip()
data=json.loads(text.removeprefix('window.LAB_DATA = ').removesuffix(';'))
for c in data['cases']:
    times,values,meta,source=load(c['id'])
    stride=max(1,meta['grid']//64)
    indices=[0,150,300,600]
    c['frames']=[{'time':float(times[i]),'fields':[[float(format(float(v),'.7g')) for v in field[::stride,::stride].T.ravel()] for field in values[i]]} for i in indices]
    if c['id']=='case5':c['spectrum']={'k':source['animation_k'].tolist(),'frames':source['animation_spectra'][:,indices].T.tolist()}
    c['source']='Dense full-trajectory numerical run'
data['note']='Full trajectories recorded as PNG sequences and encoded without temporal interpolation.'
path.write_text('window.LAB_DATA = '+json.dumps(data,ensure_ascii=False,separators=(',',':'))+';\n',encoding='utf-8')
print('Initial previews updated from all nine dense trajectories.')
