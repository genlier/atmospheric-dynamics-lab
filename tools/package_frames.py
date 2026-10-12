"""Package computed PNG sequences with their original integration times."""
from pathlib import Path
import json,zipfile,argparse
from render_dense_animations import load

ROOT=Path(__file__).resolve().parents[1];PUBLIC=ROOT/'docs'/'animations'
OUT=ROOT/'outputs'/'frame_archives';OUT.mkdir(parents=True,exist_ok=True)
parser=argparse.ArgumentParser();parser.add_argument('labels',nargs='+');args=parser.parse_args()
for label in args.labels:
    times,_,meta,_=load(label)
    info=PUBLIC/(label+'.json');manifest=json.loads(info.read_text())
    with zipfile.ZipFile(OUT/(label+'-png.zip'),'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as archive:
        archive.write(info,arcname=label+'.json')
        archive.writestr('times.csv','frame,time,step\n'+'\n'.join(f'{i:04d},{time:.17g},{round(time/meta["dt"])}' for i,time in enumerate(times)))
        for field in range(len(manifest['variants'])):
            folder=ROOT/'outputs'/'dense'/'frames'/f'{label}-{field}'
            images=sorted(folder.glob('frame_*.png'));assert len(images)==601
            for image in images:archive.write(image,arcname=f'field-{field}/'+image.name)
    print(f'Packaged {label}: actual PNGs and integration times.',flush=True)
