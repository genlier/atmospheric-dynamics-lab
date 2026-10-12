'use strict';
(function(){
  let n=64,initial=new Float64Array(n*n),brush='positive',drawing=false,lastPoint=null;
  let worker=null,frames=[],result=null,position=0,animation=null,range=[Infinity,-Infinity],currentValues=null;
  const editor=$('initial-field'),canvas=$('live-flow');
  let exporting=false;
  const presets={vortex:'40*exp(-((x+1)^2+y^2)/0.16)+40*exp(-((x-1)^2+y^2)/0.16)',grid:'-2*cos(x)*cos(y)',dipole:'20*exp(-((x+0.65)^2+y^2)/0.25)-20*exp(-((x-0.65)^2+y^2)/0.25)',blank:'0'};
  function status(text,error=false){$('simulation-status').textContent=text;$('simulation-status').classList.toggle('error',error);}
  function dirty(message){stopLive();status(frames.length?t('stale'):message);}
  function validateValues(values){if(!values.every(v=>typeof v==='number'&&Number.isFinite(v)&&Math.abs(v)<=1000))throw Error(t('invalidValues'));}
  function preview(){const [lo,hi]=paintField(editor,initial,n);$('initial-range').textContent=`${plain(Math.min(...initial))} … ${plain(Math.max(...initial))}`;}
  function applyExpression(){const f=FieldExpression.compile($('initial-expression').value),a=new Float64Array(n*n);for(let y=0;y<n;y++)for(let x=0;x<n;x++)a[x+y*n]=f(-Math.PI+2*Math.PI*x/n,-Math.PI+2*Math.PI*y/n);validateValues(a);initial=a;preview();dirty(t('formulaLoaded'));}
  function loadPreset(){try{$('initial-expression').value=presets[$('preset').value];applyExpression();}catch(e){status(t('error')+': '+translatedError(e.message),true);}}
  function point(event){const r=editor.getBoundingClientRect();return {x:-Math.PI+2*Math.PI*Math.max(0,Math.min(1,(event.clientX-r.left)/r.width)),y:Math.PI-2*Math.PI*Math.max(0,Math.min(1,(event.clientY-r.top)/r.height))};}
  function stamp(p,weight=1){const width=Number($('brush-width').value),strength=Number($('brush-strength').value);if(!Number.isFinite(strength)||strength<0.1||strength>100)return;for(let y=0;y<n;y++)for(let x=0;x<n;x++){let dx=Math.abs(-Math.PI+2*Math.PI*x/n-p.x),dy=Math.abs(-Math.PI+2*Math.PI*y/n-p.y);dx=Math.min(dx,2*Math.PI-dx);dy=Math.min(dy,2*Math.PI-dy);const g=Math.exp(-(dx*dx+dy*dy)/(width*width)),i=x+y*n;if(brush==='erase')initial[i]*=1-Math.min(1,g*weight);else initial[i]=Math.max(-1000,Math.min(1000,initial[i]+(brush==='positive'?1:-1)*strength*g*weight));}}
  editor.addEventListener('pointerdown',event=>{if(worker||event.button>0)return;event.preventDefault();drawing=true;lastPoint=point(event);editor.setPointerCapture(event.pointerId);stamp(lastPoint);preview();dirty(t('drawing'));});
  editor.addEventListener('pointermove',event=>{const p=point(event);$('editor-probe').textContent=`x ${p.x.toFixed(2)} · y ${p.y.toFixed(2)}`;if(!drawing)return;const dx=p.x-lastPoint.x,dy=p.y-lastPoint.y,count=Math.max(1,Math.ceil(Math.hypot(dx,dy)/(2*Math.PI/n)));for(let i=1;i<=count;i++)stamp({x:lastPoint.x+dx*i/count,y:lastPoint.y+dy*i/count},0.2);lastPoint=p;preview();});
  const finishDrawing=()=>{drawing=false;lastPoint=null;};editor.addEventListener('pointerup',finishDrawing);editor.addEventListener('pointercancel',finishDrawing);editor.addEventListener('lostpointercapture',finishDrawing);
  document.querySelectorAll('[data-brush]').forEach(button=>button.onclick=()=>{brush=button.dataset.brush;document.querySelectorAll('[data-brush]').forEach(b=>{b.classList.toggle('active',b===button);b.setAttribute('aria-pressed',String(b===button));});});
  $('clear-initial').onclick=()=>{initial.fill(0);preview();dirty(t('drawing'));};$('load-preset').onclick=loadPreset;$('apply-expression').onclick=()=>{try{applyExpression();}catch(e){status(t('error')+': '+translatedError(e.message),true);}};
  $('sim-grid').onchange=()=>{n=Number($('sim-grid').value);initial=new Float64Array(n*n);loadPreset();dirty(t('gridChanged'));};
  document.querySelectorAll('.parameter-grid input,.parameter-grid select').forEach(input=>input.addEventListener('change',()=>dirty(t('ready'))));
  function saveJSON(value,name){const blob=new Blob([JSON.stringify(value)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  $('export-initial').onclick=()=>saveJSON({schema:'bve-initial-v1',n,domain:[-Math.PI,Math.PI],layout:'x-fastest, y-increasing',values:Array.from(initial)},'initial-vorticity.json');
  $('import-field').onchange=async event=>{try{const file=event.target.files[0];if(!file)return;if(file.size>5*1024*1024)throw Error(t('invalidField'));const obj=JSON.parse(await file.text()),size=obj.n,values=obj.values||obj.initial;if(![32,64].includes(size)||!Array.isArray(values)||values.length!==size*size)throw Error(t('invalidField'));validateValues(values);n=size;initial=Float64Array.from(values);$('sim-grid').value=n;preview();dirty(t('imported'));}catch(e){status(t('error')+': '+translatedError(e.message),true);}finally{event.target.value='';}};
  function stopLive(){if(animation)cancelAnimationFrame(animation);animation=null;$('live-play').textContent=t('play');}
  function drawLive(){
    if(!frames.length)return;position=Math.max(0,Math.min(frames.length-1,position));const i=Math.round(position),a=frames[i];
    currentValues=a.values;
    const bounds=paintField(canvas,currentValues,result.n,$('live-fixed').checked?range:null);
    $('result-empty').hidden=true;$('live-scale-min').innerHTML=number(bounds[0]);$('live-scale-max').innerHTML=number(bounds[1]);$('live-time').textContent=`t = ${plain(a.time)}`;
    $('live-energy').innerHTML=number(a.energy,5);$('live-enstrophy').innerHTML=number(a.enstrophy,5);$('live-frame').textContent=`${i+1} / ${frames.length}`;$('live-timeline').value=position;drawEnergy();
  }
  function drawEnergy(){const el=$('live-chart'),c=el.getContext('2d'),W=el.width,H=el.height,L=72,R=20,T=20,B=30,ys=frames.map(f=>f.energy);let lo=0,hi=Math.max(...ys)*1.08||1;const end=result.parameters.T;
    c.clearRect(0,0,W,H);c.font='13px Consolas, monospace';c.strokeStyle='#294045';c.fillStyle='#93acae';
    for(let i=0;i<3;i++){const y=T+(H-T-B)*i/2;c.beginPath();c.moveTo(L,y);c.lineTo(W-R,y);c.stroke();c.textAlign='right';c.fillText(plain(hi*(1-i/2)),L-8,y+4);}
    c.textAlign='left';c.fillText('0',L,H-8);c.textAlign='right';c.fillText(`t = ${plain(end)}`,W-R,H-8);c.strokeStyle='#75d5da';c.lineWidth=2;c.beginPath();frames.forEach((f,i)=>{const x=L+f.time/end*(W-L-R),y=H-B-f.energy/hi*(H-T-B);if(i)c.lineTo(x,y);else c.moveTo(x,y);});c.stroke();
    const f=frames[Math.round(position)];c.fillStyle='#c9ec83';c.beginPath();c.arc(L+f.time/end*(W-L-R),H-B-f.energy/hi*(H-T-B),3,0,2*Math.PI);c.fill();
  }
  function playLive(){if(animation){stopLive();return;}if(frames.length<2)return;if(position>=frames.length-1)position=0;$('live-play').textContent=t('pause');let last=performance.now();const tick=now=>{const elapsed=Math.max(0,Math.min(now-last,100));last=now;position=Math.min(frames.length-1,position+elapsed/10000*(frames.length-1)*Number($('playback-speed').value));drawLive();if(position>=frames.length-1){stopLive();return;}animation=requestAnimationFrame(tick);};animation=requestAnimationFrame(tick);}
  $('live-play').onclick=playLive;$('live-timeline').oninput=event=>{stopLive();position=Math.round(Number(event.target.value));drawLive();};$('live-fixed').onchange=drawLive;
  canvas.addEventListener('pointermove',event=>{if(!currentValues)return;const r=canvas.getBoundingClientRect(),N=result.n,x=Math.min(N-1,Math.max(0,Math.floor((event.clientX-r.left)/r.width*N))),y=N-1-Math.min(N-1,Math.max(0,Math.floor((event.clientY-r.top)/r.height*N)));$('live-probe').textContent=`x ${(-Math.PI+2*Math.PI*x/N).toFixed(2)} · y ${(-Math.PI+2*Math.PI*y/N).toFixed(2)} · ${plain(currentValues[x+y*N])}`;});
  function lock(running){$('run-simulation').disabled=running;$('cancel-simulation').disabled=!running;document.querySelector('.initial-panel').inert=running;document.querySelectorAll('.parameter-grid input,.parameter-grid select').forEach(el=>el.disabled=running);$('use-reference').disabled=running;}
  function finish(message,error=false){if(worker)worker.terminate();worker=null;lock(false);$('live-play').disabled=frames.length<2;$('live-timeline').disabled=frames.length<2;$('download-result').disabled=!frames.length;$('export-png').disabled=exporting||!frames.length;status(message,error);}
  $('simulation-form').onsubmit=event=>{
    event.preventDefault();if(worker)return;try{
      if(!$('simulation-form').reportValidity())return;validateValues(initial);
      const parameters={Re:Number($('sim-re').value),beta:Number($('sim-beta').value),T:Number($('sim-duration').value),maxDt:Number($('sim-dt').value),count:Number($('sim-frames').value)};
      worker=new Worker('solver-worker.js?v=3');stopLive();frames=[];position=0;range=[Infinity,-Infinity];result={schema:'bve-result-v1',n,domain:[-Math.PI,Math.PI],layout:'x-fastest, y-increasing',initial:Array.from(initial),parameters,removedMean:0,steps:0,complete:false};lock(true);$('simulation-progress').value=0;$('live-play').disabled=true;$('live-timeline').disabled=true;$('download-result').disabled=true;status(t('progress'));
      worker.onmessage=event=>{const data=event.data;
        if(data.type==='started'){result.removedMean=data.mean;}
        if(data.type==='progress'){result.steps=data.steps;$('simulation-progress').value=data.time/parameters.T;status(`${t('progress')} · t ${plain(data.time)} · ${data.steps} ${t('steps')}`);}
        if(data.type==='frame'){frames.push(data.frame);result.steps=data.steps;for(const v of data.frame.values){range[0]=Math.min(range[0],v);range[1]=Math.max(range[1],v);}$('live-timeline').max=frames.length-1;position=frames.length-1;$('simulation-progress').value=data.frame.time/parameters.T;drawLive();}
        if(data.type==='error')finish(t('runFailed')+': '+translatedError(data.message),true);
        if(data.type==='done'){result.complete=true;result.steps=data.steps;finish(`${t('done')} · ${frames.length} ${t('frames')} · ${data.steps} ${t('steps')} · ${t('initialMean')}: ${plain(result.removedMean)}`);position=0;drawLive();playLive();}
      };
      worker.onerror=()=>finish(t('workerError'),true);worker.postMessage({initial:Array.from(initial),n,...parameters});
    }catch(e){finish(t('error')+': '+translatedError(e.message),true);}
  };
  $('cancel-simulation').onclick=()=>finish(t('stopped'));
  $('download-result').onclick=()=>{if(!result||!frames.length)return;saveJSON({...result,frames:frames.map(f=>({...f,values:Array.from(f.values)}))},'vorticity-evolution.json');};
  $('export-png').onclick=async()=>{
    if(exporting||!frames.length)return;exporting=true;$('export-png').disabled=true;stopLive();
    const savedFrames=frames.slice(),savedResult=result,bounds=range.slice(),fixed=$('live-fixed').checked,output=document.createElement('canvas');output.width=640;output.height=640;const files=[];
    try{for(let i=0;i<savedFrames.length;i++){
      const f=savedFrames[i];paintField(output,f.values,savedResult.n,fixed?bounds:null);const c=output.getContext('2d');c.fillStyle='#091416db';c.fillRect(10,10,300,30);c.fillStyle='#ffffff';c.font='16px monospace';c.fillText(`t = ${plain(f.time)} · ${i+1}/${savedFrames.length}`,18,31);
      const blob=await new Promise(resolve=>output.toBlob(resolve,'image/png'));if(!blob)throw Error('PNG export failed');files.push({name:`frame_${String(i).padStart(4,'0')}.png`,bytes:new Uint8Array(await blob.arrayBuffer())});$('export-png').textContent=`PNG ${i+1}/${savedFrames.length}`;
    }
    files.push({name:'times.csv',bytes:new TextEncoder().encode('frame,time\n'+savedFrames.map((f,i)=>`${i},${f.time}`).join('\n'))});files.push({name:'simulation.json',bytes:new TextEncoder().encode(JSON.stringify(savedResult))});
    const blob=FrameArchive.zipStore(files),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='computed-png-frames.zip';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
    }catch(error){status(t('error')+': '+translatedError(error.message),true);}finally{exporting=false;$('export-png').disabled=!!worker||!frames.length;$('export-png').textContent='PNG ZIP ↓';}
  };
  $('use-reference').onclick=()=>{const c=window.LAB_DATA.cases[selected];n=c.displayGrid;initial=Float64Array.from(c.frames[0].fields[0]);$('sim-grid').value=n;$('sim-re').value=c.Re;$('sim-beta').value=c.beta;$('sim-duration').value=Math.min(30,c.frames[c.frames.length-1].time);preview();dirty(t('fieldLoaded')+' · '+t('useReferenceNote'));document.getElementById('live-lab').scrollIntoView({behavior:'smooth'});};
  document.addEventListener('visibilitychange',()=>{if(document.hidden)stopLive();});
  loadPreset();
})();





