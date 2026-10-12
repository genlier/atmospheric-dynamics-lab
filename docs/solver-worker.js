'use strict';
importScripts('solver.js');
self.onmessage=function(event){
  try{
    const {initial,n,Re,beta,T,maxDt,count}=event.data;
    if(!Number.isFinite(T)||T<=0||T>30||!Number.isFinite(maxDt)||maxDt<0.00001||maxDt>0.1||!Number.isInteger(count)||count<21||count>241)throw Error('Invalid duration, step or frame count');
    const model=new self.SpectralBVE.BVE(initial,n,Re,beta);
    self.postMessage({type:'started',mean:model.removedMean});
    let progressAt=performance.now();
    for(let index=0;index<count;index++){
      const target=T*index/(count-1);
      while(target-model.time>1e-12){
        if(model.steps>=250000)throw Error('Step limit reached: use a shorter duration or a larger maximum time step');
        model.step(maxDt,target-model.time);
        if(performance.now()-progressAt>150){self.postMessage({type:'progress',time:model.time,steps:model.steps});progressAt=performance.now();}
      }
      const frame=model.snapshot();self.postMessage({type:'frame',index,frame,steps:model.steps},[frame.values.buffer]);
    }
    self.postMessage({type:'done',steps:model.steps});
  }catch(error){self.postMessage({type:'error',message:error.message});}
};
