'use strict';
(function(root){
  class FFT {
    constructor(n){
      this.n=n;this.reverse=new Uint16Array(n);this.cos=new Float64Array(n/2);this.sin=new Float64Array(n/2);
      const bits=Math.log2(n);
      for(let i=0;i<n;i++){let v=i,r=0;for(let b=0;b<bits;b++){r=(r<<1)|(v&1);v>>=1;}this.reverse[i]=r;}
      for(let i=0;i<n/2;i++){this.cos[i]=Math.cos(2*Math.PI*i/n);this.sin[i]=Math.sin(2*Math.PI*i/n);}
      this.rowR=new Float64Array(n);this.rowI=new Float64Array(n);
    }
    line(re,im,inverse){
      const n=this.n;
      for(let i=0;i<n;i++){const j=this.reverse[i];if(j>i){[re[i],re[j]]=[re[j],re[i]];[im[i],im[j]]=[im[j],im[i]];}}
      for(let size=2;size<=n;size*=2){const half=size/2,stride=n/size;
        for(let base=0;base<n;base+=size)for(let j=0;j<half;j++){
          const c=this.cos[j*stride],s=this.sin[j*stride]*(inverse?1:-1),a=base+j,b=a+half;
          const r=c*re[b]-s*im[b],v=s*re[b]+c*im[b],ar=re[a],ai=im[a];
          re[a]=ar+r;im[a]=ai+v;re[b]=ar-r;im[b]=ai-v;
        }
      }
      if(inverse)for(let i=0;i<n;i++){re[i]/=n;im[i]/=n;}
    }
    transform(re,im,inverse=false){
      const n=this.n,r=this.rowR,v=this.rowI;
      for(let y=0;y<n;y++){for(let x=0;x<n;x++){r[x]=re[x+y*n];v[x]=im[x+y*n];}this.line(r,v,inverse);for(let x=0;x<n;x++){re[x+y*n]=r[x];im[x+y*n]=v[x];}}
      for(let x=0;x<n;x++){for(let y=0;y<n;y++){r[y]=re[x+y*n];v[y]=im[x+y*n];}this.line(r,v,inverse);for(let y=0;y<n;y++){re[x+y*n]=r[y];im[x+y*n]=v[y];}}
    }
  }
  class BVE {
    constructor(initial,n,Re,beta){
      if(![32,64].includes(n)||initial.length!==n*n||!Number.isFinite(Re)||Re<=0||!Number.isFinite(beta))throw Error('Invalid grid or parameters');
      this.n=n;this.size=n*n;this.Re=Re;this.beta=beta;this.fft=new FFT(n);this.re=Float64Array.from(initial);this.im=new Float64Array(this.size);this.time=0;this.steps=0;
      if(!this.re.every(Number.isFinite))throw Error('Initial field contains non-finite values');
      const mean=this.re.reduce((a,b)=>a+b,0)/this.size;this.removedMean=mean;for(let i=0;i<this.size;i++)this.re[i]-=mean;
      this.fft.transform(this.re,this.im);this.k2=new Float64Array(this.size);this.kx=new Float64Array(this.size);this.ky=new Float64Array(this.size);this.mask=new Uint8Array(this.size);
      for(let y=0;y<n;y++)for(let x=0;x<n;x++){const i=x+y*n,kx=x<n/2?x:x-n,ky=y<n/2?y:y-n;this.kx[i]=kx;this.ky[i]=ky;this.k2[i]=kx*kx+ky*ky;this.mask[i]=Math.abs(kx)<=n/3&&Math.abs(ky)<=n/3?1:0;}
      this.work=Array.from({length:4},()=>[new Float64Array(this.size),new Float64Array(this.size)]);
      this.stageR=new Float64Array(this.size);this.stageI=new Float64Array(this.size);
      this.ks=Array.from({length:4},()=>[new Float64Array(this.size),new Float64Array(this.size)]);
      this.lastSpeed=0;
    }
    rhs(re,im,outR,outI){
      const [u,v,zx,zy]=this.work;
      for(let i=0;i<this.size;i++){
        const inv=this.k2[i]? -1/this.k2[i]:0,m=this.mask[i],pr=re[i]*inv,pi=im[i]*inv;
        u[0][i]=this.ky[i]*pi*m;u[1][i]=-this.ky[i]*pr*m;
        v[0][i]=-this.kx[i]*pi*m;v[1][i]=this.kx[i]*pr*m;
        zx[0][i]=-this.kx[i]*im[i]*m;zx[1][i]=this.kx[i]*re[i]*m;
        zy[0][i]=-this.ky[i]*im[i]*m;zy[1][i]=this.ky[i]*re[i]*m;
      }
      for(const w of this.work)this.fft.transform(w[0],w[1],true);
      let speed=0;
      for(let i=0;i<this.size;i++){outR[i]=-(u[0][i]*zx[0][i]+v[0][i]*zy[0][i]);outI[i]=0;speed=Math.max(speed,Math.abs(u[0][i])+Math.abs(v[0][i]));}
      this.lastSpeed=speed;this.fft.transform(outR,outI);
      for(let i=0;i<this.size;i++){
        const inv=this.k2[i]?-1/this.k2[i]:0;
        outR[i]=outR[i]*this.mask[i]+this.beta*this.kx[i]*im[i]*inv;
        outI[i]=outI[i]*this.mask[i]-this.beta*this.kx[i]*re[i]*inv;
      }
    }
    step(maxDt,remaining){
      const [k1,k2,k3,k4]=this.ks;this.rhs(this.re,this.im,...k1);
      const dt=Math.min(maxDt,remaining,0.25*(2*Math.PI/this.n)/Math.max(this.lastSpeed,1e-12),0.25/Math.max(Math.abs(this.beta),1e-12));
      for(let stage=1;stage<4;stage++){
        const previous=this.ks[stage-1],h=stage===3?dt:dt/2;
        for(let i=0;i<this.size;i++){this.stageR[i]=this.re[i]+h*previous[0][i];this.stageI[i]=this.im[i]+h*previous[1][i];}
        this.rhs(this.stageR,this.stageI,...this.ks[stage]);
      }
      for(let i=0;i<this.size;i++){
        const decay=Math.exp(-this.k2[i]*dt/this.Re);
        this.re[i]=(this.re[i]+dt/6*(k1[0][i]+2*k2[0][i]+2*k3[0][i]+k4[0][i]))*decay;
        this.im[i]=(this.im[i]+dt/6*(k1[1][i]+2*k2[1][i]+2*k3[1][i]+k4[1][i]))*decay;
        if(!Number.isFinite(this.re[i])||!Number.isFinite(this.im[i]))throw Error('Numerical overflow: reduce the amplitude or maximum time step');
      }
      this.time+=dt;this.steps++;
    }
    snapshot(){
      const r=this.re.slice(),v=this.im.slice();this.fft.transform(r,v,true);
      let energy=0,enstrophy=0;
      for(let i=0;i<this.size;i++){const power=this.re[i]**2+this.im[i]**2;energy+=this.k2[i]?power/this.k2[i]:0;enstrophy+=r[i]**2;}
      return {time:this.time,values:r,energy:energy/(2*this.size**2),enstrophy:enstrophy/(2*this.size)};
    }
  }
  root.SpectralBVE={FFT,BVE};
  if(typeof module!=='undefined')module.exports=root.SpectralBVE;
})(typeof self!=='undefined'?self:globalThis);
