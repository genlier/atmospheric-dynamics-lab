'use strict';
(function(root){
  const table=Uint32Array.from({length:256},(_,i)=>{let c=i;for(let j=0;j<8;j++)c=c&1?0xedb88320^(c>>>1):c>>>1;return c>>>0;});
  function crc32(bytes){let crc=0xffffffff;for(const byte of bytes)crc=table[(crc^byte)&255]^(crc>>>8);return (crc^0xffffffff)>>>0;}
  function header(size){const a=new Uint8Array(size);return [a,new DataView(a.buffer)];}
  function zipStore(files){const parts=[],central=[];let offset=0,total=0;const encode=new TextEncoder();
    for(const file of files){const name=encode.encode(file.name),bytes=file.bytes,crc=crc32(bytes);const [a,v]=header(30);v.setUint32(0,0x04034b50,true);v.setUint16(4,20,true);v.setUint16(6,0x800,true);v.setUint16(12,0x21,true);v.setUint32(14,crc,true);v.setUint32(18,bytes.length,true);v.setUint32(22,bytes.length,true);v.setUint16(26,name.length,true);parts.push(a,name,bytes);
      const [b,w]=header(46);w.setUint32(0,0x02014b50,true);w.setUint16(4,20,true);w.setUint16(6,20,true);w.setUint16(8,0x800,true);w.setUint16(14,0x21,true);w.setUint32(16,crc,true);w.setUint32(20,bytes.length,true);w.setUint32(24,bytes.length,true);w.setUint16(28,name.length,true);w.setUint32(42,offset,true);central.push(b,name);total+=46+name.length;offset+=30+name.length+bytes.length;
    }
    const [end,e]=header(22);e.setUint32(0,0x06054b50,true);e.setUint16(8,files.length,true);e.setUint16(10,files.length,true);e.setUint32(12,total,true);e.setUint32(16,offset,true);return new Blob([...parts,...central,end],{type:'application/zip'});
  }
  root.FrameArchive={zipStore,crc32};
  if(typeof module!=='undefined')module.exports=root.FrameArchive;
})(typeof self!=='undefined'?self:globalThis);
