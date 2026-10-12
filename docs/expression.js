'use strict';
(function(root){
  const functions={sin:Math.sin,cos:Math.cos,tan:Math.tan,tanh:Math.tanh,exp:Math.exp,log:Math.log,sqrt:Math.sqrt,abs:Math.abs};
  function compile(source){
    source=source.replace(/−/g,'-');
    if(source.length>800)throw Error('Expression too long');
    const tokens=[];let offset=0;
    while(offset<source.length){
      if(/\s/.test(source[offset])){offset++;continue;}
      const m=/^(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[a-zA-Z]+|[+\-*/^()])/.exec(source.slice(offset));
      if(!m)throw Error('Unexpected character at '+(offset+1));tokens.push(m[0]);offset+=m[0].length;
    }
    let index=0;
    function primary(){const token=tokens[index++];
      if(token==='('){const value=sum();if(tokens[index++]!==')')throw Error('Missing closing parenthesis');return value;}
      if(token==='x')return (x,y)=>x;if(token==='y')return (x,y)=>y;if(token==='pi')return ()=>Math.PI;
      if(Object.hasOwn(functions,token)){if(tokens[index++]!=='(')throw Error('Function requires parentheses');const value=sum();if(tokens[index++]!==')')throw Error('Missing closing parenthesis');return (x,y)=>functions[token](value(x,y));}
      if(token!==undefined&&/^(?:\d|\.)/.test(token)){const value=Number(token);if(!Number.isFinite(value))throw Error('Invalid number');return ()=>value;}
      throw Error('Expected a number, x, y or function');
    }
    function power(){const left=primary();if(tokens[index]==='^'){index++;const right=unary();return (x,y)=>left(x,y)**right(x,y);}return left;}
    function unary(){if(tokens[index]==='-'){index++;const v=unary();return (x,y)=>-v(x,y);}if(tokens[index]==='+'){index++;return unary();}return power();}
    function product(){let left=unary();while(['*','/'].includes(tokens[index])){const op=tokens[index++],a=left,b=unary();left=op==='*'?(x,y)=>a(x,y)*b(x,y):(x,y)=>a(x,y)/b(x,y);}return left;}
    function sum(){let left=product();while(['+','-'].includes(tokens[index])){const op=tokens[index++],a=left,b=product();left=op==='+'?(x,y)=>a(x,y)+b(x,y):(x,y)=>a(x,y)-b(x,y);}return left;}
    const evaluate=sum();if(index!==tokens.length)throw Error('Unexpected token');return evaluate;
  }
  root.FieldExpression={compile};
  if(typeof module!=='undefined')module.exports=root.FieldExpression;
})(typeof self!=='undefined'?self:globalThis);
