'use strict';
const EPS=0n;
const ELEMENTS=new Set('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split(' '));
function gcd(a,b){a=BigInt(a);b=BigInt(b);while(b){let t=a%b;a=b;b=t}return a<0n?-a:a}
function lcm(a,b){a=BigInt(a);b=BigInt(b);return a===0n||b===0n?0n:(a/gcd(a,b))*b}
class F{constructor(n,d=1n){n=BigInt(n);d=BigInt(d);if(d===0n)throw Error('Division by zero');if(d<0n){n=-n;d=-d}const g=gcd(n,d);this.n=n/g;this.d=d/g}add(o){o=toF(o);return new F(this.n*o.d+o.n*this.d,this.d*o.d)}sub(o){o=toF(o);return new F(this.n*o.d-o.n*this.d,this.d*o.d)}mul(o){o=toF(o);return new F(this.n*o.n,this.d*o.d)}div(o){o=toF(o);return new F(this.n*o.d,this.d*o.n)}neg(){return new F(-this.n,this.d)}isZero(){return this.n===0n}toString(){return this.d===1n?this.n.toString():this.n+'/'+this.d}}
function toF(x){return x instanceof F?x:new F(x)}
function normalize(s){return s.replace(/→|⟶|⟹/g,'->').replace(/⇌|⇄/g,'->').replace(/=/g,'->').replace(/•/g,'·').trim()}
function readNumber(s,i){let j=i;while(j<s.length&&/[0-9]/.test(s[j]))j++;return {value:j===i?1:Number(s.slice(i,j)),i:j}}
function parseFormula(formula){
 if(formula.replace(/\s+/g,'')==='e-'||formula.replace(/\s+/g,'')==='e')return {formula,atoms:{},charge:-1};
 let s=formula.replace(/\s+/g,''),i=0; if(!s)throw Error('Empty chemical formula.');
 let charge=0;
 const chargeMatch=s.match(/(?:\^([0-9]*)([+-])|^([A-Z][a-z]?)([0-9]+)([+-])$|([+-]))$/); if(chargeMatch){const sign=chargeMatch[2]||chargeMatch[5]||chargeMatch[6];const mag=chargeMatch[2]?Number(chargeMatch[1]||1):(chargeMatch[5]?Number(chargeMatch[4]):1);charge=(sign==='+'?1:-1)*mag;s=chargeMatch[5]?chargeMatch[3]:s.slice(0,chargeMatch.index)}
 if(/\^/.test(s))throw Error('Charge notation must look like SO4^2- or NH4+.');
 function group(stop){
  const out={};
  while(i<s.length){
   if(stop&&s[i]===stop){i++;return out}
   if(s[i]==='('||s[i]==='['||s[i]==='{'){
    const open=s[i++],close=open==='('?')':open==='['?']':'}';
    const inner=group(close);const num=readNumber(s,i);i=num.i;
    for(const e in inner)out[e]=(out[e]||0)+inner[e]*num.value;
    continue;
   }
   if(s[i]===')'||s[i]===']'||s[i]==='}')throw Error('Unmatched closing bracket in '+formula);
   if(!/[A-Z]/.test(s[i]))throw Error('Expected an element symbol near "'+s.slice(i)+'" in '+formula);
   let e=s[i++];while(i<s.length&&/[a-z]/.test(s[i]))e+=s[i++];if(!ELEMENTS.has(e))throw Error('Unknown element symbol: '+e);
   const num=readNumber(s,i);i=num.i;out[e]=(out[e]||0)+num.value;
  }
  if(stop)throw Error('Missing closing bracket in '+formula);return out;
 }
 const atoms=group(null);if(!Object.keys(atoms).length)throw Error('No elements found in '+formula);return {formula,atoms,charge};
}
function splitHydrate(raw){return raw.split(/[·.]/).filter(Boolean)}
function parseCompound(token){
 let raw=token.trim(); if(!raw)throw Error('Empty compound.');
 let m=raw.match(/^(\d+)\s*/),given=m?Number(m[1]):null; if(m)raw=raw.slice(m[0].length); if(!raw)throw Error('A coefficient cannot be used without a formula.');
 const total={},parts=splitHydrate(raw);let charge=0;for(const p of parts){const q=p.match(/^(\d+)(?=[A-Z\(\[])/);let mult=q?Number(q[1]):1;const body=q? p.slice(q[0].length):p;const parsed=parseFormula(body);for(const e in parsed.atoms)total[e]=(total[e]||0)+parsed.atoms[e]*mult;charge+=parsed.charge*mult}
 return {raw,given,atoms:total,charge}
}
function parseSide(side){const parts=side.split(/\s+\+\s+|\+\+(?=\s*[A-Z(\[]|\s*e-)|\+(?=\s*(?:[A-Z(\[]|e-))/);if(parts.some(x=>!x.trim()))throw Error('Empty compound near +.');return parts.map(parseCompound)}
function matrixFor(L,R){const elements=[...new Set([...L,...R].flatMap(c=>Object.keys(c.atoms)))];const hasCharge=L.some(c=>c.charge)||R.some(c=>c.charge);if(hasCharge)elements.push('Charge');const cols=[...L,...R];const M=elements.map(e=>cols.map((c,j)=>{const v=e==='Charge'?c.charge:(c.atoms[e]||0);return j<L.length?v:-v}));return {M,elements,cols}}
function rref(A){A=A.map(row=>row.map(toF));let r=0,piv=[];for(let c=0;c<A[0].length&&r<A.length;c++){let p=r;while(p<A.length&&A[p][c].isZero())p++;if(p===A.length)continue;[A[p],A[r]]=[A[r],A[p]];const inv=A[r][c];for(let j=c;j<A[0].length;j++)A[r][j]=A[r][j].div(inv);for(let k=0;k<A.length;k++)if(k!==r&&!A[k][c].isZero()){const f=A[k][c];for(let j=c;j<A[0].length;j++)A[k][j]=A[k][j].sub(f.mul(A[r][j]))}piv.push(c);r++}return {A,piv}}
function nullBasis(M){const rr=rref(M),n=M[0].length,free=[];for(let c=0;c<n;c++)if(!rr.piv.includes(c))free.push(c);if(!free.length)return null;return free.map(f=>{const v=Array.from({length:n},()=>new F(0n));v[f]=new F(1n);for(let r=0;r<rr.piv.length;r++)v[rr.piv[r]]=rr.A[r][f].neg();return v})}
function integerVector(v){let L=1n;for(const x of v)L=lcm(L,x.d);let a=v.map(x=>x.n*(L/x.d));let g=0n;for(const x of a)g=gcd(g,x);a=a.map(x=>x/g);if(a.every(x=>x<0n))a=a.map(x=>-x);return a}
function findPositiveVector(basis){
 if(basis.length===1){const v=integerVector(basis[0]);return v.every(x=>x>0n)?v:null}
 // Search small combinations for underdetermined systems; deterministic and bounded.
 const k=basis.length, range=8;function rec(pos,coef){if(pos===k){const v=basis[0].map(()=>new F(0n));for(let q=0;q<k;q++)for(let i=0;i<v.length;i++)v[i]=v[i].add(basis[q][i].mul(coef[q]));if(v.every(x=>x.n>0n))return integerVector(v);return null}for(let x=-range;x<=range;x++){if(x===0)continue;coef[pos]=x;const got=rec(pos+1,coef);if(got)return got}return null}return rec(0,[])
}
function balanceEquation(input){
 const s=normalize(input);const parts=s.split('->');if(parts.length!==2)throw Error('Use exactly one reaction arrow (for example ->).');const L=parseSide(parts[0]),R=parseSide(parts[1]);if(!L.length||!R.length)throw Error('Both sides need at least one compound.');const {M,elements}=matrixFor(L,R);const basis=nullBasis(M);if(!basis)throw Error('This equation has no non-zero conservation solution.');const co=findPositiveVector(basis);if(!co)throw Error('The equation is underdetermined or cannot produce one positive integer coefficient set with the current constraints.');const original=[...L,...R].map(c=>c.given);return {L,R,co,elements,original,input:s,multi:basis.length>1};
}
function verify(b){return b.elements.map(e=>{const l=b.L.reduce((s,c,i)=>s+b.co[i]*BigInt(e==='Charge'?c.charge:(c.atoms[e]||0)),0n);const r=b.R.reduce((s,c,i)=>s+b.co[b.L.length+i]*BigInt(e==='Charge'?c.charge:(c.atoms[e]||0)),0n);return {e,l,r,ok:l===r}})}
function formulaHTML(raw){return esc(raw).replace(/\^(\d*[+-])/g,'<sup>$1</sup>').replace(/(?<=[A-Za-z0-9)\]}])([+-])$/,'<sup>$1</sup>').replace(/([A-Za-z)\]}])(\d+)/g,'$1<sub>$2</sub>')}
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const MASS={};'H1.008 He4.0026 Li6.94 Be9.0122 B10.81 C12.011 N14.007 O15.999 F18.998 Ne20.18 Na22.99 Mg24.305 Al26.982 Si28.085 P30.974 S32.06 Cl35.45 Ar39.948 K39.098 Ca40.078 Sc44.956 Ti47.867 V50.942 Cr51.996 Mn54.938 Fe55.845 Co58.933 Ni58.693 Cu63.546 Zn65.38 Ga69.723 Ge72.63 As74.922 Se78.971 Br79.904 Kr83.798 Rb85.468 Sr87.62 Y88.906 Zr91.224 Nb92.906 Mo95.95 Ru101.07 Rh102.91 Pd106.42 Ag107.87 Cd112.41 In114.82 Sn118.71 Sb121.76 Te127.6 I126.9 Xe131.29 Cs132.91 Ba137.33 La138.91 Ce140.12 W183.84 Pt195.08 Au196.97 Hg200.59 Tl204.38 Pb207.2 Bi208.98 Th232.04 U238.03'.split(' ').forEach(t=>{const m=t.match(/^([A-Z][a-z]?)([\d.]+)$/);MASS[m[1]]=parseFloat(m[2])});
function molar(c){let t=0;for(const e in c.atoms){if(!MASS[e])return null;t+=MASS[e]*c.atoms[e]}return t>0?t:null}
function checkMsg(b){if(!b.original.every(x=>x!==null))return'';const ok=b.original.every((x,i)=>BigInt(x)===b.co[i]);return `<p><span class="status ${ok?'ok':'bad'}">${ok?'YOUR COEFFICIENTS ARE CORRECT':'YOUR COEFFICIENTS DIFFER FROM THE BALANCED SET'}</span></p>`}
function renderStoich(b){const box=$('#stoich'),sp=[...b.L,...b.R],n=b.L.length,mm=sp.map(molar);if(!mm.some(Boolean)){box.classList.add('hidden');return}box.classList.remove('hidden');
 const rows=sp.map((c,i)=>`<tr><td>${formulaHTML(c.raw)}</td><td>${b.co[i]}</td><td>${mm[i]?mm[i].toFixed(3):'—'}</td><td>${i<n?'Reactant':'Product'}</td><td>${mm[i]?`<input class="input" type="number" min="0" step="any" data-i="${i}" placeholder="grams" aria-label="Grams of ${esc(c.raw)}">`:'—'}</td><td data-out="${i}" class="muted">—</td></tr>`).join('');
 box.innerHTML=`<h2>Stoichiometry &amp; limiting reagent</h2><p class="muted">Enter grams for one or more species. The limiting reagent, leftover reactants and theoretical masses update live.</p><div class="scroll"><table><thead><tr><th>Species</th><th>Coef.</th><th>Molar mass (g/mol)</th><th>Role</th><th>Given (g)</th><th>Result</th></tr></thead><tbody>${rows}</tbody></table></div><p id="stoichNote" class="hint"></p>`;
 const calc=()=>{const g=sp.map((c,i)=>{const el=box.querySelector(`input[data-i="${i}"]`),v=el?parseFloat(el.value):NaN;return v>0?v:0}),cs=b.co.map(Number);let ext=null,lim=-1;
  for(let i=0;i<n;i++)if(g[i]){const x=g[i]/mm[i]/cs[i];if(ext===null||x<ext){ext=x;lim=i}}
  if(ext===null)for(let i=n;i<sp.length;i++)if(g[i]){ext=g[i]/mm[i]/cs[i];lim=i;break}
  sp.forEach((c,i)=>{const o=box.querySelector(`[data-out="${i}"]`);if(ext===null||!mm[i]){o.textContent='—';return}const m=ext*cs[i]*mm[i];o.textContent=i<n?(g[i]?`uses ${m.toFixed(3)} g · left ${Math.max(0,g[i]-m).toFixed(3)} g`:`needs ${m.toFixed(3)} g`):`${m.toFixed(3)} g ${g[i]?'(target)':'theoretical'}`});
  $('#stoichNote').textContent=ext===null?'':(lim<n?`Limiting reagent: ${sp[lim].raw} · reaction extent ${ext.toPrecision(4)} mol`:`Masses computed from the target amount of ${sp[lim].raw}.`)};
 box.querySelectorAll('input').forEach(i=>i.oninput=calc)}
function equationHTML(b){const n=b.L.length;const left=b.L.map((c,i)=>(b.co[i]===1n?'':`<span class="co">${b.co[i]}</span> `)+formulaHTML(c.raw)).join(' + ');const right=b.R.map((c,i)=>{const j=n+i;return (b.co[j]===1n?'':`<span class="co">${b.co[j]}</span> `)+formulaHTML(c.raw)}).join(' + ');return left+` <span class="arrow">→</span> `+right}
function stepsHTML(b,checks){const n=b.L.length,rows=b.elements.map(e=>{const L=b.L.map((c,i)=>`${b.co[i]}×${c.atoms[e]||0}`).join(' + '),R=b.R.map((c,i)=>`${b.co[n+i]}×${c.atoms[e]||0}`).join(' + ');const total=checks.find(x=>x.e===e);return `<div class="step"><div class="num">•</div><div><b>${e} conservation</b><span class="muted">${L} = ${R} → ${total.l.toString()} = ${total.r.toString()}</span></div></div>`}).join('');return `<div class="steps"><div class="step"><div class="num">1</div><div><b>Parse the reaction</b><span class="muted">Each compound is converted into element counts. Parentheses, brackets and hydrate groups are expanded.</span></div></div><div class="step"><div class="num">2</div><div><b>Build conservation equations</b><span class="muted">Every element must have the same total count on both sides. Ionic charge is also conserved when supplied.</span></div></div><div class="step"><div class="num">3</div><div><b>Solve for coefficient ratios</b><span class="muted">StoichLab solves the homogeneous linear system using exact rational arithmetic.</span></div></div><div class="step"><div class="num">4</div><div><b>Reduce to the smallest integers</b><span class="muted">The coefficient vector is converted to whole numbers and divided by its greatest common divisor.</span></div></div>${rows}</div>`}
const FIXED={F:-1,O:-2,H:1,Li:1,Na:1,K:1,Rb:1,Cs:1,Be:2,Mg:2,Ca:2,Sr:2,Ba:2,Al:3,Zn:2,Ag:1};
function oxState(c){const els=Object.keys(c.atoms);if(!els.length)return null;if(els.length===1)return {[els[0]]:c.charge/c.atoms[els[0]]};
 const hasOF=c.atoms.O||c.atoms.F,st={},unk=[];let sum=0;
 for(const e of els){let v=FIXED[e];if(v===undefined&&['Cl','Br','I'].includes(e)&&!hasOF)v=-1;if(v===undefined)unk.push(e);else{st[e]=v;sum+=v*c.atoms[e]}}
 if(unk.length>1)return null;if(unk.length===1)st[unk[0]]=(c.charge-sum)/c.atoms[unk[0]];else if(sum!==c.charge)return null;return st}
function redoxRows(b){const side=sp=>{const m={};sp.forEach(c=>{const s=oxState(c);if(!s)return;for(const e in s)(m[e]=m[e]||new Set()).add(s[e])});return m},l=side(b.L),r=side(b.R),rows=[];
 for(const e in l)if(r[e]){const a=[...l[e]],z=[...r[e]];if(a.length===z.length&&a.every(x=>z.some(y=>Math.abs(x-y)<1e-9)))continue;rows.push({e,a,z,k:z.every(y=>a.every(x=>y>x))?'oxidized':z.every(y=>a.every(x=>y<x))?'reduced':'disproportionation'})}return rows}
function redoxHTML(b){const rows=redoxRows(b),f=v=>(v>0?'+':v<0?'−':'')+Math.abs(Number.isInteger(v)?v:+v.toFixed(2)),agents=k=>[...new Set(rows.filter(x=>x.k===k).flatMap(x=>b.L.filter(c=>c.atoms[x.e]).map(c=>c.raw)))].map(formulaHTML).join(', ');
 if(!rows.length)return '<h3 style="margin-top:16px">Redox analysis</h3><p class="muted">No oxidation-state change detected among resolvable species — likely not a redox reaction (e.g. acid–base or precipitation).</p>';
 return `<h3 style="margin-top:16px">Redox analysis</h3><div class="chips">${rows.map(x=>`<span class="chip">${x.e}: ${x.a.map(f).join(', ')} → ${x.z.map(f).join(', ')} · ${x.k}</span>`).join('')}</div><p class="muted" style="margin-top:10px">${agents('oxidized')?'Reducing agent: '+agents('oxidized')+'. ':''}${agents('reduced')?'Oxidizing agent: '+agents('reduced')+'.':''} Oxidation states use standard textbook rules and are skipped where ambiguous (peroxides, hydrides, mixed anions).</p>`}
