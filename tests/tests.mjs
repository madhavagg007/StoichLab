// Run from repo root: node tests/tests.mjs  (tests the pure chemistry engine in js/engine.js)
import fs from 'fs';
const html=fs.readFileSync(new URL('../js/engine.js',import.meta.url),'utf8');
const code=html.replace("'use strict';",'');
const {balanceEquation,verify,parseFormula,molar,oxState,redoxRows}=new Function(code+';return {balanceEquation,verify,parseFormula,molar,oxState,redoxRows}')();
let pass=0,fail=0;const t=(n,f)=>{try{f();pass++;console.log('✓',n)}catch(e){fail++;console.log('✗',n,'-',e.message)}};
const eq=(a,b)=>{if(JSON.stringify(a)!==JSON.stringify(b))throw Error(JSON.stringify(a)+' != '+JSON.stringify(b))};
const co=s=>balanceEquation(s).co.map(Number);
[['H2 + O2 -> H2O',[2,1,2]],['C3H8 + O2 -> CO2 + H2O',[1,5,3,4]],['Fe + O2 -> Fe2O3',[4,3,2]],['Ca(OH)2 + HCl -> CaCl2 + H2O',[1,2,1,2]],['Al + HCl -> AlCl3 + H2',[2,6,2,3]],['KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2',[2,16,2,2,8,5]],['Fe^3+ + e- -> Fe^2+',[1,1,1]],['MnO4^- + Fe^2+ + H+ -> Mn^2+ + Fe^3+ + H2O',[1,5,8,1,5,4]],['CuSO4·5H2O -> CuSO4 + H2O',[1,1,5]]].forEach(([q,e])=>t('balance '+q,()=>eq(co(q),e)));
t('all balances verify',()=>{const b=balanceEquation('C3H8 + O2 -> CO2 + H2O');if(!verify(b).every(x=>x.ok))throw Error('unbalanced')});
t('nested brackets',()=>eq(parseFormula('Al2(SO4)3').atoms,{Al:2,S:3,O:12}));
t('Co vs CO',()=>{eq(parseFormula('Co').atoms,{Co:1});eq(parseFormula('CO').atoms,{C:1,O:1})});
t('molar mass H2O',()=>{const m=molar({atoms:{H:2,O:1}});if(Math.abs(m-18.015)>0.01)throw Error(m)});
t('oxidation state Mn in KMnO4',()=>eq(oxState({atoms:{K:1,Mn:1,O:4},charge:0}).Mn,7));
t('redox detected',()=>{if(redoxRows(balanceEquation('Fe + O2 -> Fe2O3')).length!==2)throw Error('expected Fe and O')});
t('acid-base not redox',()=>eq(redoxRows(balanceEquation('HCl + NaOH -> NaCl + H2O')).length,0));
['H2 + O2','Xx + O2 -> XxO','H2 + O2 -> -> H2O','Ca(OH + O2 -> CaO','H2 -> O2'].forEach(q=>t('rejects '+q,()=>{try{balanceEquation(q)}catch{return}throw Error('should throw')}));
console.log(`\n${pass} passed, ${fail} failed`);process.exit(fail?1:0)
