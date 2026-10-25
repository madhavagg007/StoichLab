'use strict';
const $=s=>document.querySelector(s), historyKey='chembalance-history-v2';
function reactionAnalysis(b){
 const left=b.L.map(x=>x.atoms), right=b.R.map(x=>x.atoms);
 const has=(side,el)=>side.some(a=>a[el]);
 const allL=left.reduce((o,a)=>{for(const k in a)o[k]=(o[k]||0)+a[k];return o},{}), allR=right.reduce((o,a)=>{for(const k in a)o[k]=(o[k]||0)+a[k];return o},{});
 let type='General reaction';
 if(b.L.length===2&&b.R.length===1) type='Combination / synthesis';
 else if(b.L.length===1&&b.R.length>=2) type='Decomposition';
 else if(has(left,'C')&&has(left,'H')&&has(right,'C')&&has(right,'O')&&has(left,'O')) type='Possible combustion';
 else if(b.L.length===2&&b.R.length===2) type='Double displacement / exchange';
 const charge=verify(b).find(x=>x.e==='Charge');
 const conserved=b.elements.filter(e=>e!=='Charge').length;
 return `<div class="grid"><div class="metric"><strong>${type}</strong><span>Heuristic reaction class</span></div><div class="metric"><strong>${conserved}</strong><span>Conserved elements</span></div><div class="metric"><strong>${charge? (charge.ok?'Yes':'No'):'N/A'}</strong><span>Charge conserved</span></div><div class="metric"><strong>${b.co.reduce((a,x)=>a+Number(x),0)}</strong><span>Total coefficient sum</span></div></div><p class="muted" style="margin-top:12px">Reaction classification is heuristic; the balancing result itself is determined from conservation constraints.</p>`;
}
function saveHistory(input,display){let h=JSON.parse(localStorage.getItem(historyKey)||'[]').filter(x=>x.input!==input);h.unshift({input,display,time:new Date().toLocaleString()});localStorage.setItem(historyKey,JSON.stringify(h.slice(0,12)));renderHistory()}
function renderHistory(){const h=JSON.parse(localStorage.getItem(historyKey)||'[]'),box=$('#history');if(!h.length){box.innerHTML='<p class="muted">No saved equations yet.</p>';return}box.innerHTML=h.map((x,i)=>`<button data-hist="${i}">${esc(x.display)}<small>${esc(x.time)}</small></button>`).join('');box.querySelectorAll('[data-hist]').forEach(btn=>btn.onclick=()=>{const x=h[Number(btn.dataset.hist)];$('#eq').value=x.input;run()})}
function toast(msg){const t=$('#toast');t.textContent=msg;t.classList.add('show');clearTimeout(window.__toast);window.__toast=setTimeout(()=>t.classList.remove('show'),1800)}
function renderError(err){$('#result').innerHTML=`<div class="resultHead"><h2>Could not balance</h2><span class="status bad">CHECK INPUT</span></div><p class="error" style="color:var(--red);font-weight:800">${esc(err.message)}</p><p class="muted">Try one of the examples above, or check element capitalization, brackets, coefficients, and the reaction arrow.</p>`;$('#analysis').classList.add('hidden');$('#explain').classList.add('hidden');$('#verification').classList.add('hidden');$('#stoich').classList.add('hidden')}
function run(){const raw=$('#eq').value.trim();if(!raw){renderError(new Error('Enter a chemical equation first.'));return}try{const b=balanceEquation(raw),checks=verify(b);if(!checks.every(x=>x.ok))throw Error('Internal verification failed.');const display=equationHTML(b);$('#result').innerHTML=`<div class="resultHead"><div><span class="status ok">BALANCED & VERIFIED</span></div><div class="toolbar"><button class="btn secondary" id="copyBtn">Copy</button><button class="btn secondary" id="shareBtn">Share</button></div></div><div class="equation" aria-label="Balanced equation">${display}</div>${checkMsg(b)}${b.multi?'<p class="muted"><span class="status warn">MULTIPLE SOLUTIONS</span> This reaction can be balanced in several independent ways; the smallest all-positive set is shown.</p>':''}<div class="grid"><div class="metric"><strong>${b.co.length}</strong><span>Species</span></div><div class="metric"><strong>${b.elements.filter(e=>e!=='Charge').length}</strong><span>Elements conserved</span></div><div class="metric"><strong>${Math.max(...b.co.map(x=>Number(x)))}</strong><span>Largest coefficient</span></div><div class="metric"><strong>${b.co.every(x=>x===1n)?'1:1':'Reduced'}</strong><span>Coefficient state</span></div></div>`;
 $('#analysis').classList.remove('hidden');$('#analysisBody').innerHTML=reactionAnalysis(b)+redoxHTML(b);
 $('#explain').classList.remove('hidden');$('#explain').innerHTML=`<h2>How StoichLab solved it</h2><p class="muted">The result below is derived from conservation constraints rather than a memorized reaction.</p>${stepsHTML(b,checks)}`;
 $('#verification').classList.remove('hidden');$('#verification').innerHTML=`<h2>Conservation audit</h2><div class="scroll"><table><thead><tr><th>Quantity</th><th>Left</th><th>Right</th><th>Result</th></tr></thead><tbody>${checks.map(x=>`<tr><td>${x.e}</td><td>${x.l}</td><td>${x.r}</td><td class="check">✓ conserved</td></tr>`).join('')}</tbody></table></div>`;
 renderStoich(b);
 $('#copyBtn').onclick=()=>{const plain=b.L.map((c,i)=>(b.co[i]===1n?'':b.co[i]+' ')+c.raw).join(' + ')+' -> '+b.R.map((c,i)=>(b.co[b.L.length+i]===1n?'':b.co[b.L.length+i]+' ')+c.raw).join(' + ');navigator.clipboard?.writeText(plain).then(()=>toast('Balanced equation copied')).catch(()=>toast('Copy unavailable in this browser'))};
 $('#shareBtn').onclick=async()=>{const text='StoichLab: '+display.replace(/<[^>]+>/g,'');if(navigator.share)try{await navigator.share({title:'StoichLab Lab',text});}catch{}else{navigator.clipboard?.writeText(location.href.split('#')[0]+'#'+encodeURIComponent(raw));toast('Share link copied')}};
 saveHistory(raw,display.replace(/<[^>]+>/g,''));document.title='StoichLab · '+display.replace(/<[^>]+>/g,'');
 }catch(err){renderError(err)} }
$('#balanceBtn').onclick=run;
$('#eq').addEventListener('keydown',e=>{if((e.ctrlKey||e.metaKey)&&e.key==='Enter'){e.preventDefault();run()}});
$('#clearInput').onclick=()=>{$('#eq').value='';$('#eq').focus();toast('Equation input cleared')};
$('#insertArrow').onclick=()=>{const el=$('#eq'),start=el.selectionStart,end=el.selectionEnd;el.setRangeText(' → ',start,end,'end');el.focus()};
$('#ptable').innerHTML=[...ELEMENTS].map(e=>`<button type="button" class="btn secondary" title="${e}${MASS[e]?' · '+MASS[e]+' g/mol':''}">${e}</button>`).join('');$('#ptable').onclick=ev=>{const t=ev.target.closest('button');if(!t)return;const el=$('#eq');el.setRangeText(t.textContent,el.selectionStart,el.selectionEnd,'end');el.focus()};
$('#loadTemplate').onclick=()=>{$('#eq').value='Reactant1 + Reactant2 -> Product1 + Product2';$('#eq').focus();$('#eq').select();toast('Template loaded — replace the placeholders')};
document.querySelectorAll('[data-example]').forEach(b=>b.onclick=()=>{$('#eq').value=b.dataset.example;$('#eq').focus();run()});
$('#clearHistory').onclick=()=>{localStorage.removeItem(historyKey);renderHistory();toast('History cleared')};
const themeKey='chembalance-theme-v3'; const savedTheme=localStorage.getItem(themeKey); if(savedTheme==='light')document.body.classList.add('light'); $('#themeBtn').onclick=()=>{document.body.classList.toggle('light');localStorage.setItem(themeKey,document.body.classList.contains('light')?'light':'dark')}; const VT={balance:'Balance',stoich:'Stoichiometry',redox:'Redox analysis',library:'Reaction library',history:'History',about:'About & settings'};
function show(v){document.querySelectorAll('.view').forEach(x=>x.classList.toggle('on',x.id==='v-'+v));document.querySelectorAll('.nb').forEach(x=>x.classList.toggle('on',x.dataset.v===v));$('#vtitle').textContent=VT[v];scrollTo(0,0)}
document.querySelector('.nav').onclick=e=>{const b=e.target.closest('.nb');if(b)show(b.dataset.v)};
const mo=new MutationObserver(()=>{[['stoich','stoich'],['analysis','redox']].forEach(([id,v])=>document.querySelector(`.nb[data-v="${v}"]`).classList.toggle('has',!$('#'+id).classList.contains('hidden')))});
['stoich','analysis'].forEach(id=>mo.observe($('#'+id),{attributes:true,attributeFilter:['class']}));

let libCat='All';
function renderLib(){const q=$('#libQ').value.toLowerCase(),cats=['All',...new Set(LIB.map(x=>x[0]))];
 $('#libCats').innerHTML=cats.map(c=>`<span class="chip${c===libCat?' on':''}" data-c="${esc(c)}" role="button" tabindex="0">${esc(c)}</span>`).join('');
 $('#lib').innerHTML=LIB.filter(x=>(libCat==='All'||x[0]===libCat)&&(x[1]+x[2]+x[0]).toLowerCase().includes(q)).map(x=>`<button class="lc" type="button" data-eq="${esc(x[2])}"><small>${esc(x[0]).toUpperCase()}</small><b>${esc(x[1])}</b><code>${esc(x[2])}</code></button>`).join('')||'<div class="empty">No reactions match.</div>'}
$('#libQ').oninput=renderLib;$('#libCats').onclick=e=>{const c=e.target.closest('.chip');if(c){libCat=c.dataset.c;renderLib()}};
$('#lib').onclick=e=>{const b=e.target.closest('.lc');if(b){$('#eq').value=b.dataset.eq;show('balance');$('#balanceBtn').click()}};
$('#sEl').textContent=ELEMENTS.size;$('#sLib').textContent=LIB.length;renderLib();
const cnt=()=>{$('#sHist').textContent=$('#history').querySelectorAll('button').length};new MutationObserver(cnt).observe($('#history'),{childList:true,subtree:true});cnt();
$('#history').addEventListener('click',e=>{if(e.target!==$('#history'))show('balance')});
$('#aboutClear').onclick=()=>$('#clearHistory').click();
addEventListener('keydown',e=>{if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='k'){e.preventDefault();show('library');$('#libQ').focus()}});
let dip;addEventListener('beforeinstallprompt',e=>{e.preventDefault();dip=e;$('#installBtn').hidden=false});$('#installBtn').onclick=()=>dip&&dip.prompt();
if('serviceWorker' in navigator&&/^https?:/.test(location.protocol))navigator.serviceWorker.register('sw.js').catch(()=>{});if(location.hash.length>1){try{$('#eq').value=decodeURIComponent(location.hash.slice(1))}catch{}}renderHistory();run();
