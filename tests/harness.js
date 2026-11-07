// Loads SangerScope's <script> (index.html) into a Node vm context with a synthetic gene and a mocked NCBI.
// No network: VariantValidator / annotation calls are stubbed out.
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(require('path').join(__dirname, '..', 'index.html'), 'utf8');
let src = html.slice(html.lastIndexOf('<script>') + 8, html.lastIndexOf('</script>'));
src = src.replace(/\nbind\(\);\s*$/, '\n');
const el = () => new Proxy(function () {}, { get: (t, k) => k === Symbol.toPrimitive ? () => '' : el(), set: () => true, apply: () => el() });
const ctx = { console, setTimeout, clearTimeout, Math, JSON, Date, Promise, Map, Set, Int32Array, Uint8Array, Uint16Array, Int16Array, Int8Array, Float32Array, DataView, ArrayBuffer, TextDecoder,
  document: el(), localStorage: { getItem: () => null, setItem() {}, removeItem() {}, key() {}, length: 0 }, fetch: () => Promise.reject(new Error('offline')),
  navigator: {}, requestAnimationFrame: () => 0, window: el() };
vm.createContext(ctx);
vm.runInContext(src, ctx, { filename: 'index.js' });
vm.runInContext(`flushNom = async () => {}; annotate = () => {}; ST.auto = false;`, ctx);

const COMPL = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const rc = s => [...s].reverse().map(c => COMPL[c] || 'N').join('');

// Build a gene in transcript orientation: [flank][ex1][intron]...[exN][flank]
// spec.parts: array of strings alternating exon,intron,exon... ; cds given as mRNA 1-based [cdsS,cdsE]
function loadGene({ flank5, flank3, exons, introns, cdsS, cdsE, strand = 1, chrStart = 1000000 }) {
  const FL = 250;
  if (flank5.length !== FL || flank3.length !== FL) throw new Error('flanks must be 250');
  let T = flank5; const exT = []; // transcript-orientation exon coords (0-based in T)
  exons.forEach((e, k) => { exT.push([T.length, T.length + e.length - 1]); T += e; if (k < introns.length) T += introns[k]; });
  T += flank3;
  const L = T.length, start = chrStart, stop = chrStart + L - 1;   // gene+flanks occupies chr start..stop
  const plus = strand > 0 ? T : rc(T);
  // chromosome coords of exons (begin<=end) in transcript order
  const gex = exT.map(([a, b], k) => strand > 0 ? { begin: start + a, end: start + b, order: k + 1 } : { begin: stop - b, end: stop - a, order: k + 1 });
  const mrna = exons.join('');
  const cds = mrna.slice(cdsS - 1, cdsE);
  let prot = ''; for (let i = 0; i + 2 < cds.length; i += 3) prot += vm.runInContext(`tr1(${JSON.stringify(cds.slice(i, i + 3))})`, ctx);
  const gene = { product: { symbol: 'SIM', description: 'simulated', gene_id: '1', transcripts: [{ accession_version: 'NM_000001.1', select_category: 'MANE_SELECT', cds: { range: [{ begin: String(cdsS), end: String(cdsE) }] }, protein: { accession_version: 'NP_000001.1' },
    genomic_locations: [{ genomic_accession_version: 'NC_000007.14', sequence_name: 'Primary Assembly GRCh38', genomic_range: { orientation: strand > 0 ? 'plus' : 'minus' }, exons: gex }] }] } };
  ctx.__mock = { plus, chrStart: start, mrna, trans: prot.replace(/\*$/, '') };
  vm.runInContext(`ncbi = async (url, type) => {
    if (/efetch.*rettype=fasta/.test(url)) { const m = url.match(/seq_start=(\\d+)&seq_stop=(\\d+)&strand=(\\d)/); const a = +m[1] - __mock.chrStart, b = +m[2] - __mock.chrStart; let s = __mock.plus.slice(a, b + 1); if (m[3] === '2') s = revcomp(s); return '>x\\n' + s; }
    if (/rettype=gb/.test(url)) return 'LOCUS x\\n  /translation="' + __mock.trans + '"\\nORIGIN\\n' + __mock.mrna.toLowerCase() + '\\n//';
    throw new Error('unexpected ' + url);
  };`, ctx);
  ctx.__prod = gene.product; ctx.__prod._txs = gene.product.transcripts;
  return vm.runInContext(`PLUS = null; PROD = __prod; loadTranscript(0).then(() => { renderGene = () => {}; return R; })`, ctx).then(() => ({ T, plus, L, exT, mrna, cds, prot, start, stop, strand }));
}
// stub UI renderers
vm.runInContext(`for (const f of ['renderGene','setGeneState','loadScope','refresh','renderRef','refreshVars','renderVariants','renderSummary','renderCoverage','renderAlignment','renderFlags','renderReads','toast']) globalThis[f] = () => {};`, ctx);

function run(code) { return vm.runInContext(code, ctx); }
function setReads(reads) {
  ctx.__reads = reads;
  run(`READS.length = 0; for (const r of __reads) addRead(r); computeVariants();`);
  return run(`({vars: VARS.map(v => ({hgvs: v.hgvs, kind: v.kind, s: v.s, e: v.e, ref: v.ref, alt: v.alt, zyg: v.zyg, p: v.p, type: v.type, loc: v.loc, effect: v.effect, flags: v.flags.map(f => f.t), nsup: v.sup.length, g: (()=>{try{return genomicOf(v)}catch(e){return {err: e.message}}})() })),
    reads: READS.map(r => ({name: r.name, status: r.res && r.res.status, ident: r.res && r.res.ident, rev: r.res && r.res.rev, mixed: r.res && r.res.mixInfo, het: r.res && r.res.het && {shift: r.res.het.shift, hgvs: r.res.het.hgvs}}))})`);
}
module.exports = { ctx, run, loadGene, setReads, rc };
