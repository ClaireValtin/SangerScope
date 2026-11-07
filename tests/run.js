#!/usr/bin/env node
/* SangerScope variant-calling regression tests — run with:  node tests/run.js  [--quick]
 *
 * 1. Homozygous: random SNV / del / ins / dup / delins / MNV / repeat changes anywhere in a synthetic 5-exon gene
 *    (both strands), read as forward or reverse FASTA; c., p., g. and VCF are compared with an independent implementation.
 * 2. Heterozygous: the same kinds of change as simulated two-allele AB1 traces (50 % and 30 % allele 2).
 * 3. Scenarios: hand-built cases for bugs fixed in the past.
 * Exits non-zero if anything regresses. */
const H = require('./harness'), T_ = require('./truth'), { makeAB } = require('./trace'); const { rc } = H;
const QUICK = process.argv.includes('--quick');
let failures = 0; const fail = m => { failures++; console.log('  FAIL ' + m); };

function expectG(P, P2, start) {
  let p = 0; while (p < P.length && P[p] === P2[p]) p++;
  let s = 0; while (s < P.length - p && s < P2.length - p && P[P.length - 1 - s] === P2[P2.length - 1 - s]) s++;
  const ref = P.slice(p, P.length - s), alt = P2.slice(p, P2.length - s), a = p, b = P.length - s - 1, g = i => start + i;
  let hg;
  if (!ref) { const k = alt.length; hg = P.slice(a - k, a) === alt ? (k === 1 ? `g.${g(a - 1)}dup` : `g.${g(a - k)}_${g(a - 1)}dup`) : `g.${g(a - 1)}_${g(a)}ins${alt}`; }
  else if (!alt) hg = a === b ? `g.${g(a)}del` : `g.${g(a)}_${g(b)}del`;
  else if (ref.length === 1 && alt.length === 1) hg = `g.${g(a)}${ref}>${alt}`;
  else hg = (a === b ? `g.${g(a)}` : `g.${g(a)}_${g(b)}`) + `delins${alt}`;
  let s2 = 0; while (s2 < P.length && s2 < P2.length && P[P.length - 1 - s2] === P2[P2.length - 1 - s2]) s2++;
  let p2 = 0; while (p2 < P.length - s2 && p2 < P2.length - s2 && P[p2] === P2[p2]) p2++;
  let vr = P.slice(p2, P.length - s2), va = P2.slice(p2, P2.length - s2), vpos = p2;
  if (!vr || !va) { vpos--; vr = P[vpos] + vr; va = P[vpos] + va; }
  return { hg: `NC_000007.14:${hg}`, vcf: { pos: g(vpos), ref: vr, alt: va } };
}
const REPS = ['AAAAAAAA', 'GAGGAGGAGGAG', 'CTGCTGCTG', 'CAGCAGCAGCAGCAG', 'TTTTTTT', 'CACACACACA'];
function randomChange(T, i, kinds) {
  const kind = T_.pick(kinds);
  if (kind === 'snv') return T.slice(0, i) + T_.pick([...'ACGT'].filter(x => x !== T[i])) + T.slice(i + 1);
  if (kind === 'del') { const k = 1 + T_.ri(T_.pick([3, 6, 12, 25])); return T.slice(0, i) + T.slice(i + k); }
  if (kind === 'ins') { const k = 1 + T_.ri(T_.pick([3, 6, 12, 25])); return T.slice(0, i) + T_.randSeq(k) + T.slice(i); }
  if (kind === 'dup') { const k = 1 + T_.ri(T_.pick([3, 6, 12])); return T.slice(0, i + k) + T.slice(i, i + k) + T.slice(i + k); }
  if (kind === 'delins') { const k = 1 + T_.ri(6), m = 1 + T_.ri(6); return T.slice(0, i) + T_.randSeq(m) + T.slice(i + k); }
  if (kind === 'mnv') { const k = 2 + T_.ri(3); let x = ''; for (let q = 0; q < k; q++) x += T_.pick([...'ACGT'].filter(y => y !== T[i + q])); return T.slice(0, i) + x + T.slice(i + k); }
  const r = T_.pick(REPS), hits = []; let q = -1; while ((q = T.indexOf(r, q + 1)) >= 0) if (q > 300 && q < T.length - 300) hits.push(q);
  if (!hits.length) return null;
  const at = T_.pick(hits), unit = /^[AT]/.test(r) ? r[0] : r.startsWith('CA') && r[2] === 'C' ? 'CA' : r.slice(0, 3), u = unit.repeat(1 + T_.ri(2));
  return kind === 'repdel' ? T.slice(0, at) + T.slice(at + u.length) : T.slice(0, at) + u + T.slice(at);
}

async function homozygous(seed, strand, N) {
  const spec = T_.makeGene(seed), G = await H.loadGene({ ...spec, strand }), C = T_.makeCoord(spec, G.exT), T = G.T;
  if (!H.run('R.gates.every(g => g.ok)')) fail(`seed ${seed}: reference gates failed`);
  let n = 0, split = 0;
  for (let t = 0; t < N; t++) {
    const r = T_.ri(10), ex = r < 6 ? T_.pick(G.exT.slice(1)) : null;
    const i = ex ? ex[0] + T_.ri(ex[1] - ex[0] + 1) : r < 8 ? T_.pick(T_.pick(G.exT)) - 6 + T_.ri(13) : 260 + T_.ri(T.length - 520);
    const T2 = randomChange(T, i, ['snv', 'snv', 'del', 'del', 'ins', 'dup', 'delins', 'mnv', 'repdel', 'repins']); if (!T2) continue;
    const exp = T_.expectC(T, T2, C); if (!exp.hgvs) continue;
    let read = T2.slice(Math.max(0, exp.s - 260), exp.s + 260); const rev = T_.rand() < 0.5; if (rev) read = rc(read);
    const out = H.setReads([{ name: 'r', seq: read }]); n++;
    const tag = `hom seed ${seed} strand ${strand} ${rev ? 'R' : 'F'}: ${exp.hgvs}`;
    if (out.reads[0].status !== 'ok') { fail(`${tag}: read ${out.reads[0].status}`); continue; }
    const v = out.vars.find(x => x.hgvs === exp.hgvs);
    if (!v || out.vars.length !== 1) {
      // a complex delins may be split by the aligner into pieces ≥ 2 nt apart (valid HGVS) — then the pieces must carry the in-cis note
      if (exp.kind === 'delins' && out.vars.length > 1 && out.vars.every(x => x.flags.some(f => f.startsWith('in cis with')))) { split++; continue; }
      fail(`${tag}: got [${out.vars.map(x => x.hgvs).join(', ')}]`); continue;
    }
    const inEx = j => C.mpos(j) > 0, cds = j => C.mpos(j) >= spec.cdsS && C.mpos(j) <= spec.cdsE;
    const exonic = exp.kind === 'ins' ? inEx(exp.s) && inEx(exp.e) && C.mpos(exp.e) === C.mpos(exp.s) + 1 : inEx(exp.s) && inEx(exp.e) && C.mpos(exp.e) - C.mpos(exp.s) === exp.e - exp.s;
    const m1 = exp.kind === 'ins' ? C.mpos(exp.e) : exp.kind === 'dup' ? C.mpos(exp.e) + 1 : C.mpos(exp.s), m2 = exp.kind === 'ins' || exp.kind === 'dup' ? m1 - 1 : C.mpos(exp.e);
    if (exonic && [exp.s, exp.e].some(cds) && m1 <= spec.cdsE) {
      const mrna2 = G.mrna.slice(0, m1 - 1) + (exp.kind === 'del' ? '' : exp.alt) + G.mrna.slice(m2);
      const ep = T_.expectP(G.mrna, mrna2, spec.cdsS);
      if (!(v.p === ep || (ep === 'p.(=)' && /^p\.\([A-Z][a-z]{2}\d+=\)$/.test(v.p)))) fail(`${tag}: p. expected ${ep}, got ${v.p} (${v.type})`);
    }
    const eg = expectG(G.plus, strand > 0 ? T2 : rc(T2), G.start);
    if (v.g.g !== eg.hg) fail(`${tag}: g. expected ${eg.hg}, got ${v.g.g}`);
    if (v.g.vcf.pos !== eg.vcf.pos || v.g.vcf.ref !== eg.vcf.ref || v.g.vcf.alt !== eg.vcf.alt) fail(`${tag}: VCF expected ${eg.vcf.pos} ${eg.vcf.ref}>${eg.vcf.alt}, got ${v.g.vcf.pos} ${v.g.vcf.ref}>${v.g.vcf.alt}`);
  }
  console.log(`  homozygous  seed ${seed} strand ${strand > 0 ? '+' : '−'}: ${n} variants (${split} complex delins split with in-cis note)`);
}

async function heterozygous(seed, strand, N, frac2, minRate) {
  const spec = T_.makeGene(seed), G = await H.loadGene({ ...spec, strand }), C = T_.makeCoord(spec, G.exT), T = G.T;
  const st = { snv: [0, 0], indel: [0, 0] }; const miss = [];
  for (let t = 0; t < N; t++) {
    const T2 = randomChange(T, 400 + T_.ri(T.length - 800), ['snv', 'snv', 'del', 'ins', 'dup', 'repdel', 'repins', 'delins']); if (!T2) continue;
    const exp = T_.expectC(T, T2, C); if (!exp.hgvs || exp.kind === 'delins' && T2.length === T.length) continue;   // equal-length het MNV: phase unknown, reported as SNVs
    const rev = T_.rand() < 0.5, NL = 480, before = T_.pick([60, 120, 200]); let a1, a2;
    if (!rev) { const c0 = exp.s - before; if (c0 < 0) continue; a1 = T.slice(c0, c0 + NL); a2 = T2.slice(c0, c0 + NL); }
    else { const c1 = exp.e + before, d = T2.length - T.length; if (c1 - NL < 0) continue; a1 = rc(T.slice(c1 - NL, c1)); a2 = rc(T2.slice(c1 + d - NL, c1 + d)); }
    const out = H.setReads([{ name: 'r', ab: makeAB(a1, a2, { frac2 }) }]);
    const k = exp.kind === 'snv' ? 'snv' : 'indel'; st[k][1]++;
    const v = out.vars.find(x => x.hgvs === exp.hgvs);
    if (v && v.zyg === 'het' && out.vars.length === 1) st[k][0]++; else miss.push(`${exp.hgvs} ${rev ? 'R' : 'F'} → [${out.vars.map(x => x.hgvs).join(', ')}]`);
  }
  const rate = k => st[k][1] ? st[k][0] / st[k][1] : 1;
  console.log(`  heterozygous seed ${seed} strand ${strand > 0 ? '+' : '−'} allele-2 ${frac2 * 100}%: SNV ${st.snv[0]}/${st.snv[1]}, indel ${st.indel[0]}/${st.indel[1]}`);
  for (const k of ['snv', 'indel']) if (rate(k) < minRate[k]) fail(`het ${k} rate ${(rate(k) * 100).toFixed(1)}% < ${minRate[k] * 100}% (seed ${seed}): ${miss.slice(0, 5).join('; ')}`);
}

async function scenarios() {
  const spec = T_.makeGene(777); await H.loadGene({ ...spec, strand: 1 }); const T = H.run('R.seq');
  const idx = c => H.run(`cToIndex(${JSON.stringify(c)})`), sub = (S, i, x) => S.slice(0, i) + x + S.slice(i + x.length);
  const fa = (S, c) => S.slice(c - 250, c + 250);
  const het = (A1, A2, c, before = 120) => makeAB(A1.slice(c - before, c - before + 480), A2.slice(c - before, c - before + 480));
  const one = (name, reads, check) => { const o = H.setReads(reads); const r = check(o); if (r !== true) fail(`scenario "${name}": ${r} — got [${o.vars.map(v => `${v.hgvs} ${v.p}`).join(', ')}]`); };
  const has = (o, h, p) => o.vars.some(v => v.hgvs === h && (p == null || v.p === p));
  let x = idx('c.31');   // codon 11 GTT (Val) → TAT (Tyr)
  one('two substitutions in one codon → one delins, joint protein effect', [{ name: 'F', seq: fa(sub(T, x, T[x] === 'T' ? 'GA' : 'TA'), x) }], o => o.vars.length === 1 && o.vars[0].kind === 'delins' && /Val11/.test(o.vars[0].p) || 'not merged');
  x = idx('c.150');
  one('extra N in an insertion is not turned into a base', [{ name: 'F', seq: fa(T.slice(0, x + 1) + 'N' + T.slice(x + 1), x) }], o => o.vars.length === 0 || 'variant invented from N');
  one('in-frame insertion with a stop is stop gained', [{ name: 'F', seq: fa(T.slice(0, x + 1) + 'GCCTAA' + T.slice(x + 1), x) }], o => o.vars.length === 1 && o.vars[0].type === 'stop gained' && /Ter\)$/.test(o.vars[0].p) || 'not stop gained');
  x = idx('c.70');
  one('dup of an exon\'s last base gets a protein effect', [{ name: 'F', seq: fa(T.slice(0, x + 1) + T[x] + T.slice(x + 1), x) }], o => o.vars.length === 1 && o.vars[0].hgvs === 'c.70dup' && o.vars[0].type === 'frameshift' || 'no frameshift');
  x = idx('c.36');
  one('het delins resolved by deconvolution', [{ name: 'F', ab: het(T, T.slice(0, x) + 'TG' + T.slice(x + 3), x) }], o => has(o, 'c.36_38delinsTG') || 'wrong het allele');
  x = idx('c.40');
  one('unresolvable compound het indel is reported, not dropped', [{ name: 'F', ab: het(T.slice(0, x) + T.slice(x + 1), T.slice(0, x) + T.slice(x + 4), x) }], o => o.reads[0].mixed.some(m => m.beyondTrim) || 'no mixed-trace warning');
  const s2 = idx('c.120'), A1 = sub(T, s2, T[s2] === 'A' ? 'G' : 'A');
  one('het indel F + R; hom SNV downstream found by the opposite read', [{ name: 'F', ab: het(A1, A1.slice(0, x) + A1.slice(x + 2), x, 100) }, { name: 'R', ab: makeAB(rc(A1.slice(s2 + 150 - 480, s2 + 150)), rc(A1.slice(0, x).concat(A1.slice(x + 2)).slice(s2 + 148 - 480, s2 + 148))) }],
    o => has(o, 'c.40_41del') && has(o, `c.120A>G`) || 'missing');
  one('start lost', [{ name: 'F', seq: fa(sub(T, idx('c.1'), 'G'), idx('c.1')) }], o => o.vars[0]?.p === 'p.(Met1?)' || 'not start lost');
  x = idx('c.424');
  one('stop lost', [{ name: 'F', seq: fa(sub(T, x, 'C'), x) }], o => o.vars[0]?.type === 'stop lost' || 'not stop lost');
  console.log('  scenarios done');
}

(async () => {
  const N = QUICK ? 300 : 1500, NH = QUICK ? 80 : 200;
  console.log('SangerScope simulation tests');
  for (const seed of QUICK ? [777] : [777, 11]) for (const strand of [1, -1]) await homozygous(seed, strand, N);
  for (const seed of QUICK ? [777] : [777, 11, 33]) for (const strand of [1, -1]) for (const f of [0.5, 0.3]) await heterozygous(seed, strand, NH, f, { snv: 1, indel: 0.95 });
  await scenarios();
  console.log(failures ? `\n${failures} failure(s)` : '\nall passed');
  process.exit(failures ? 1 : 0);
})();
