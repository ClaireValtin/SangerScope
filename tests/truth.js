// Independent reference implementation of the expected answers (c., p., synthetic gene model).
const CODON = {}; { const b = 'TCAG', a = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; let n = 0; for (const x of b) for (const y of b) for (const z of b) CODON[x+y+z] = a[n++]; }
const AA3 = {A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',Q:'Gln',E:'Glu',G:'Gly',H:'His',I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',T:'Thr',W:'Trp',Y:'Tyr',V:'Val','*':'Ter'};
let rng = 12345; const rand = () => { rng = (rng * 1103515245 + 12345) & 0x7fffffff; return rng / 0x7fffffff; };
const seedRng = s => { rng = s; };
const ri = n => Math.floor(rand() * n), pick = a => a[ri(a.length)];
const randSeq = n => { let s = ''; for (let i = 0; i < n; i++) s += 'ACGT'[ri(4)]; return s; };
function randCodons(n) { const sense = Object.keys(CODON).filter(c => CODON[c] !== '*'); let s = ''; for (let i = 0; i < n; i++) s += pick(sense); return s; }

// gene model in transcript orientation
function makeGene(seed) {
  seedRng(seed);
  // intron: GT + random with a homopolymer & repeat + AG ; avoid accidental stops not relevant
  const intron = n => 'GTAAGT' + randSeq(20) + 'TTTTTTT' + randSeq(15) + 'CACACACACA' + randSeq(n) + 'TTTCCTTTTCAG';
  const cdsBody = [
    randCodons(30), 'AAAAAAAA' /* homopolymer in frame-agnostic */, randCodons(20), 'GAGGAGGAGGAG' /* Glu x4 */, randCodons(25),
    'CTGCTGCTG', randCodons(20), 'CAGCAGCAGCAGCAG', randCodons(30)];
  // make CDS: ATG + body (must be multiple of 3 overall) + stop
  let body = cdsBody.join(''); while (body.length % 3) body += 'C';
  // also guard no internal stop
  let cds = 'ATG' + body; for (let i = 0; i < cds.length; i += 3) if (CODON[cds.slice(i, i + 3)] === '*') cds = cds.slice(0, i) + 'CTG' + cds.slice(i + 3);
  cds += 'TGA';
  const utr5 = 'GCCGCC' + randSeq(30) + 'GCCACC', utr3 = randSeq(40) + 'AATAAA' + randSeq(30);
  // split CDS into 4 exons at different phases
  const cut = [70, 171, 302];   // phases 1, 0, 2 (70%3=1,171%3=0,302%3=2)
  const pieces = [cds.slice(0, cut[0]), cds.slice(cut[0], cut[1]), cds.slice(cut[1], cut[2]), cds.slice(cut[2])];
  const exons = [utr5 + pieces[0], pieces[1], pieces[2], pieces[3] + utr3];
  // first exon UTR-only noncoding exon prepended to test 5'UTR introns
  exons.unshift('ACGTTGCAGGCTAGCTAGGCT' + randSeq(20));
  const introns = [intron(60), intron(80), intron(70), intron(90)];
  const cdsS = exons[0].length + utr5.length + 1, cdsE = cdsS + cds.length - 1;
  return { flank5: randSeq(250), flank3: randSeq(250), exons, introns, cdsS, cdsE };
}

// coordinates in transcript-orientation T (0-based)
function makeCoord(g, exT) {
  const exs = exT.map(([s, e]) => ({ s, e })); let m = 1; for (const x of exs) { x.mS = m; m += x.e - x.s + 1; x.mE = m - 1; }
  const mrnaLen = m - 1;
  const mToC = mm => mm < g.cdsS ? String(mm - g.cdsS) : mm > g.cdsE ? '*' + (mm - g.cdsE) : String(mm - g.cdsS + 1);
  const mpos = i => { for (const x of exs) if (i >= x.s && i <= x.e) return x.mS + i - x.s; return 0; };
  const c = i => {
    if (i < exs[0].s) return mToC(1 - (exs[0].s - i));
    const last = exs[exs.length - 1]; if (i > last.e) return mToC(mrnaLen + i - last.e);
    const mm = mpos(i); if (mm) return mToC(mm);
    for (let k = 0; k < exs.length - 1; k++) if (i > exs[k].e && i < exs[k + 1].s) {
      const du = i - exs[k].e, dd = exs[k + 1].s - i; // intron length odd: middle gets +
      return du <= dd ? mToC(exs[k].mE) + '+' + du : mToC(exs[k + 1].mS) + '-' + dd;
    }
  };
  return { c, mpos, exs, mrnaLen };
}

// expected c. from ref T and mutated T2
function expectC(T, T2, C) {
  let p = 0; while (p < T.length && p < T2.length && T[p] === T2[p]) p++;
  let s = 0; while (s < T.length - p && s < T2.length - p && T[T.length - 1 - s] === T2[T2.length - 1 - s]) s++;
  const ref = T.slice(p, T.length - s), alt = T2.slice(p, T2.length - s), a = p, b = T.length - s - 1;
  const R = (x, y) => x === y ? `c.${C.c(x)}` : `c.${C.c(x)}_${C.c(y)}`;
  if (!ref && !alt) return { hgvs: null };
  if (!ref) { const k = alt.length; if (T.slice(a - k, a) === alt) return { hgvs: `${R(a - k, a - 1)}dup`, kind: 'dup', s: a - k, e: a - 1, alt };
    return { hgvs: `c.${C.c(a - 1)}_${C.c(a)}ins${alt}`, kind: 'ins', s: a - 1, e: a, alt }; }
  if (!alt) return { hgvs: `${R(a, b)}del`, kind: 'del', s: a, e: b, ref };
  if (ref.length === 1 && alt.length === 1) return { hgvs: `c.${C.c(a)}${ref}>${alt}`, kind: 'snv', s: a, e: a, ref, alt };
  // HGVS: a delins whose alt is the reverse complement would be inv; ignore
  return { hgvs: `${R(a, b)}delins${alt}`, kind: 'delins', s: a, e: b, ref, alt };
}
const tr = s => { let p = ''; for (let i = 0; i + 2 < s.length; i += 3) { const a = CODON[s.slice(i, i + 3)] || 'X'; p += a; if (a === '*') break; } return p; };
// expected p. from mRNAs
function expectP(mrna, mrna2, cdsS) {
  const P = tr(mrna.slice(cdsS - 1)), M = tr(mrna2.slice(cdsS - 1));
  const A = k => AA3[P[k]] + (k + 1), three = x => [...x].map(c => AA3[c] || 'Xaa').join('');
  if (mrna2.slice(cdsS - 1, cdsS + 2) !== 'ATG') return 'p.(Met1?)';
  if (P === M) { // synonymous
    return 'p.(=)';
  }
  let i = 0; while (P[i] === M[i]) i++;
  const fs = (mrna2.length - mrna.length) % 3 !== 0;
  const stop = M.indexOf('*') < 0 ? '?' : M.indexOf('*') - i + 1;
  if (M[i] === '*') return `p.(${A(i)}Ter)`;
  if (P[i] === '*') return `p.(Ter${i + 1}${AA3[M[i]]}extTer${stop})`;
  if (fs) return `p.(${A(i)}${AA3[M[i]]}fsTer${stop})`;
  { const d3 = (mrna2.length - mrna.length) / 3, st = M.indexOf('*'); if (st >= 0 && st < P.length - 1 + d3) return `p.(${A(i)}delins${three(M.slice(i, st + 1))})`; }
  // in-frame: prefix-first trimming gives the 3'-most description
  let s = 0; while (s < P.length - i && s < M.length - i && P[P.length - 1 - s] === M[M.length - 1 - s]) s++;
  const del = P.slice(i, P.length - s), ins = M.slice(i, M.length - s);
  if (!del && !ins) return 'p.(=)';
  if (!del) { const k = ins.length; if (P.slice(i - k, i) === ins) return k === 1 ? `p.(${A(i - 1)}dup)` : `p.(${A(i - k)}_${A(i - 1)}dup)`; return `p.(${A(i - 1)}_${A(i)}ins${three(ins)})`; }
  if (!ins) return del.length === 1 ? `p.(${A(i)}del)` : `p.(${A(i)}_${A(i + del.length - 1)}del)`;
  if (del.length === 1 && ins.length === 1) return `p.(${A(i)}${AA3[ins]})`;
  return del.length === 1 ? `p.(${A(i)}delins${three(ins)})` : `p.(${A(i)}_${A(i + del.length - 1)}delins${three(ins)})`;
}
module.exports = { makeGene, makeCoord, expectC, expectP, rand, ri, pick, randSeq, seedRng, CODON, tr };
