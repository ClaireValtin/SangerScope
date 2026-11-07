// Synthetic AB1 object (the shape parseAB1 returns): superposition of two alleles, frac2 = fraction of allele 2.
// Gaussian peaks, random heights, uniform noise; basecalls = tallest channel; Q18 where the alleles differ, Q8 at the ends.
const T_ = require('./truth');
function makeAB(a1, a2, { frac2 = 0.5, noise = 0.03, sp = 12, sigma = 2.6, lowTail = 25, lowHead = 15 } = {}) {
  const n = Math.min(a1.length, a2.length), scans = n * sp + 40;
  const tr = { A: new Float64Array(scans), C: new Float64Array(scans), G: new Float64Array(scans), T: new Float64Array(scans) };
  const ploc = new Uint16Array(n);
  for (let j = 0; j < n; j++) {
    const c = 20 + j * sp; ploc[j] = c; const h = 900 + T_.rand() * 900;
    for (const [b, f] of [[a1[j], 1 - frac2], [a2[j], frac2]]) {
      if (!'ACGT'.includes(b) || f <= 0) continue;
      for (let x = c - 8; x <= c + 8; x++) tr[b][x] += h * f * Math.exp(-((x - c) ** 2) / (2 * sigma * sigma));
    }
  }
  for (const b of 'ACGT') for (let x = 0; x < scans; x++) tr[b][x] += T_.rand() * 1800 * noise;
  const out = {}; for (const b of 'ACGT') out[b] = Int16Array.from(tr[b], v => Math.min(32000, Math.round(v)));
  let base = ''; const qual = new Uint8Array(n);
  for (let j = 0; j < n; j++) { let best = 'A'; for (const b of 'ACGT') if (out[b][ploc[j]] > out[best][ploc[j]]) best = b; base += best;
    qual[j] = j < lowHead || j >= n - lowTail ? 8 : (a1[j] === a2[j] ? 40 : 18); }
  return { base, ploc, qual, tr: out, scans };
}
module.exports = { makeAB };
