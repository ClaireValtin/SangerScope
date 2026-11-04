/*
 * SangerScope → GeneBe relay (Cloudflare Worker).
 *
 * Why: GeneBe's API only accepts browser requests from localhost and genebe.net, so a SangerScope page on
 * GitHub Pages is refused. This worker receives the request from SangerScope, forwards it to GeneBe
 * server-side (where browser CORS rules don't apply) and returns the answer with CORS headers that allow
 * SangerScope's address — and nobody else's.
 *
 * What it forwards: only the two GeneBe variant-annotation endpoints, i.e. genomic coordinates such as
 * chr11-5226820-C-T. No reads, traces or patient data ever pass through it. It stores and logs nothing.
 *
 * Settings (Cloudflare dashboard → the worker → Settings → Variables and Secrets):
 *   ALLOWED_ORIGINS  (text, optional)    comma-separated page origins allowed to use the relay.
 *                                        Default: https://clairevaltin.github.io
 *   GENEBE_USER      (secret, optional)  GeneBe account e-mail  } fallback for users who entered no key of their
 *   GENEBE_KEY       (secret, optional)  GeneBe API key         } own. Not recommended: anyone who can call the relay
 *                                                                  (the Origin check stops browsers, not scripts) would
 *                                                                  spend this account's GeneBe quota.
 * Keys typed in SangerScope's settings are passed through to GeneBe unchanged and never stored or logged here.
 */
const UPSTREAM = 'https://api.genebe.net';
const PATHS = /^\/cloud\/api-public\/v1\/variants?$/;   // /variant (GET, one variant) and /variants (POST, batch)

export default {
  async fetch(request, env) {
    const allowed = (env.ALLOWED_ORIGINS || 'https://clairevaltin.github.io').split(',').map(s => s.trim()).filter(Boolean);
    const origin = request.headers.get('origin') || '';
    const ok = allowed.includes(origin);
    const cors = ok ? {
      'access-control-allow-origin': origin,
      'access-control-allow-methods': 'GET, POST, OPTIONS',
      'access-control-allow-headers': 'content-type, authorization',
      'access-control-max-age': '86400',
      'vary': 'Origin',
    } : {'vary': 'Origin'};
    const reply = (status, text) => new Response(text, {status, headers: {...cors, 'content-type': 'text/plain'}});

    if (!ok) return reply(403, 'Origin not allowed');
    if (request.method === 'OPTIONS') return new Response(null, {status: 204, headers: cors});
    if (request.method !== 'GET' && request.method !== 'POST') return reply(405, 'Method not allowed');

    const url = new URL(request.url);
    if (!PATHS.test(url.pathname)) return reply(404, 'Not found');

    const headers = {'accept': 'application/json'};
    if (request.method === 'POST') headers['content-type'] = 'application/json';
    // The user's own key (SangerScope → Settings → GeneBe account) always wins; the worker's key is only a fallback.
    if (request.headers.get('authorization')) headers.authorization = request.headers.get('authorization');
    else if (env.GENEBE_USER && env.GENEBE_KEY) headers.authorization = 'Basic ' + btoa(`${env.GENEBE_USER}:${env.GENEBE_KEY}`);

    let res;
    try {
      res = await fetch(UPSTREAM + url.pathname + url.search, {
        method: request.method, headers,
        body: request.method === 'POST' ? await request.text() : undefined,
      });
    } catch (e) {
      return reply(502, 'GeneBe unreachable: ' + e.message);
    }
    return new Response(res.body, {status: res.status, headers: {...cors, 'content-type': res.headers.get('content-type') || 'application/json'}});
  },
};
