# GeneBe relay

GeneBe's API only answers browser requests from `localhost` and `genebe.net`, so SangerScope on
GitHub Pages can't reach it directly. `genebe-relay.js` is a Cloudflare Worker that forwards
SangerScope's GeneBe requests server-side and returns the answer with CORS headers for your
site only.

- Forwards only the two GeneBe variant endpoints, i.e. genomic coordinates (e.g. `chr11-5226820-C-T`).
  No reads, traces or patient data pass through it. It stores and logs nothing.
- Rejects requests from any other website (403) and any other path (404).
- Free plan: 100,000 requests/day.

## Deploy (about 5 minutes, no installs)

1. Sign up / log in at https://dash.cloudflare.com (free plan).
2. **Workers & Pages → Create → Create Worker** (start from "Hello World"). Name it e.g.
   `sangerscope-genebe` → **Deploy**.
3. **Edit code** → replace everything with the contents of `genebe-relay.js` → **Deploy**.
4. Optional, in the worker's **Settings → Variables and Secrets**:
   - `ALLOWED_ORIGINS` (text): only if your site isn't `https://clairevaltin.github.io`
     (e.g. a custom domain). Comma-separated, no trailing slash.
   - `GENEBE_USER` and `GENEBE_KEY` (secret): your GeneBe account e-mail and API key, for a
     higher daily limit. They stay on Cloudflare and never reach the browser.
5. Copy the worker URL (`https://sangerscope-genebe.<your-account>.workers.dev`) and either
   - put it in `GB_RELAY_DEFAULT` in `index.html` (works for everyone who opens the site), or
   - paste it in SangerScope → Settings → **GeneBe relay URL** (this browser only).
