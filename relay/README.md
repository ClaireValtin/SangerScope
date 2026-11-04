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
   - Leave out `GENEBE_USER` / `GENEBE_KEY` unless you want a shared fallback key. Each user's own
     key from SangerScope → Settings → GeneBe account is passed through and always takes priority;
     a fallback key would be usable by anyone who scripts requests to the relay.
4b. **Turn on two-factor authentication** for the Cloudflare account (My Profile → Authentication).
    The free plan has no payment method, so it can't run up charges: over the limit it just stops.
5. Copy the worker URL (`https://sangerscope-genebe.<your-account>.workers.dev`) and either
   - put it in `GB_RELAY_DEFAULT` in `index.html` (works for everyone who opens the site), or
   - paste it in SangerScope → Settings → **GeneBe relay URL** (this browser only).
