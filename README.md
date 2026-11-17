# SangerScope

Sanger sequencing analysis in the browser, from raw traces to annotated variants.

**Open it:** https://clairevaltin.github.io/SangerScope/

## What it does

- **Load a gene** – fetches the MANE Select reference sequence (GRCh38, RefSeq) and runs integrity checks on it
- **Add reads** – drop AB1 or FASTA files, or paste sequences; several at once
- **Align and call variants** – SNVs, deletions, insertions, duplications and delins, including heterozygous calls from AB1 traces
- **Coverage map** of exons with reads, the full CDS or the whole gene
- **Annotate** each variant with VariantValidator, ClinVar, GeneBe, gnomAD v4 and MobiDetails
- **Inspect traces** by clicking any base in the alignment
- **Export** the results to Excel

**Privacy:** reads, traces and the patient name never leave your browser. Only variant descriptions and genomic coordinates are sent to the annotation services.

## How to use

Open the link above. Nothing to install.

To run it locally, use `python3 serve.py`. This is needed for GeneBe annotation without the relay, because GeneBe only accepts requests from `localhost`.

## Repository contents

- `index.html` – the whole app
- `relay/genebe-relay.js` – Cloudflare Worker that lets the hosted page reach GeneBe
- `tests/` – variant-calling regression tests (`node tests/run.js`)
- `serve.py` – local server

## Disclaimer

**This tool is NOT validated for clinical or diagnostic use.** It is provided for research and educational purposes only, without any warranty.

## License & citation

© 2026 Rehioui Soufiyan. All rights reserved.

You're welcome to explore this code. If you'd like to use, adapt, or share it, please get in touch first: rehiouisoufiyan@gmail.com

If this work helps you, please cite it or link back to this repository.
