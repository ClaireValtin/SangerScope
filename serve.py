#!/usr/bin/env python3
"""Serve SangerScope at http://localhost and open it in the browser.

GeneBe's API only accepts browser requests from pages opened at http://localhost
(not 127.0.0.1, not a file:// page, not GitHub Pages). Run:

    python3 serve.py            # port 8000
    python3 serve.py 8765       # another port

Only this folder is served, and only to this computer.
"""
import functools
import http.server
import os
import sys
import threading
import webbrowser

port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
root = os.path.dirname(os.path.abspath(__file__))
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=root)
server = http.server.ThreadingHTTPServer(("localhost", port), handler)
url = f"http://localhost:{port}/index.html"
print(f"SangerScope running at {url}  (Ctrl+C to stop)")
threading.Timer(0.5, lambda: webbrowser.open(url)).start()
try:
    server.serve_forever()
except KeyboardInterrupt:
    pass
