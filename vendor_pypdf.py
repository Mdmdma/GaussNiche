#!/usr/bin/env python3
"""Vendor pure-python wheels into ./pylib using ONLY the standard library.

No pip required: fetch the py3-none-any wheel from PyPI (stdlib urllib, honours
http(s)_proxy env vars) and unzip it (stdlib zipfile). Default packages: pypdf
and its only runtime dep typing_extensions.
"""
import json, os, sys, tempfile, urllib.request, zipfile

dest = "pylib"
pkgs = sys.argv[1:] or ["pypdf", "typing_extensions"]
os.makedirs(dest, exist_ok=True)
for pkg in pkgs:
    meta = json.load(urllib.request.urlopen(f"https://pypi.org/pypi/{pkg}/json", timeout=60))
    url = next(u["url"] for u in meta["urls"]
               if u["packagetype"] == "bdist_wheel" and u["filename"].endswith("py3-none-any.whl"))
    whl = os.path.join(tempfile.gettempdir(), os.path.basename(url))
    urllib.request.urlretrieve(url, whl)
    with zipfile.ZipFile(whl) as z:
        z.extractall(dest)
    print("vendored", os.path.basename(url), "->", dest)
