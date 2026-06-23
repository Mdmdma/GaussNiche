#!/usr/bin/env python3
"""Add a clickable outline (bookmarks) to a PDF.

Usage: add_bookmarks.py <in.pdf> <bookmarks.tsv> <out.pdf>
The TSV has one 'page<TAB>title' per line (page is 1-based).
"""
import sys
from pypdf import PdfWriter

src, tsv, dst = sys.argv[1], sys.argv[2], sys.argv[3]
# Clone the original document faithfully (append() rebuilds it in a way some
# viewers flag as "damaged"); clone_from preserves the structure.
writer = PdfWriter(clone_from=src)
n = len(writer.pages)
added = 0
with open(tsv) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line or "\t" not in line:
            continue
        page_s, title = line.split("\t", 1)
        try:
            pg = int(page_s) - 1          # TSV is 1-based; pypdf page index is 0-based
        except ValueError:
            continue
        if 0 <= pg < n:
            writer.add_outline_item(title, pg)
            added += 1
with open(dst, "wb") as out:
    writer.write(out)
print(f"[add_bookmarks] {added} outline items -> {dst}")
