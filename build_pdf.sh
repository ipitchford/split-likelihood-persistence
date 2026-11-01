#!/usr/bin/env bash
# Rebuild the PDF from manuscript.md (Pandoc + XeLaTeX + xdvipdfmx).
set -euo pipefail
cd -- "$(dirname -- "$0")"
pandoc manuscript.md \
  --from=markdown+tex_math_single_backslash \
  --standalone \
  --include-in-header=paper_header.tex \
  --variable=geometry:margin=20mm \
  --variable=papersize:a4 --variable=fontsize:10pt \
  --variable=colorlinks:true \
  --resource-path=. --output=split_likelihood_persistence.tex
xelatex -no-pdf -interaction=nonstopmode -halt-on-error split_likelihood_persistence.tex > pdf_build.log
xelatex -no-pdf -interaction=nonstopmode -halt-on-error split_likelihood_persistence.tex >> pdf_build.log
xdvipdfmx -q -o split_likelihood_persistence.pdf split_likelihood_persistence.xdv >> pdf_build.log 2>&1
