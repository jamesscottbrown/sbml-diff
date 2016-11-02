#!/usr/bin/env sh
sbml-diff.py --hide-rules toggle.xml repressilator.xml | dot -Tpdf -o comparison_default.pdf
sbml-diff.py toggle.xml repressilator.xml | dot -Tpdf -o comparison_default_hidden_rules.pdf

sbml-diff.py --hide-rules --cartoon toggle.xml repressilator.xml | dot -Tpdf -o comparison_cartoon.pdf

sbml-diff.py --abstract -e=X,Y,Z toggle.xml repressilator.xml | dot -Tpdf -o comparison_abstract.pdf
