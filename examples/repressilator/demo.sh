#!/usr/bin/env sh
sbml-diff.py BIOMD0000000012.xml | dot -Tpdf -o default.pdf
sbml-diff.py --hide-rules BIOMD0000000012.xml  | dot -Tpdf -o hidden_rules_labelled.pdf
sbml-diff.py --label=none --hide-rules BIOMD0000000012.xml  | dot -Tpdf -o hidden_rules_unlabelled.pdf

sbml-diff.py --abstract BIOMD0000000012.xml | dot -Tpdf -o abstract.pdf
sbml-diff.py --abstract -e=X,Y,Z BIOMD0000000012.xml | dot -Tpdf -o abstract.pdf