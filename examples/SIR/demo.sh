#!/usr/bin/env sh

sbml-diff.py SIRModel1.xml SIRModel2.xml SIRModel3.xml | dot -Tpdf -o default_comparison.pdf
sbml-diff.py --abstract SIRModel1.xml SIRModel2.xml SIRModel3.xml | dot -Tpdf -o abstract_comparison.pdf


sbml-diff.py --model 1 SIRModel1.xml SIRModel2.xml SIRModel3.xml | dot -Tpdf -o separate_comparison_1.pdf
sbml-diff.py --model 2 SIRModel1.xml SIRModel2.xml SIRModel3.xml | dot -Tpdf -o separate_comparison_2.pdf
sbml-diff.py --model 3 SIRModel1.xml SIRModel2.xml SIRModel3.xml | dot -Tpdf -o separate_comparison_3.pdf
