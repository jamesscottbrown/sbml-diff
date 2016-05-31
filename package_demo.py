from sbml_diff import *

f1 = open('examples/model-comparison/mRNAselfReg1.sbml', 'r')
m1 = f1.read()

f2 = open('examples/model-comparison/mRNAselfReg2.sbml', 'r')
m2 = f2.read()

sbml_diff.diff_models([m1, m2], sbml_diff.GenerateDot(["red", "blue"], 2))
