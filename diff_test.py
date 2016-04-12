import sbml_diff
from bs4 import BeautifulSoup
from generate_dot import GenerateDot

f1 = open('examples/model-comparison/mRNAselfReg1.sbml', 'r')
html_doc1 = f1.read()
soup1 = BeautifulSoup(html_doc1, 'xml')

f2 = open('examples/model-comparison/mRNAselfReg2.sbml', 'r')
html_doc2 = f2.read()
soup2 = BeautifulSoup(html_doc2, 'xml')


sbml_diff.diff_models([soup1, soup2], GenerateDot(["red", "blue"]))
