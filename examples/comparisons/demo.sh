#!/usr/bin/env sh
sbml-diff.py assignmentRuleModel.xml modifiedAssignmentRuleModel.xml | dot -Tpdf -o modified_assignment_rule.pdf
sbml-diff.py assignmentRuleModel.xml modifiedKineticLaw.xml | dot -Tpdf -o modified_kinetic_law.pdf
sbml-diff.py subStructure.xml entireStructure.xml | dot -Tpdf -o added_elements.pdf
sbml-diff.py assignmentRuleModel.xml modifiedParameter.xml modifiedKineticLaw.xml | dot -Tpdf -o triple_comparison.pdf

sbml-diff.py --params  modifiedParameterValue.xml assignmentRuleModel.xml > parameter_comparison.txt

