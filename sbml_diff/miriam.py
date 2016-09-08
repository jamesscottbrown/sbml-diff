import sys
from collections import OrderedDict


def get_identifiers(obj):
    """
    Given a BeautifulSoup object, find all of the annotations of type "is"
    (rather than e.g. "isDerivedFrom", or "isHomologTo")
    """
    identifiers = set()

    if not obj.find("annotation"):
        return identifiers

    for annotation in obj.find("annotation").find_all("is"):
        for i in annotation.find_all("li"):
            if "resource" not in i.attrs:
                continue
            resource = i.attrs["resource"]
            identifiers.add(resource)
    return identifiers


def align_element(models, element_type):

    all_identifiers = OrderedDict()
    elements = []

    # Construct a list of identifiers for every id in the species
    for model in models:
        for tag in model.select(element_type):
            identifiers = get_identifiers(tag)
            id = tag["id"]

            if not identifiers:
                continue

            if id in all_identifiers.keys() and all_identifiers[id] != identifiers:
                print "Cannot match using MIRIAM identifiers: %s id %s has two or more sets of annotations"\
                      % (element_type, id)
                print "Set one: \n", get_identifiers(all_identifiers[id])
                print "Set two: \n", identifiers
                sys.exit()

            identifier_values = all_identifiers.values()
            if identifiers in identifier_values:
                rename_to = all_identifiers.keys()[identifier_values.index(identifiers)]
                if rename_to != id:
                    elements.append((model, id, rename_to))

            all_identifiers[id] = identifiers

    if len(all_identifiers.keys()) == 0:
        print("Cannot match using MIRIAM identifiers: no %s in any model has any identifier" % element_type)

    return elements


def align_models(models):
    """
    Try to match species/reactions using annotations, in addition to their ids.

    If models containing species or reactions with the same MIRIAM 'is' annotations, but different id's, then change the
    ids to be the same.

    This sounds simple, but there are many potential pitfalls.

    One problem is multiple annotations. Something may have more than one annotations because:
    * they are of different kinds (e.g. a uniprot, chebi and kegg identifier)
    * it is a mixture (example on p.94 of SBL v3 of a species that is a pool of GMP, GDP and GTP, represented by a Bag
    of 3 bqbiol:hasVersion qualifiers)
    * it has more than one identifier in same namespace (e.g. "urn:miriam:biomodels.db:MODEL1507170000" and
    "urn:miriam:biomodels.db:BIOMD0000000584")

    So two things may be each by annotated with a set of external identifiers, and if these partially overlap it is
    difficult to tell whether they are the same thing (e.g. GMP with references to a different set of databases), or
    distinct (GMP only, vs a mixture of GMP and GDP). Argubaly mixtures should be represented by bqbiol:hasPart
    qualifiers, but I don't know if this can be relied on.

    Therefore, we treat two things as identical if they have exactly the same set of annotations.

    For example, in BIOMD0000000612.xml, there are several distinct 'species' (Cyy_A, Ocy_I, Ocy_I_PTY) whose only
    annotation is that they are osteocytes ("urn:miriam:cl:CL%3A0000137"); merging based on this would be a disaster.

    Parameters
    ----------
    models

    Returns
    -------

    """
    species_to_rename = align_element(models, "species")
    for model, old_id, new_id in species_to_rename:
        # replace species ids in species definitions
        for species in model.find_all('species'):
            if species["id"] == old_id:
                species["id"] = new_id

        # replace species names in formula
        for ci in model.find_all("ci"):
            if ci.string.strip() == old_id:
                ci.string.replace_with(new_id)

        # replace speciesReference (reactant/product lists)
        for ref in model.find_all('speciesReference'):
            if ref["species"] == old_id:
                ref["species"] = new_id
        # replace modifierSpeciesReference (modifierSpecies lists)
        for ref in model.find_all('modifierSpeciesReference'):
            if ref["species"] == old_id:
                ref["species"] = new_id

    reactions_to_rename = align_element(models, "reaction")
    for model, old_id, new_id in reactions_to_rename:
        for species in model.find_all('reaction'):
            if species["id"] == old_id:
                species["id"] = new_id
