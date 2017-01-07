## this is some basic semantic mining file example..

#hedwig bk/ examples/chr1_clusters_fix.n3 -o 393_clusters_w_hierarchy_fix.txt -A 0.05 -a fwer -l


## main process steps, first, map gene names to ontolgy terms. Then amidst this mapping, create RDF graph and at the end save it as input.n3. Mouse ontolgy will be used to further explore connections.. Label specific genes with labels according to the mouse dataset. Genes, most significant for specific groups should be labeles as D, other N

## first generate some input based on gene analysis.candidates table consists of gene names and label.


# pseudocode for this system:
# first list the files needed for BK construction.
# Then iterate through files and obtaine subjeCT object (ontology) pairs.
# In the next step, use rdflib to construct the ontology and flush it to .n3 file.
# run hedwig on those files to further explore the rules.



import sys
import rdflib
from collections import defaultdict
## input files here..

mapping_file = 'gene_association.mgi'
targets = {}

### kegg and ontology files..
tempGO = defaultdict(list)
for line in open(mapping_file,'r').readlines():
    if not line.startswith("!"):    
        parts = line.split("\t")
        if parts[7] != "" or parts[7] != None:
            tempGO[parts[8]].append(parts[4].replace(":","_"))


## from this point on, construct a rdf graph!
KT = rdflib.Namespace('http://kt.ijs.si/hedwig#')
g = rdflib.graph.Graph()
obo_uri = "http://purl.obolibrary.org/obo/"
AMP = rdflib.Namespace(obo_uri)

# # simply go through dict and for each key, add items as predicate is_a
for id, example1 in enumerate(tempGO.keys()):
    # Write to rdf graph
    u = rdflib.term.URIRef('%sTERM%s' % (obo_uri, example1))
    g.add((u, rdflib.RDF.type, KT.is_a))
    for ex in tempGO[example1]:
        annotation_uri = rdflib.term.URIRef('%s%s' % (obo_uri, rdflib.Literal(ex)))
        g.add((u, KT.annotated_with, annotation_uri))



g.serialize(destination="WHOLE_MOUSE_BK1.n3",format='n3')

## hedwig run> 

#hedwig bk/ examples/chr1_clusters_fix.n3 -o 393_clusters_w_hierarchy_fix.txt -A 0.05 -a fwer -l
