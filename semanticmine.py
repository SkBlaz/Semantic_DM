## this is some basic semantic mining file example..

#hedwig bk/ examples/chr1_clusters_fix.n3 -o 393_clusters_w_hierarchy_fix.txt -A 0.05 -a fwer -l


## main process steps, first, map gene names to ontolgy terms. Then amidst this mapping, create RDF graph and at the end save it as input.n3. Mouse ontolgy will be used to further explore connections.. Label specific genes with labels according to the mouse dataset. Genes, most significant for specific groups should be labeles as D, other N

## first generate some input based on gene analysis.candidates table consists of gene names and label.


import sys
import rdflib
from collections import defaultdict
## input files here..


data_file = 'candidates.csv'
mapping_file = 'gene_association.mgi'
targets = {}

with open(data_file) as f:
    lines = f.readlines()
    for line in lines:
        targets[line.split(",")[0]] = line.split(",")[1].replace("\n","")

  
for k,v in targets.items():
    print(k,v)
tempGO = defaultdict(list)
for line in open(mapping_file,'r').readlines():
    if not line.startswith("!"):    
        parts = line.split("\t")[4]
        for k in targets.keys():
            if k != 'gene' and k in line:                
                tempGO[k].append(line.split("\t")[4])

#print(tempGO)

## from this point on, construct a rdf graph!


g = rdflib.graph.Graph()
KT = rdflib.Namespace('http://kt.ijs.si/hedwig#')
amp_uri = 'http://kt.ijs.si/ontology/hedwig#'
obo_uri = "http://purl.obolibrary.org/obo/"
AMP = rdflib.Namespace(amp_uri)


for id, example1 in enumerate(targets.keys()):
    # Write to rdf graph
    u = rdflib.term.URIRef('%sgene%s' % (amp_uri, example1))
    g.add((u, rdflib.RDF.type, KT.Example))
    g.add((u, KT.class_label, rdflib.Literal(targets[example1])))
    for ex in tempGO[example1]:
        annotation_uri = rdflib.term.URIRef('%s%s' % (obo_uri, rdflib.Literal(ex)))

        blank = rdflib.BNode()
        g.add((u, KT.annotated_with, blank))
        g.add((blank, KT.annotation, annotation_uri))


g.serialize(destination="backgroundONT4.n3",format='n3')

## hedwig run> 

#hedwig bk/ examples/chr1_clusters_fix.n3 -o 393_clusters_w_hierarchy_fix.txt -A 0.05 -a fwer -l
