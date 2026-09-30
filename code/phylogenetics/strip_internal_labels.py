import dendropy
t = dendropy.Tree.get(path="in.treefile", schema="newick",
                      preserve_underscores=True)
for nd in t.internal_nodes():
    nd.label = None
t.write(path="out.treefile", schema="newick",
        suppress_rooting=True, unquoted_underscores=True)