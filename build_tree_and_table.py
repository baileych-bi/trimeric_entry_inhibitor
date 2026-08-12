#! /usr/bin/env python

#requires the ete3 python package
#requires the muscle3 binary accessible in the shell's $PATH

import subprocess
import sys

from ete3 import Tree, TreeStyle, NodeStyle

HTML_TEMPLATE = """
<!DOCTYPE HTML>
<html>
<head>
<meta charset="utf-8" />
<style>
    table {
        text-align: center;
    }

    tr > td:first-child {
        font-weight: bold;
        text-align: left;
    }

    th {
        transform: rotate(-45deg);
        width: 3em;
    }
</style>
</head>
<body style="font-family:Arial">
%s
</body>
</html>
"""

#return index of the first char in string that is not in chars, or -1 if not found
def find_first_not_of_any(string, chars):
    for i, c in enumerate(string):
        if c not in chars:
            return i 
    return -1

#get map of variant => sequence from MUSCLE's .afa alignment file, does NOT remove gaps
def reread_alignment(afa_filename):
    sequences = {}
    current_sequence = None

    afa_file = open(afa_filename)
    while line := afa_file.readline():
        line = line.rstrip()
        if line.startswith('>'):
            current_sequence = line[1:]
        else:
            sequences.setdefault(current_sequence, []).append(line)
    afa_file.close()

    for var in sequences:
        sequences[var] = ''.join(sequences[var])

    return sequences

#take a gapped sequence and map the ungapped indices onto the gapped indices
#e. g., due to gaps we might find that residue 100 in WIV04 corresponds to
#residue 110 in the MSA, due to gaps in the aligned sequences; in this case
#gap_map[100] is 110
def map_prototype_to_msa(gapped_prototype):
    gap_map = []
    n_aas = sum(1 for aa in gapped_prototype if aa != '-') #count non-gaps
    #j is the index into the gapped sequence
    #the index into the ungapped sequence is implicitly len(gap_map)
    j = 0
    for _ in range(n_aas):
        while gapped_prototype[j] == '-': j += 1
        gap_map.append(j)
        j += 1
    return gap_map


if __name__ == "__main__":
    PROTOTYPE = 'WIV04' #root our phylogeny in this sequence and use it as the referece/top row in the mutation table
    POSITIONS = (375, 376, 405, 435, 437, 503, 508) #1-based amino acid positions for positions of interest
    #VARIANTS = ('WIV04', 'Alpha-1', 'Beta-1', 'Delta-1', 'BA1-1', 'BA2-1', 'BA5-1', 'BQ-1', 'XBB-1', 'XBB1.5-1', 'BA2.86-1', 'JN1-1', 'XDV1-1', 'KP2.3-1', 'LF7-1', 'MV1-1', 'LP8.1-1', 'NB1.8.1-1', 'KP3.1.1-1', 'XEC-1')
    SARS2_FASTA_FILENAME = 'SARS_CoV2_Variants.fasta'
    TEMP_AFA_FILENAME    = 'temp.afa'
    TEMP_PHY_FILENAME    = 'temp.phy'
    HTML_OUT_FILENAME    = 'variants_table.html'
    PHYL_OUT_FILENAME    = 'variants_tree.pdf'

    #align SARS2 variants using MUSCLE3
    #muscle -in temp.fasta -out temp.afa
    exit_code = subprocess.call(['muscle3', '-in', SARS2_FASTA_FILENAME,
        '-out', TEMP_AFA_FILENAME])

    #demand no less than success from subprocess.call()
    if exit_code != 0:
        print ('Call to muscle3 -align failed', file=sys.stderr)
        sys.exit(exit_code)

    #reread the now gapped sequences from new alignment
    sequences = reread_alignment(TEMP_AFA_FILENAME)

    #generate the map from the ungapped prototype
    gap_map = map_prototype_to_msa(sequences[PROTOTYPE])
    msa_positions = [gap_map[i - 1] for i in POSITIONS] #i - 1 because POSITIONS are 1-based

    #call MUSCLE3 again to generate a Newick tree from the MSA
    #muscle3 -maketree -in temp.afa -out temp.phy -cluster neighborjoining
    exit_code = subprocess.call(['muscle3', '-maketree', '-in', TEMP_AFA_FILENAME,
        '-out', TEMP_PHY_FILENAME, '-cluster', 'neighborjoining'])
    
    if exit_code != 0:
        print ('Call to muscle3 -maketree failed', file=sys.stderr)
        sys.exit(exit_code)

    #load the Newick tree using the ete3 package
    t = Tree(TEMP_PHY_FILENAME)

    #and root the tree in our prototype sequence (WIV04)
    t.set_outgroup(PROTOTYPE)

    #ete3's tree will render top-to-bottom with branch labels in the order
    #specified in order; because we rooted in PROTOTYPE, the PROTOYPE will
    #come first
    order = []
    for nd in t.traverse('postorder'):
        if nd.name:
            order.append(nd.name)

    prototype_seq = sequences[PROTOTYPE]

    #now we create the mutation table from top-to-bottom in the same order as
    #ete3's phylogentic tree
    tokens = ["<table>"]
    tokens.append('<tr><th></th>%s</tr>' % ''.join(["<th>%d</th>" % pos for pos in POSITIONS]))
    for var in order:
        syms = []
        seq = sequences[var]
        for pos in msa_positions:
            wt = prototype_seq[pos]
            sym = '.' if seq[pos] == wt and var != PROTOTYPE else seq[pos] 
            syms.append(sym)
        tokens.append('<tr><td>%s</td>%s</tr>' % (var, ''.join(["<td>%s</td>" % sym for sym in syms])))
    tokens.append('</table>')

    #we write the html file to HTML_OUT_FILENAME
    html = open(HTML_OUT_FILENAME, 'w')
    print (HTML_TEMPLATE % '\n'.join(tokens), file=html)
    html.close()

    #uncomment to remove the annoying circle that ete3 puts at every branch point
    #node_style = NodeStyle()
    #node_style["size"] = 0
    #for nd in t.traverse():
    #    nd.set_style(node_style)

    #... and we write the tree to PHYL_OUT_FILENAME
    t.render(PHYL_OUT_FILENAME)

    #it is left as an exercise for the reader to convert the table to .pdf and 
    #glue the tree to the left side of it in some vector image editor
    print ("\n***\nPhylogentic tree written to '%s' and mutation table to '%s'!" % (PHYL_OUT_FILENAME, HTML_OUT_FILENAME))
