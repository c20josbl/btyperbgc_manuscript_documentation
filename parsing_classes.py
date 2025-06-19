import os
from Bio import SeqIO
from collections import defaultdict

import sys
def parse_genbank_concatenated(concatenated_file, output_file):
    names = []
    classes = []
       # Parse concatenated GenBank file
    with open(concatenated_file, 'r') as gb_file:
        with open(output_file, "w") as f:
            gb_records = SeqIO.parse(gb_file, 'gb')
            for record in gb_records:
                names.append(record.name)
                f.write(record.name)
                f.write("\t")
                if record.name.startswith('AS.BTDB'):
                    category_set = set()
                    protocluster_present = False
                    for feature in record.features:
                        if feature.type == 'protocluster':
                            protocluster_present = True
                            category_set.add(feature.qualifiers.get('category', ['Unknown'])[0])
                    if protocluster_present:
                        if len(category_set) > 1:
                            category = 'mixed'
                            classes.append(category)
                            f.write(category)
                            f.write("\n")
                        else:
                            category = category_set.pop()
                            classes.append(category)
                            f.write(category)
                            f.write("\n")
                    if not protocluster_present:
                        f.write("Unknown")
                        f.write("\n")
                elif record.name.startswith('BGC'):
                    category_set = set()
                    protocluster_present = False
                    for feature in record.features:
                        if feature.type == 'protocluster':
                            protocluster_present = True
                            category_set.add(feature.qualifiers.get('category', ['Unknown'])[0])
                    if protocluster_present:
                        if len(category_set) > 1:
                            category = 'mixed'
                            classes.append(category)
                            f.write(category)
                            f.write("\n")
                        else:
                            category = category_set.pop()
                            classes.append(category)
                            f.write(category)
                            f.write("\n")
                    if not protocluster_present:
                        f.write("Unknown")
                        f.write("\n")
                elif record.name.startswith('BTDB'):
                    classes.append(record.annotations['structured_comment']['GECCO-Data']['cluster_type'])
                    f.write(record.annotations['structured_comment']['GECCO-Data']['cluster_type'])
                    f.write("\n")
    print(names)
    print(classes)


concatenated_genbank_file = "bacillus_total_mibig.gbk"
output_file = "biosynthetic_classes.tsv"
parse_genbank_concatenated(concatenated_genbank_file, output_file)
