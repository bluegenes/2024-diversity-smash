import argparse
import csv
import re
import screed
import ncbi_taxdump_utils 

WANT_TAXONOMY = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def main(args):

    # load taxonomy info 
    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    print(f"loading nodes file '{args.nodes_dmp}'")
    taxfoo.load_nodes_dmp(args.nodes_dmp)
    print(f"loading names file '{args.names_dmp}'")
    taxfoo.load_names_dmp(args.names_dmp)

    w = csv.writer(args.output)
    w.writerow(['ident', 'taxid'] + WANT_TAXONOMY)

    # get the taxid from gi #
    gi2taxid = {}
    with open(args.gi2tax, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            key = row[0]
            value = row[2]
            gi2taxid[key] = value

    # map full gi "acc" --> taxid
    lineages_count = 0
    with screed.open(args.ncbi_fasta) as seqfile:
        for record in seqfile:
            header = record.name
            acc = header.split(' ')[0]
            match = re.search(r'gi\|(\d+)\|', header)
            if match:
                gi_number = match.group(1)
                if gi_number in gi2taxid:
                    taxid = gi2taxid[gi_number]
                    lineages_count+=1
                else:
                    print(f"{gi_number}\tNo match found")
                    continue
                lin_dict = taxfoo.get_lineage_as_dict(taxid, WANT_TAXONOMY)
                if not lin_dict:
                    print(f"WARNING: taxid {taxid} not in taxdump files. Producing empty lineage.")
                row = [acc, taxid]
                for rank in WANT_TAXONOMY:
                    name = lin_dict.get(rank, '')
                    row.append(name)
                w.writerow(row)
    print(f"output {lineages_count} lineages")
                    

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Map numbers from one file to another based on matching IDs.')
    p.add_argument('gi2tax', help='NCBI taxid info (TSV format)')
    p.add_argument('ncbi_fasta', help='NCBI FASTA')
    p.add_argument('--nodes-dmp')
    p.add_argument('--names-dmp')
    p.add_argument('-o', '--output', help='output lineages file', type=argparse.FileType('wt'))

    args = p.parse_args()
    main(args)