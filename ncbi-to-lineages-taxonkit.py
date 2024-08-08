import argparse
import csv
import re
import screed
# import ncbi_taxdump_utils 
import pytaxonkit


WANT_TAXONOMY = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
RANK_FORMATSTR = "{k};{p};{c};{o};{f};{g};{s};{t}" # this needs to match WANT_TAXONOMY, in format specified by taxonkit reformat cmd

 
def taxonkit_get_lineages_as_dict(taxidlist, ranks=WANT_TAXONOMY, formatstr=RANK_FORMATSTR, data_dir=None):
    # get lineage, taxpath for taxids using taxonkit
    n_failed = 0
    taxinfo = {}
    try:
        tk_lineage = pytaxonkit.lineage(taxidlist, fill_missing=True, pseudo_strain=True, formatstr=formatstr, threads=2, data_dir=data_dir) # , data_dir='.') use this to use this taxonomic data. BUT doesn't have 'strain' rank
    except Exception as e:
        print(f"ERROR: Failed to retrieve lineage data with taxonkit: {e}")
        return taxinfo, len(taxidlist)
    
    for taxid in taxidlist:
        taxid_row = tk_lineage[tk_lineage['TaxID'] == taxid]
        if not taxid_row.empty:
            try:
                names = taxid_row.iloc[0]['Lineage'].split(';')
                taxpath = taxid_row.iloc[0]['LineageTaxIDs'].replace(';', '|')
            except KeyError as e:
                print(f"ERROR: KeyError for taxid {taxid}: {e}")
                n_failed += 1
                continue
        else:
            names = []
            taxpath = ''
        
        n_taxids = len(taxpath.split('|'))
        if len(names) != n_taxids or len(names) != len(ranks):
            print(f"ERROR: taxonkit lineage for taxid {taxid} has mismatched lengths")
            print(f"names: {len(names)} taxids: {n_taxids} ranks: {len(ranks)}")
            n_failed += 1
            continue
        
        taxinfo[taxid] = (taxpath, names)
    
    return taxinfo, n_failed

def write_lineages(acc2taxid, taxinfo, w):
    for acc, taxid in acc2taxid.items():
        taxpath, lin_names = taxinfo[taxid]
        if not taxpath or not lin_names:
            print(f"WARNING: taxid {taxid} not in taxdump files or produced incompatible lineage. Writing empty lineage.")
            taxpath = ''
            lin_names = [''] * len(WANT_TAXONOMY)
        row = [acc, taxid, taxpath, *lin_names]
        w.writerow(row)
        

def main(args):

    w = csv.writer(args.output)
    w.writerow(['ident', 'taxpath'] + WANT_TAXONOMY)

    # get the taxid from gi #
    gi2taxid = {}
    with open(args.gi2tax, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            key = row[0]
            value = row[2]
            gi2taxid[key] = value

    lineages_count = 0
    failed_lineages = 0
    acc2taxid = {}
    n_taxids = 0
    with screed.open(args.ncbi_fasta) as seqfile:
        for n, record in enumerate(seqfile):
            acc = record.name.split(' ')[0]
            match = re.search(r'gi\|(\d+)\|', record.name)
            if match:
                # map full gi 'ident' --> taxid
                gi_number = match.group(1)
                if gi_number in gi2taxid:
                    taxid = gi2taxid[gi_number]
                    taxid = int(taxid)
                    lineages_count += 1
                    n_taxids += 1
                    acc2taxid[acc] = taxid

                    if n_taxids >= 100000:
                        # get lineages for all taxids in acc2taxid
                        taxinfo, n_fail = taxonkit_get_lineages_as_dict(acc2taxid.values(), WANT_TAXONOMY)
                        failed_lineages += n_fail
                        if taxinfo:
                            write_lineages(acc2taxid, taxinfo, w)
                        
                        print(f"Processed {n+1} sequences")
                        n_taxids = 0
                        acc2taxid = {}
                else:
                    print(f"{gi_number}\tNo match found")
                    failed_lineages += 1
    
    # get + write remaining lineages
    taxinfo, n_fail = taxonkit_get_lineages_as_dict(acc2taxid.values(), WANT_TAXONOMY)
    failed_lineages += n_fail
    if taxinfo:
        write_lineages(acc2taxid, taxinfo, w)

    print(f"output {lineages_count} lineages")
    print(f"failed {failed_lineages} lineages")
                    

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Map numbers from one file to another based on matching IDs.')
    p.add_argument('gi2tax', help='NCBI gi --> taxid mapping (TSV format)')
    p.add_argument('ncbi_fasta', help='NCBI FASTA')
    p.add_argument('--data-dir', help='directory containing NCBI taxdump data (optional; default uses version associated with pytaxonkit)')
    p.add_argument('-o', '--output', help='output lineages file', type=argparse.FileType('wt'))

    args = p.parse_args()
    main(args)