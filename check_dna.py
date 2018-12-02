import collections
import argparse
import pyperclip


def load_dna_file():
    """
    Open dna data file
    Ignore comments and header
    Spilt each line into 5 fields
    Load contenets into a list object and return list as result
    #rsid	chromosome	position	allele1	allele2
    #['rs369202065', '1', '569388', 'G', 'G']
    """
    data = []
    with open('dna.txt', 'rt') as file:
        for line in file:
            if ('#' not in line.lower() ) and ('rsid' not in line.lower() ):
                field = line.split()
                data.append(field)
    return data


def print_chromsome_summary(data):
    c = collections.Counter()
    count = 0
    for record in data:
        chromosome = record[1]
        c.update({int(chromosome): 1})
    print("DNA base pairs reported by chromosome:")
    for key in c.keys():
        print("Chromosome {:>2} reports {:>7,} base pairs or SNPs".format(key, c[key]))


def print_base_pair_chrom(data, basePair):
    c = collections.Counter()
    count = 0
    for record in data:
        chromosome = record[1]
        allele1 = record[3]
        allele2 = record[4]
        if allele1 + "-" + allele2 == basePair:
            c.update({int(chromosome): 1})
    print("DNA base pair {} occurances in file by chromosome:".format(basePair))
    for key in c.keys():
        print("Chromosome {:>2} has {:>7,} occurances of base pair {}".format(key, c[key],basePair ))


def print_base_pair_summary(data):
    c_valid = collections.Counter()
    c_invalid = collections.Counter()

    for record in data:
        allele1 = record[3]
        allele2 = record[4]

        if allele1 in 'ACGT' and allele2 in 'ACGT':
            c_valid.update( {allele1 + "-" + allele2: 1})
        else:
            c_invalid.update({allele1 + "-" + allele2: 1})

    print("Number of occurances in file of valid base pairs:")
    for key in c_valid.keys():
        print("DNA Base Pair {} occurs {:>7,} times".format(key, c_valid[key]))
    print("Number of occurances in file of invalid base pairs:")
    for key in c_invalid.keys():
        print("Invalid Base Pair {} occurs {:>7,} times".format(key, c_invalid[key]))


def main(get_num_base_pairs, print_chrom_summary, print_snp, base_pair_summary, base_pair_chrom, summary, print_snp2):

    dnaTable = load_dna_file()
    numBasePairs = len(dnaTable)
    if get_num_base_pairs or summary:
        print ("*" * 40)
        print ("The file reports on {:,} DNA base pairs or SNPs".format(numBasePairs ))

    if base_pair_summary or summary:
        print ("*" * 40)
        print_base_pair_summary(dnaTable)

    if print_chrom_summary or summary:
        print ("*" * 40)
        print_chromsome_summary(dnaTable)

    if print_snp != "None":
        print ("*" * 40)
        for record in dnaTable:
            if record[0] == print_snp:
                print(record)

    if print_snp2 != "None":
        print ("*" * 40)
        chr,pos = print_snp2.split(",")
        for record in dnaTable:
            if record[1] == chr and record[2] == pos:
                print(record)

    if base_pair_chrom != "None":
        print ("*" * 40)
        print_base_pair_chrom(dnaTable, base_pair_chrom )

    print ("*" * 40)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--summary',
                        action='store_true',
                        default=False,
                        help="Print summary")
    parser.add_argument('--numBasePairs',
                        action='store_true',
                        default=False,
                        help="Print total number of DNA base pairs reported in file")
    parser.add_argument('--basePairSummary',
                        action='store_true',
                        default=False,
                        help="Report on number of occurances of all base pair combinations")
    parser.add_argument('--printChromSummary',
                        action='store_true',
                        default=False,
                        help="Print number of DNA base pairs by chromosome")
    parser.add_argument('-S','--snpDetail',
                        type=str,
                        default = "None",
                        help="print the detail from file given SNP reference")
    parser.add_argument('--chromPos',
                        type=str,
                        default="None",
                        help="print the detail from file given chrom number and position separated by comma")
    parser.add_argument('-P', '--paste', action='store_true',
                        help="paste SNP ref from clip board and get SNP detail from file", default=False)
    parser.add_argument('--basePairChrom', type=str,
                        help="for a given base pair e.g A-A, print count by chromosome", default="None")

    args = parser.parse_args()

    if args.paste:
        args.snpDetail = pyperclip.paste()

    main( args.numBasePairs, args.printChromSummary,
          args.snpDetail, args.basePairSummary,
          args.basePairChrom, args.summary, args.chromPos)
