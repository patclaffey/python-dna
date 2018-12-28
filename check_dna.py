from dnatools.exceptions import *
import dnatools.exceptions
import collections
import argparse
import pyperclip
import configparser
import webbrowser
import requests
import bs4
import re
import logging
import os

#class configFileError (Exception):
    #pass
    #print ("Error: Configuration file not found")

LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s %(message)s'
LOG_LEVEL = logging.DEBUG


def load_dna_file(file_in):
    """
    Open dna data file
    Ignore comments and header
    Spilt each line into 5 fields
    Load contenets into a list object and return list as result
    #rsid	chromosome	position	allele1	allele2
    #['rs369202065', '1', '569388', 'G', 'G']
    """

    data = []
    with open(file_in, 'rt') as file:
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


def print_header(numBasePairs, subject):
    print ("The subject is {}".format(subject))
    print ("The file reports on {:,} DNA base pairs or SNPs".format(numBasePairs))

def check_snp_view(dnaTable, snp_list):
    for snp in snp_list:
        #print("The snp value from file is {}".format(snp))
        print ("*" * 40)
        found = False
        for record in dnaTable:
            if record[0] == snp:
                #print(record)
                print('{} (position {}-{}) reported in file with values {} and {}'.format(*record) )
                found = True
        if found == False:
            print('{} not reported in dna file, values not know'.format(snp) )

        webbrowser.open(config['WEB']['ncbi'] + snp)

def check_snp_auto(dnaTable, snp_list):
    for snp in snp_list:
        url = config['WEB']['ncbi'] + snp
        print ("*" * 40)
        found = False
        for record in dnaTable:
            if record[0] == snp:
                print("{:>19}:  {}".format('SNP Identifier',snp))
                print("{:>19}:  {},{}".format('Observed Alleles', record[3], record[4]))
                print("{:>19}:  {}".format('Chromosome', record[1]))
                #print('{} (position {}-{}) reported in file with values {} and {}'.format(*record) )
                found = True
        if found == False:
            print('{} not reported in dna file, values not know'.format(snp) )

        res = requests.get(url)
        #print ("status code is " + str(res.status_code))
        #print (res.text[:50])
        #print (res.request.headers)
        #print (res.headers)
        page = bs4.BeautifulSoup(res.text, features="html.parser")
        #print("page title is " + page.title.string)

        keyPos = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dt")[1].getText().strip()
        valuePos0 = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dd")[1].find_all("span")[0].getText().strip()
        valuePos1 = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dd")[1].find_all("span")[1].getText().strip()
        keyAlleles = page.find_all('dl',class_='usa-width-one-half')[0].find_all("dt")[2].getText().strip()
        valueAlleles = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dd")[2].getText().strip()
        keyVarType = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dt")[3].getText().strip()
        valueVarType = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dd")[3].find("span").getText().strip()
        keyFreq = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dt")[4].getText().strip()
        valueFreq = page.find_all('dl', class_='usa-width-one-half')[0].find_all("dd")[4].find('div').getText().strip()

        keyGene = page.find_all('dl', class_='usa-width-one-half')[1].find_all("dt")[1].getText().strip()
        valueGene = page.find_all('dl', class_='usa-width-one-half')[1].find_all("dd")[1].getText().strip()
        keyPub = page.find_all('dl', class_='usa-width-one-half')[1].find_all("dt")[2].getText().strip()
        valuePub = page.find_all('dl', class_='usa-width-one-half')[1].find_all("dd")[2].getText().strip()

        pat = re.compile(r"\s{2,}|\\n") # match 2 or more spaces or a new line character


        print ("{:>19}:  {} {}".format(keyPos, valuePos0, valuePos1))
        print ("{:>19}:  {}".format("SNP " + keyAlleles, valueAlleles))
        print ("{:>19}:  {}".format(keyVarType, valueVarType))
        print ("{:>19}:  {}".format(keyFreq, pat.sub(' ',valueFreq)))
        print ("{:>19}:  {}".format(keyGene, valueGene))
        print ("{:>19}:  {}".format(keyPub, pat.sub(' ',valuePub)))

def validate_config_file_exists(config_file_name):
    if not os.access(config_file_name, os.F_OK):
        raise configFileError(config_file_name)

def get_dna_file(config):
    # get data from config file
    curr_config_section = 'PROVIDER'
    curr_config_option = 'dna_file'
    if not config.has_section(curr_config_section):
        raise configSectionError(curr_config_section)
    if not config.has_option(curr_config_section,curr_config_option):
        raise configOptionError(curr_config_section,curr_config_option)
    dna_file_name = config[curr_config_section][curr_config_option]
    return dna_file_name

def validate_dna_file_exists(dna_file_name):
    if not os.access(dna_file_name, os.F_OK):
        raise configFileError(dna_file_name)

def main(args):

    # validate config file exists and read it
    config_file_name = args.config
    validate_config_file_exists(config_file_name)
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(config_file_name)

    # get data from config file
    dna_file_name = get_dna_file(config)
    validate_dna_file_exists(dna_file_name)

    subject = config['PERSON']['name']

    dnaTable = load_dna_file(dna_file_name)
    numBasePairs = len(dnaTable)

    # these are the argument values past in at run time
    summary = args.summary  # print summary dna report
    get_num_base_pairs = args.numBasePairs
    print_chrom_summary = args.printChromSummary
    print_snp = args.snpDetail
    base_pair_summary = args.basePairSummary
    base_pair_chrom = args.basePairChrom
    print_snp2 = args.chromPos
    snp_web_view = args.snpWebView
    snp_web_auto = args.snpWebAuto

    if get_num_base_pairs or summary:
        print ("*" * 40)
        print_header(numBasePairs, subject)

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

    if base_pair_chrom != "None":
        print ("*" * 40)
        print_base_pair_chrom(dnaTable, base_pair_chrom )

    if snp_web_view != "None":
        snp_list = config[snp_web_view].values()  # note that snp_group is passed in here
        print ("*" * 40)
        check_snp_view(dnaTable, snp_list)

    if snp_web_auto != "None":
        snp_list = config[snp_web_auto].values()  # note that snp_group is passed in here
        print ("*" * 40)
        check_snp_auto(dnaTable, snp_list)

    print ("*" * 40)




def validate_config_file_section(config):
    if not config.has_section("PROVIDERxxyy"):
        raise configSectionError('PROVIDERxxyy')




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c',
                        help='name of config file e.g. dna.ini', default = 'dna.ini')
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
    parser.add_argument('--snpWebView', type=str,
                        help="check a list of SNP ids with a given group heading in config file", default="None")
    parser.add_argument('--snpWebAuto', type=str,
                        help="check a list of SNP ids with a given group heading in config file", default="None")


    args = parser.parse_args()



    if args.paste:
        args.snpDetail = pyperclip.paste()


    try:
        main(args)
    except Exception as exc:
        logging.exception("Unexpected error - cannot run program for DNA analysis")
        exit(1)

