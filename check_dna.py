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

app_name = "check_dna_app"
LOG_LEVEL = logging.DEBUG
DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
# experiment with format
# Exception tuple (à la sys.exc_info) or, if no exception has occurred, None.
LOG_FORMAT = '%(asctime)s %(levelname)s %(module)s %(funcName)s %(message)s'
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=DATE_FORMAT)

def load_dna_file(file_in):
    """
    Open dna data file
    Ignore comments and header
    Spilt each line into 5 fields
    Load contenets into a list object and return list as result
    #rsid	chromosome	position	allele1	allele2
    #['rs369202065', '1', '569388', 'G', 'G']
    """

    header = []
    data = []
    with open(file_in, 'rt') as file:
        for line in file:
            if ('#' in line.lower() ) or ('rsid' in line.lower() ):
                header.append(line.strip())
            else:
                field = line.split()
                data.append(field)
    return header, data


def print_chromsome_summary(file_out, data):
    c = collections.Counter()
    count = 0
    for record in data:
        chromosome = record[1]
        c.update({int(chromosome): 1})
    file_out.write("DNA base pairs reported by chromosome:"  + '\n')
    for key in c.keys():
        file_out.write("Chromosome {:>2} reports {:>7,} base pairs or SNPs \n".format(key, c[key]))


def run_base_pair_by_chrom(data, basePair):
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


def print_base_pair_summary(file_out,data):
    c_valid = collections.Counter()
    c_invalid = collections.Counter()
    for record in data:
        allele1 = record[3]
        allele2 = record[4]

        if allele1 in 'ACGT' and allele2 in 'ACGT':
            c_valid.update( {allele1 + "-" + allele2: 1})
        else:
            c_invalid.update({allele1 + "-" + allele2: 1})

    file_out.write("Number of occurances in file of valid base pairs:\n" )
    for key in c_valid.keys():
        file_out.write("DNA Base Pair {} occurs {:>7,} times\n".format(key, c_valid[key]))
    file_out.write("Number of occurances in file of invalid base pairs:\n")
    for key in c_invalid.keys():
        file_out.write("Invalid Base Pair {} occurs {:>7,} times\n".format(key, c_invalid[key]))



def print_header_obs(numBasePairs, subject):
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

def get_dna_file_name(config):
    # get data from config file
    curr_config_section = 'PROVIDER'
    curr_config_option = 'dna_file'
    if not config.has_section(curr_config_section):
        raise configSectionError(curr_config_section)
    if not config.has_option(curr_config_section,curr_config_option):
        raise configOptionError(curr_config_section,curr_config_option)
    dna_file_name = config[curr_config_section][curr_config_option]

    if not os.access(dna_file_name, os.F_OK):
        raise dnaFileError(dna_file_name)

    return dna_file_name


'''def validate_dna_file_exists(dna_file_name):
    if not os.access(dna_file_name, os.F_OK):
        raise dnaFileError(dna_file_name)'''

def report_header(file_out, config, header, dnaTable):
    logging.info('Start of report_header function')
    AncestryFlag = 'AncestryDNA' in header[0]

    assert AncestryFlag, 'DNA file provider not supported'
    numBasePairs = len(dnaTable)
    subject = "Not Available"
    try:
        subject = config['PERSON']['name']
    except KeyError as error:
        pass

    file_out.write ("The subject is {}\n".format(subject))
    file_out.write("The file reports on {:,} DNA base pairs or SNPs\n".format(numBasePairs))
    logging.info('run_summary function done')

def run_summary(config, header, dnaTable):
    logging.info('Start of run_summary function')
    file_out_text = os.path.join(os.getcwd(), 'out', 'summary.txt')
    logging.info(file_out_text)
    file_out = open( file_out_text,'w')
    file_out.write("*" * 40 + '\n')
    report_header(file_out, config, header, dnaTable)
    file_out.write ("*" * 40 + '\n')
    print_base_pair_summary(file_out,dnaTable)
    file_out.write ("*" * 40 + '\n')
    print_chromsome_summary(file_out,dnaTable)
    file_out.write ("*" * 40 + '\n')
    file_out.close()
    logging.info('run_summary function done')


def run_coordinate_snp_detail(dnaTable, coordinateSnpDetail):
    print ("*" * 40)
    chr, pos = coordinateSnpDetail.split(",")
    for record in dnaTable:
        if record[1] == chr and record[2] == pos:
            print(record)
    print ("*" * 40)


def run_snp_id_detail(dnaTable, snpId):
    print ("*" * 40)
    for record in dnaTable:
        if record[0] == snpId:
            print(record)
    print ("*" * 40)

def get_config(config_file_name):
    validate_config_file_exists(config_file_name)
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(config_file_name)
    return config


def main(args):
    logging.info('Starting Main function')
    logging.info('get name of configuration file')
    # terminate program if it does not exist
    config_file_name = args.config
    config = get_config(config_file_name)

    logging.info('get name of dna data file from config file')
    # terminate program if it does not exist
    dna_file_name = get_dna_file_name(config)

    logging.info('load the dna file data into Python list')
    header, dnaTable = load_dna_file(dna_file_name)

    # get the argument values given at run time
    summaryReport = args.summary
    coordinateSnp = args.coordinateSnp
    snpId = args.snpId
    basePairByChrom = args.basePairByChrom
    snp_web_view = args.snpBrowse
    snp_web_auto = args.snpWebScrap

    if not summaryReport and coordinateSnp == "None" and  snpId == "None" and \
        basePairByChrom == "None" and snp_web_view == "None" and snp_web_auto == "None":
        logging.info('No arguments supplied - by default call summary report')
        run_summary(config, header, dnaTable)
        logging.info('Summary report done')

    if summaryReport:
        logging.info('About to run summary report')
        run_summary(config, dnaTable)
        logging.info('Summary report done')

    if coordinateSnp != "None":
        run_coordinate_snp_detail(dnaTable, coordinateSnp)

    if snpId != "None":
        run_snp_id_detail(dnaTable, snpId)

    if basePairByChrom != "None":
        run_base_pair_by_chrom(dnaTable, basePairByChrom )

    if snp_web_view != "None":
        snp_list = config[snp_web_view].values()  # note that snp_group is passed in here
        check_snp_view(dnaTable, snp_list)

    if snp_web_auto != "None":
        snp_list = config[snp_web_auto].values()  # note that snp_group is passed in here
        check_snp_auto(dnaTable, snp_list)


if __name__ == "__main__":

    logging.info('Create the argparse command line parser')
    parser = argparse.ArgumentParser()
    logging.info('Add the command line arugments to the parser')
    parser.add_argument('--config', '-c',
                        help='name of config file e.g. dna.ini', default = 'dna.ini')
    parser.add_argument('--summary',
                        action='store_true', default=False,
                        help="Print summary")
    parser.add_argument('--coordinateSnp',
                        type=str, default="None",
                        help="print line from dna file for a given dna coordinate "\
                             "e.g. 1,569388. "\
                             "A coordinate is specified by chromosone number, "\
                            "followed by a comma, followed by position number")
    parser.add_argument('-S', '--snpId',
                        type=str, default="None",
                        help="print the detail from file given SNP reference")
    parser.add_argument('-P', '--paste', action='store_true',
                        help="paste SNP ref from clip board and get SNP detail from file", default=False)
    parser.add_argument('--basePairByChrom', type=str,
                        help="for a given base pair e.g A-A, print count by chromosome", default="None")
    parser.add_argument('--snpBrowse', type=str,
                        help="check a list of SNP ids with a given group heading in config file", default="None")
    parser.add_argument('--snpWebScrap', type=str,
                        help="check a list of SNP ids with a given group heading in config file", default="None")
    logging.info('Parse the command line arugments')
    args = parser.parse_args()

    logging.info('Support paste of SNP Id')
    if args.paste:
        args.snpDetail = pyperclip.paste()

    logging.info('Call main function with the command line arguments')
    try:
        main(args)
        logging.info('Main function done')
    except configFileError as error_message:
        logging.error(error_message)
    except configSectionError as error_message:
        logging.error(error_message)
    except configOptionError as error_message:
        logging.error(error_message)
    except dnaFileError as error_message:
        logging.error(error_message)
    except Exception as exc:
        logging.exception("Unexpected error - cannot run program for DNA analysis")
        exit(1)
