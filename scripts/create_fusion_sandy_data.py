#!/usr/bin/env python
import os
import sys
import logging
import argparse
scriptDir = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(scriptDir, '../../simulator/'))
from simulator.parser import fusion_tsv_parser
from simulator.template_generator import generate_sandy_template
from simulator.simulation import sandy
import glob

def main():
    """Main function to parse arguments and process VCF entries."""
    parser = argparse.ArgumentParser(
        prog=__file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="VCF Simulator"
    )
    parser.add_argument('--input', required=True, help="Input file, it can be tsv")
    parser.add_argument('--fusiondb', required=True, help="Fusion database, like fusiongdb2")
    parser.add_argument('--rl', required=True, type=int, help="Read length")
    parser.add_argument('--rc', required=True, type=int, help="Read count")
    parser.add_argument('--output', required=True, type=str, help="output folder of template and fastq")
    args = parser.parse_args()
    
    input_file = args.input
    output = args.output
    fusiondb = args.fusiondb
    rl = args.rl
    rc = args.rc

    template_dir = os.path.abspath(f'{output}/template')
    fastq_dir = os.path.abspath(f'{output}/fastq')
    
    if not os.path.exists(template_dir):
        os.makedirs(template_dir)  
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)  

    if input_file.endswith('.tsv'): 
        pairs = fusion_tsv_parser(input_file)

    templates = generate_sandy_template(pairs, fusiondb, template_dir)    
    sandy(templates, rl, rc, fastq_dir)


# def output_samplesheet(fastq_dir):
#     samplesheet_path = os.path.join(fastq_dir, 'samplesheet.csv')
#     r1_files = glob.glob(f'{fastq_dir}/*_R1.fq.gz')
#     r2_files = glob.glob(f'{fastq_dir}/*_R2.fq.gz')
#     r1_prefix = [os.path.basename(r1_file).replace('_R1.fq.gz', '') for r1_file in r1_files]

#     with open(samplesheet_path, 'a') as f:
#         f.write('set,type,sample,fastq_1,fastq_2\n')
#         for idx, (p, r1, r2) in enumerate(zip(r1_prefix, r1_files, r2_files), start=1):
#             f.write(f'batch{idx},case,{p},{r1},{r2}\n')
  

if __name__ == '__main__':
    main()