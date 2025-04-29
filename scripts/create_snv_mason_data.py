#!/usr/bin/env python
import argparse
import logging

import os
import sys
scriptDir = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(scriptDir, '../../simulator/'))
from simulator.parser import vcf_parser, snv_tsv_parser
from simulator.template_generator import generate_mason_template
from simulator.simulation import mason
import glob

def main():
    """Main function to parse arguments and process VCF entries."""
    parser = argparse.ArgumentParser(
        prog=__file__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="VCF Simulator"
    )
    parser.add_argument('--input', required=True, help="Input file, it can be tsv or vcf file")
    parser.add_argument('--ref_name', required=True, help="Reference genome name, like hg19")
    parser.add_argument('--rl', required=True, type=int, help="Read length")
    parser.add_argument('--rc', required=True, type=int, help="Read count")

    parser.add_argument('--tmp_count', required=True, type=int, help="template number")
    parser.add_argument('--tmp_shift', required=True, type=int, help="template shift")

    parser.add_argument('--fsm', required=True, choices=['normal', 'uniform'], help="Fragment size mode: choose 'normal' or 'uniform'")
    # Parameters for "normal" mode
    parser.add_argument('--fms', type=int, help="Fragment mean size (for normal mode)")
    parser.add_argument('--fstd', type=int, help="Fragment size standard deviation (for normal mode)")

    # Parameters for "uniform" mode
    parser.add_argument('--fmins', type=int, help="Fragment minimum size (for uniform mode)")
    parser.add_argument('--fmaxs', type=int, help="Fragment maximum size (for uniform mode)")
    parser.add_argument('--output', required=True, type=str, help="output folder of template and fastq")
    args = parser.parse_args()

    # Validate mode-specific parameters
    if args.fsm == "normal":
        if args.fms is None or args.fstd is None:
            parser.error("For normal mode, both --fms and --fstd must be provided.")
    elif args.fsm == "uniform":
        if args.fmins is None or args.fmaxs is None:
            parser.error("For uniform mode, both --fmins and --fmaxs must be provided.")

    if args.fsm == "normal":
        template_length = args.fms
    else:
        template_length = args.fmaxs
    
    
    template_num = args.tmp_count
    input_file = args.input
    output = args.output
    ref_name = args.ref_name
    shift = args.tmp_shift
    fsm = args.fsm 
    rl = args.rl
    rc = args.rc

    template_dir = os.path.abspath(f'{output}/template')
    fastq_dir = os.path.abspath(f'{output}/fastq')
    
    if not os.path.exists(template_dir):
        os.makedirs(template_dir)  
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)  

    if input_file.endswith('.vcf'): 
        coords = vcf_parser(input_file)
    if input_file.endswith('.tsv'): 
        coords = snv_tsv_parser(input_file)

    templates = generate_mason_template(coords, ref_name, template_length, template_num=template_num, shift=shift, template_dir=template_dir)
    if args.fsm == "normal":
        mason(templates, rl, rc, fsm, fastq_dir, fms=args.fms, fstd=args.fstd)
    else:
        mason(templates, rl, rc, fsm, fastq_dir, fmaxs=args.fmaxs, fmins=args.fmins)
    output_samplesheet(fastq_dir)

def output_samplesheet(fastq_dir):
    samplesheet_path = os.path.join(fastq_dir, 'samplesheet.csv')
    r1_files = glob.glob(f'{fastq_dir}/*_R1.fq.gz')
    r2_files = glob.glob(f'{fastq_dir}/*_R2.fq.gz')
    r1_prefix = [os.path.basename(r1_file).replace('_R1.fq.gz', '') for r1_file in r1_files]

    with open(samplesheet_path, 'a') as f:
        f.write('set,type,sample,fastq_1,fastq_2\n')
        for idx, (p, r1, r2) in enumerate(zip(r1_prefix, r1_files, r2_files), start=1):
            f.write(f'batch{idx},case,{p},{r1},{r2}\n')
  

if __name__ == '__main__':
    main()
