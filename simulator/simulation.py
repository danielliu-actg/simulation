import subprocess
import glob
import os
import random
import statistics
from Bio import SeqIO
from simulator.utils import load_config
config = load_config()


def mason(templates, rl, rc, fsm, fastq_dir, fms=None, fstd=None, fmaxs=None, fmins=None):
    """
    Reads simulation using mason_simulator

    Parameters
    ----------
    templates : list 
                each list has three items 
                    prefix: prefix of ref and alt sequences
                    vf: variant frequency
                    output_dir: output directory of ref and alt sequences
    rl        : int 
                sequencing read length
    rc        : int
                sequencing read count
    fsm       : string 
                ['normal', 'uniform']
    fastq_dir : str
                output directory of fastq files
    fms       : int
                Fragment mean size (for normal mode)
    fstd      : int
                Fragment size standard deviation (for normal mode)
    fmaxs     : int
                Fragment maximum size (for uniform mode)
    fmins     : int
                Fragment minimum size (for uniform mode)                
    Returns
    -------
    {prefix}_ref_R1.fq, {prefix}_ref_R2.fq
    {prefix}_alt_R1.fq, {prefix}_alt_R2.fq

    Example
    -------
    >>> mason([['chr1_4396579_TG_T', 0.3, './'], ['chr1_4396587_G_GTAATCA', 0.1, './']], 330, 30, 150, 1000, './fastq')
    """
    processes = []
    for template in templates:
        prefix, vf, template_dir = template
        ref_count = round(rc * (1 - vf))
        alt_count = round(rc * vf)

        if fsm == "normal":
            ref_mason_cmd = (
                f"mason_simulator -ir {template_dir}/{prefix}_ref.fasta --n {ref_count} "
                f"-o {fastq_dir}/{prefix}_ref_R1.fq -or {fastq_dir}/{prefix}_ref_R2.fq --fragment-mean-size {fms} "
                f"--fragment-size-std-dev {fstd} --illumina-read-length {rl} "
                f"--illumina-prob-mismatch 0 --illumina-prob-insert 0 --illumina-prob-deletion 0 --illumina-prob-mismatch-begin 0 --illumina-prob-mismatch-end 0"
            )

            alt_mason_cmd = (
                f"mason_simulator -ir {template_dir}/{prefix}_alt.fasta --n {alt_count} "
                f"-o {fastq_dir}/{prefix}_alt_R1.fq -or {fastq_dir}/{prefix}_alt_R2.fq --fragment-mean-size {fms} "
                f"--fragment-size-std-dev {fstd} --illumina-read-length {rl} "
                f"--illumina-prob-mismatch 0 --illumina-prob-insert 0 --illumina-prob-deletion 0 --illumina-prob-mismatch-begin 0 --illumina-prob-mismatch-end 0"
            )    
            output_logs([ref_mason_cmd, alt_mason_cmd], fastq_dir)
        elif fsm == "uniform":
            ref_mason_cmd = (
                f"mason_simulator -ir {template_dir}/{prefix}_ref.fasta --n {ref_count} "
                f"-o {fastq_dir}/{prefix}_ref_R1.fq -or {fastq_dir}/{prefix}_ref_R2.fq --fragment-max-size {fmaxs} "
                f"--fragment-min-size {fmins} --illumina-read-length {rl} "
                f"--illumina-prob-mismatch 0 --illumina-prob-insert 0 --illumina-prob-deletion 0 --illumina-prob-mismatch-begin 0 --illumina-prob-mismatch-end 0"
            )

            alt_mason_cmd = (
                f"mason_simulator -ir {template_dir}/{prefix}_alt.fasta --n {alt_count} "
                f"-o {fastq_dir}/{prefix}_alt_R1.fq -or {fastq_dir}/{prefix}_alt_R2.fq --fragment-max-size {fmaxs} "
                f"--fragment-min-size {fmins} --illumina-read-length {rl} "
                f"--illumina-prob-mismatch 0 --illumina-prob-insert 0 --illumina-prob-deletion 0 --illumina-prob-mismatch-begin 0 --illumina-prob-mismatch-end 0"
            )             
            output_logs([ref_mason_cmd, alt_mason_cmd], fastq_dir)       


        processes.append(subprocess.Popen(ref_mason_cmd, shell=True))
        processes.append(subprocess.Popen(alt_mason_cmd, shell=True))

        # Wait for all processes to complete
        for process in processes:
            process.wait()

        # Get list of matching files
        # 
        r1_files = glob.glob(os.path.join(fastq_dir, f'{prefix}_ref_R1.fq')) + glob.glob(os.path.join(fastq_dir, f'{prefix}_alt_R1.fq'))
        r2_files = glob.glob(os.path.join(fastq_dir, f'{prefix}_ref_R2.fq')) + glob.glob(os.path.join(fastq_dir, f'{prefix}_alt_R2.fq'))

        # Define output files for concatenated reads
        r1_output = os.path.join(fastq_dir, f'{prefix}_R1.fq')
        r2_output = os.path.join(fastq_dir, f'{prefix}_R2.fq')

        # Perform concatenation
        concatenate_files(r1_files, r1_output)
        concatenate_files(r2_files, r2_output)
        
        # Shuffle and rename the FASTQ reads (each record is 4 lines)
        shuffle_and_rename_fastq(r1_output, r1_output, 'R1')
        shuffle_and_rename_fastq(r2_output, r2_output, 'R2')
        
        # Compress the output files
        os.system(f'gzip {r1_output}')
        os.system(f'gzip {r2_output}')

def output_logs(cmds, output_dir):
    with open(f'{output_dir}/cmds.log', 'a') as f:
        f.write('\n'.join(cmds) + '\n')

def concatenate_files(input_files, output_file):
    if not input_files:
        print(f"No matching files found for {output_file}")
        return

    try:
        with open(output_file, "w") as out_f:
            for file in input_files:
                with open(file, "r") as in_f:
                    out_f.write(in_f.read())  # Efficient binary write
        print(f"Successfully created {output_file}")
    except Exception as e:
        print(f"Error concatenating files: {e}")

def shuffle_and_rename_fastq(input_file, output_file, type):
    """
    Reads a FASTQ file, shuffles the records, renames the header with a serial number,
    and writes the result to the output file.
    """
    # Read all lines from the file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Check that the file length is a multiple of 4 (standard FASTQ format)
    if len(lines) % 4 != 0:
        print(f"Warning: FASTQ file format error in {input_file}")

    # Group every 4 lines into a record
    records = [lines[i:i+4] for i in range(0, len(lines), 4)]

    # Shuffle the records randomly
    random.shuffle(records)

    # Rename each record's header with a serial ID (e.g. @read_1, @read_2, …)
    for i, record in enumerate(records, start=1):
        if type == 'R1':
            read = "1"
        else:
            read = "2"
        # Generate and assign the new header.
        record[0] = generate_illumina_header(i, read)

    # Write the shuffled and renamed records back to the output file
    with open(output_file, "w") as f:
        for record in records:
            f.writelines(record)
    print(f"Shuffled and renamed reads in {input_file}")

# Define a function to generate the Illumina FASTQ header.
def generate_illumina_header(serial_number, read,
                             instrument="NB000001", # "NB[0-9]{6}$": ["NextSeq"],
                             run_number="1",
                             flowcell_id="H0001BGXY", # "H[A-Z,0-9]{4}BGXY$" : (["NextSeq"], "High output flow cell"),
                             lane="1",
                             tile="1101",
                             is_filtered="N",
                             control_number="0",
                             index_sequence="NNNNNN"):
    # Compute x_pos and y_pos based on the serial number.
    # x_pos increments after every 9999 reads; y_pos cycles from 0001 to 9999.
    x = ((serial_number - 1) // 9999) + 1
    y = ((serial_number - 1) % 9999) + 1
    
    # Optional: Check to ensure x and y do not exceed 9999.
    if x > 9999:
        raise ValueError("Serial number exceeds maximum allowed for a 4-digit x_pos (9999)")
    
    # Format x and y as 4-digit numbers.
    x_pos = f"{x:04d}"
    y_pos = f"{y:04d}"
    
    header = (f"@{instrument}:{run_number}:{flowcell_id}:{lane}:{tile}:"
              f"{x_pos}:{y_pos} {read}:{is_filtered}:{control_number}:{index_sequence}\n")
    return header

def compute_fasta_stats(fasta_file):
    # Extract lengths of all sequences
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    
    if not lengths:
        return None, None  # No sequences found
    
    # Calculate mean length
    mean_length = statistics.mean(lengths)
    # For a single sequence, set standard deviation to 0
    std_length = statistics.stdev(lengths) if len(lengths) > 1 else 0
    
    return int(mean_length), int(std_length)

def sandy(templates, rl, rc, fastq_dir):
    """
    Reads simulation using mason_simulator

    Parameters
    ----------
    templates    : list 
                   each list has two items [[base_prefix, template_dir], [base_prefix2, template_dir]]
    rl           : int
                   read-mean size
    rc           : int
                   read counts
    fastq_dir    : string 
                   output directory of fastq files               
    Returns
    -------
    {base_prefix}_R1_001.fastq
    {base_prefix}_R2_001.fastq

    Example
    -------
    >>> sandy([['MARCH2_OR7G2', 'template_dir'], ['FAM172A_PRR16', 'template_dir']], 150, 5000 './fastq')
    """
    for template in templates:
        base_prefix, template_dir = template
        mean_length, std_length = compute_fasta_stats(f'{template_dir}/{base_prefix}.fasta')
        cmd = f'sandy transcriptome --prefix {base_prefix} -m {rl} --fragment-mean {mean_length} --fragment-stdd {std_length} -s 1 -O fastq -n {rc} --jobs 8 --sequencing-type paired-end -o {fastq_dir} {template_dir}/{base_prefix}.fasta'
        os.system(cmd)
        output_logs([cmd], fastq_dir)       
        shuffle_and_rename_fastq(f'{fastq_dir}/{base_prefix}_R1_001.fastq', f'{fastq_dir}/{base_prefix}_R1_001.fastq', 'R1')
        shuffle_and_rename_fastq(f'{fastq_dir}/{base_prefix}_R2_001.fastq', f'{fastq_dir}/{base_prefix}_R2_001.fastq', 'R2')
        os.system(f'gzip {fastq_dir}/{base_prefix}_R1_001.fastq')
        os.system(f'gzip {fastq_dir}/{base_prefix}_R2_001.fastq')
