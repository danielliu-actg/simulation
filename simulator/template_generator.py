import re 
import os 
import sys
import requests
scriptDir = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(os.path.join(scriptDir, '../../simulator/'))
from Bio import SeqIO
from simulator.utils import load_config
config = load_config()

def generate_mason_template(coords, ref_name='hg19', template_length=330, template_num=3, shift=20, template_dir="."):
    """
    Generate reference and alternate template sequences from variant coordinates and a reference genome.
    For each variant, generate multiple template sequences shifted by a specified distance.

    Args:
        coords (list): List of variant coordinates in the format 
                       [['chr1_798959_G>A', 0.5421], ['chr1_4396579_TG>T', 0.1081], ...].
        ref_name (str): Reference genome, like hg19 or hg38.
        template_length (int): Total length of the generated template sequence.
                               (Should be >= 300 to cover paired-end reads of length 150 bp.)
        template_num (int): Number of templates to generate per variant (default is 3).
                            The templates will be evenly spaced around the variant.
        shift (int): The base pair distance between consecutive templates (default is 20).
                     (For 3 templates the offsets will be [-shift, 0, +shift].)
        template_dir (str): Directory where the FASTA files will be written.

    Returns:
        list: A list of entries, each in the format [template_prefix, vf, template_dir].
              (If multiple templates are generated for a variant, each will have a unique prefix.)
    """
    templates = []
    # Compile a regex to parse variant IDs of the form "chr1_798959_G>A"
    pattern = re.compile(r'^([^_]+)_([^_]+)_([^>]+)>(.+)$')

    # Load the reference genome from the FASTA file.
    ref_genome = SeqIO.to_dict(SeqIO.parse(config[ref_name], "fasta"))

    # Calculate the center of the template.
    center = template_length // 2

    # Create a list of offsets. For example, with template_num=3 and shift=20, offsets will be [-20, 0, 20].
    median_index = template_num // 2
    offsets = [(i - median_index) * shift for i in range(template_num)]

    for coord in coords:
        variant_id, vf = coord
        # Create a base prefix for file naming.
        base_prefix = variant_id.replace(">", "_")
        match = pattern.search(variant_id)
        if not match:
            print(f"Skipping invalid variant format: {variant_id}")
            continue

        chrom, pos_str, ref_allele, alt_allele = match.groups()
        pos = int(pos_str)

        # Get the chromosome sequence from the reference genome.
        chrom_seq = ref_genome[chrom].seq

        for i, offset in enumerate(offsets):
            # Calculate the template start position (1-indexed).
            template_start = pos - center + offset

            # Convert to 0-indexed coordinates.
            slice_start = max(0, template_start - 1)

            # Extract upstream sequence: from template start (0-indexed) to just before the variant position.
            upstream = chrom_seq[slice_start: pos - 1]

            # Determine the end of the downstream window.
            downstream_end = slice_start + template_length

            # Extract downstream sequence for the reference allele.
            downstream = chrom_seq[pos - 1 + len(ref_allele): downstream_end]

            # Assemble the reference sequence.
            ref_seq = (str(upstream) + ref_allele + str(downstream))[:template_length]

            # Adjust downstream extraction for deletion variants.
            # When the variant is a deletion, len(ref_allele) > len(alt_allele)
            # so we extend the downstream region to compensate.
            if len(ref_allele) > len(alt_allele):
                deletion_offset = len(ref_allele) - len(alt_allele)
                downstream_alt = chrom_seq[pos - 1 + len(ref_allele): downstream_end + deletion_offset]
            else:
                downstream_alt = downstream

            # Assemble the alternate sequence.
            alt_seq = (str(upstream) + alt_allele + str(downstream_alt))[:template_length]

            # Create a unique prefix for this template.
            template_prefix = f"{base_prefix}_T{i+1}"
            write_snv_fasta(template_prefix, ref_seq, alt_seq, template_dir)

        templates.append([base_prefix, vf, template_dir])

    
    merge_fasta(templates)
    return templates

def merge_fasta(templates):
    """
    Merge multiple template sequences (reference and alternative) to FASTA files.
    
    Args:
        templates (list): List of templates [base_prefix, vf, template_dir].
    """
    for template in templates:
        base_prefix, vf, template_dir = template
        cmd_ref = f'cat {template_dir}/{base_prefix}_T*_ref.fasta > {template_dir}/{base_prefix}_ref.fasta'
        cmd_alt = f'cat {template_dir}/{base_prefix}_T*_alt.fasta > {template_dir}/{base_prefix}_alt.fasta'
        os.system(cmd_ref)
        os.system(cmd_alt)

def write_snv_fasta(template_prefix, ref_seq, alt_seq, template_dir="."):
    """
    Write the reference and alternate sequences to FASTA files.
    

    Args:
        prefix (str): Prefix for the filenames.
        ref_seq (str): The reference sequence.
        alt_seq (str): The alternate sequence.
        template_dir (str): Directory (abspath) where the FASTA files will be saved.
    """
    ref_filename = os.path.join(template_dir, f"{template_prefix}_ref.fasta")
    alt_filename = os.path.join(template_dir, f"{template_prefix}_alt.fasta")
    
    with open(ref_filename, "w") as ref_file:
        ref_file.write(f">{template_prefix}_ref\n{ref_seq}\n")
    with open(alt_filename, "w") as alt_file:
        alt_file.write(f">{template_prefix}_alt\n{alt_seq}\n")

def generate_sandy_template(pairs, fusiondb='fusiongdb2', template_dir="."):
    """
    Generate fusion gene template sequences from gene pairs in FusionGDB2.
    
    For each gene pair, the code creates a FASTA file that:
      - Contains one FASTA record per matching fusion from FusionGDB2 
        (headers include the target transcript, partner transcript, and the gene names).
      - Appends one extra record that is the target gene's transcript sequence, 
        as fetched from the Ensembl REST API.
    
    Args:
        pairs (list): List of gene pairs in the format [['MARCH2', 'OR7G2'], ...].
        fusiondb (str): Fusion database configuration key (used to access config[fusiondb]).
        template_dir (str): Directory where the FASTA files will be written.
    
    Returns:
        list: A list of entries, each in the format [template_prefix, template_dir].
    """
    # Ensure output directory exists.
    if not os.path.exists(template_dir):
        os.makedirs(template_dir)

    fusions = []
    # Note: Assumes that 'config' is available in the global scope and maps fusiondb keys to file paths.
    fusiondb_file = config[fusiondb]

    with open(fusiondb_file, 'r') as infile:
        lines = infile.readlines()

    # Process each gene pair
    for target, partner in pairs:
        matched_lines = []
        # Use a dictionary to keep track of unique target transcripts for this pair.
        # key: transcript id, value: target gene (redundant here, but left for clarity)
        target_transcripts = {}

        # Iterate over the lines from the fusion database file.
        for line in lines:
            columns = line.strip().split("\t")
            if len(columns) < 16:
                continue  # Skip malformed lines

            # Assume columns[6] contains the target gene and columns[10] the partner gene.
            if (target == columns[7]) and (partner == columns[11]):
                matched_lines.append(line)
        
        if not matched_lines:
            print(f"Warning: No matches found for {target} and {partner}")
            continue
        
        # Create output FASTA file for this gene pair.
        base_prefix = f"{target}_{partner}"
        output_filename = f"{base_prefix}.fasta"
        output_path = os.path.join(template_dir, output_filename)

        with open(output_path, 'w') as outfile:
            # Write out each fusion record.
            for line in matched_lines:
                columns = line.strip().split("\t")
                if len(columns) < 16:
                    continue

                # Extract transcript IDs.
                transcript_target = columns[1]  # Target transcript ID.
                transcript_partner = columns[2]  # Partner transcript ID.
                sample_id = columns[6]
                sequence = columns[-1].strip()   # Fusion gene sequence.

                # Add target transcript to local dictionary if not already included.
                if transcript_target not in target_transcripts:
                    target_transcripts[transcript_target] = target

                # Construct FASTA header and write record.
                fasta_header = f">{transcript_target}|{transcript_partner}|{sample_id}|{target}|{partner}"
                outfile.write(f"{fasta_header}\n{sequence}\n")
            
            # Append one extra record: the target gene transcript sequence.
            for transcript, gene in target_transcripts.items():
                try:
                    seq = fetch_sequence(transcript)
                except Exception as e:
                    print(f"Error fetching sequence for transcript {transcript}: {str(e)}")
                    continue
                # Construct a header indicating that this is the target transcript.
                fasta_header = f">{gene}|{transcript}_target"
                outfile.write(f"{fasta_header}\n{seq}\n")
        
        fusions.append([base_prefix, template_dir])
    
    return fusions

def fetch_sequence(transcript_id):
    """
    Fetch the cdna sequence for a given transcript id from the GRCh37 Ensembl REST API.

    Args:
        transcript_id (str): The Ensembl transcript ID.
    
    Returns:
        str: The fetched transcript sequence.
    
    Raises:
        requests.HTTPError: If the HTTP request fails.
    """
    url = f"https://grch37.rest.ensembl.org/sequence/id/{transcript_id}?type=cdna"
    headers = {"Content-Type": "text/plain"}
    response = requests.get(url, headers=headers)
    if response.ok:
        return response.text.strip()
    else:
        response.raise_for_status()
