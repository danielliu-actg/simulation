import vcfpy


def vcf_parser(file_path):
    """
    Parse VCF file and extract relevant variant information,
    filtering out records with a reference genotype ('./.').

    Parameters
    ----------
    file_path : str
        Path to the VCF file.

    Returns
    -------
    list
        A list of records in the format:
        [[chr1_798959_G>A, 0.5421], [chr1_4396579_TG>T, 0.1081], ...]
    """
    coords = []
    reader = vcfpy.Reader.from_path(file_path)
    
    for record in reader:
        # Process only if the FILTER field is "PASS" or empty.
        if "PASS" in record.FILTER or not record.FILTER:
            call = record.calls[0]
            genotype = call.data.get('GT')
            # Skip the record if the genotype is reference or missing ('./.')
            if genotype == "./.":
                continue

            # Retrieve the variant ID; if the ID field is empty, you may create one from other fields.
            # Here we assume record.ID is a list and we pop one element.
            # (You may adjust this logic as needed.)
            variant_id = record.ID.pop() if record.ID else None
            if not variant_id:
                # Fallback: create an ID from CHROM, POS, REF, ALT
                variant_id = f"{record.CHROM}_{record.POS}_{record.REF}>{record.ALT[0].value}"

            # Extract the variant frequency 'VF' value.
            # Assuming call.data['VF'] is a list, we pop one element and convert it to float.
            vf_field = call.data.get('VF')
            vf_value = float(vf_field.pop()) if vf_field else 0.0

            coords.append([variant_id, vf_value])
    return coords

def snv_tsv_parser(file_path):
    """
    Parse TSV file and extract VariantID and VF (allele or variant frequency)

    Parameters
    ----------
    file_path : string
                path of the TSV file

    Returns
    -------
    [list]
        [[chr1_798959_G>A, 0.5421], [chr1_4396579_TG>T, 0.1081]...]

    Examples
    --------
    >>> snv_tsv_parser(test.tsv)
    [[chr1_798959_G>A, 0.5421], [chr1_4396579_TG>T, 0.1081]...]
    """
    coords = []
    with open(file_path, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                continue  # Skip comment lines
            
            parts = line.split("\t")
            if len(parts) != 2:  # Ensure only two tab-separated values exist
                print(f"Skipping malformed line: {line}")  # Optional debug message
                continue
            
            id, vf = parts
            try:
                vf = float(vf)  # Convert to float safely
            except ValueError:
                print(f"Skipping line with invalid float value: {line}")  # Optional debug message
                continue
            
            coords.append([id, vf])  # Append the cleaned values
    
    return coords

def fusion_tsv_parser(file_path):
    """
    Parse TSV file and extract VariantID and VF (allele or variant frequency)

    Parameters
    ----------
    file_path : string
                path of the TSV file

    Returns
    -------
    [list]
        [[targetA, papartnerA], [targetB, papartnerC]...]

    Examples
    --------
    >>> fusion_tsv_parser(test.tsv)
    [[BRCA1, BECN1], [EGFR, ERBB3]...]
    """
    pairs = []
    with open(file_path, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                continue  # Skip comment lines
            
            parts = line.split("\t")
            if len(parts) != 2:  # Ensure only two tab-separated values exist
                print(f"Skipping malformed line: {line}")  # Optional debug message
                continue
            target, partner = parts
            pairs.append([target, partner])  # Append the cleaned values
    
    return pairs