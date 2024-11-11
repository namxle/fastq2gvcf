import re

def sort_gtf(input_file, output_file):
    # Open the input file and read lines
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Separate headers and data lines
    headers = [line for line in lines if line.startswith('#')]
    data_lines = [line for line in lines if not line.startswith('#')]

    # Define a function to extract sort keys: chromosome and start position
    def sort_key(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[3])
        
        # Use natural sorting for chromosomes (e.g., chr1, chr2, chr10)
       
        chrom_key = chrom

        return (chrom_key, start)

    # Sort data lines by chromosome and start position
    sorted_data = sorted(data_lines, key=sort_key)

    # Write headers and sorted data to the output file
    with open(output_file, 'w') as out_file:
        out_file.writelines(headers)  # Write header lines first
        out_file.writelines(sorted_data)  # Write sorted data lines

# Usage
input_file = 'assets/GCF_009035845.1_ASM903584v1_genomic.gtf'
output_file = 'assets/sorted_output.gtf'
sort_gtf(input_file, output_file)
