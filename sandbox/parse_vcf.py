import pandas as pd
import sys


# Define function to parse VCF file
def parse_vcf(vcf_file):
    # Store data in a list to create DataFrame later
    vcf_data = []
    headers = []

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("##"):
                # Skip metadata lines
                continue
            elif line.startswith("#"):
                # Column headers
                headers = line.strip().split("\t")
            else:
                # Variant data
                vcf_data.append(line.strip().split("\t"))

    # Create DataFrame
    vcf_df = pd.DataFrame(vcf_data, columns=headers)
    return vcf_df


def main():
    # Get the path to the VCF file as a command line argument
    vcf_file_path = sys.argv[1]

    # Parse the VCF file and load it into a DataFrame
    vcf_df = parse_vcf(vcf_file_path)

    print(vcf_df.shape)

    # Display the first few rows to understand the structure
    print(vcf_df.head())


if __name__ == "__main__":
    main()
