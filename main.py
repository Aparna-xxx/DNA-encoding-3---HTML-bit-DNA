import random

# Base mapping dictionary (all uppercase)
BASE_MAP = {
    "0": "A", "1": "G",
    "A": "C", "G": "T",
    "C": "A", "T": "G"
}

# Universal primers
FORWARD_PRIMER = "CTACACGACGCTCTTCCGATCT"
REVERSE_PRIMER = "AGATCGGAAGAGCGGTTCAGCA"

# Parameters
N = 12  # Segment length in bytes (96 bp)
L = 19  # Address length in bits (19 bp)
SEED = 2  # Fixed random seed for reproducibility

# Initialize variables for homopolymer avoidance
u1, u2, u3 = "", "", ""

def bytes_to_base_pair(byte):
    """Convert an 8-bit byte to an 8-base DNA sequence."""
    global u1, u2, u3
    binary = f"{byte:08b}"  # Convert byte to binary string
    dna_sequence = ""
    for bit in binary:
        base = BASE_MAP[bit]
        if random.random() < 0.5:
            base = BASE_MAP[base]  # Choose synonym randomly
        # Avoid homopolymer runs
        if base == u1 == u2 == u3:
            base = BASE_MAP[base]
        # Shift last three bases
        u1, u2, u3 = u2, u3, base
        dna_sequence += base
    return dna_sequence

def add_to_base_pair(number):
    """Convert a 19-bit integer to a 19-base DNA sequence."""
    global u1, u2, u3
    binary = f"{number:019b}"  # Convert number to 19-bit binary string
    dna_sequence = ""
    for bit in binary:
        base = BASE_MAP[bit]
        if random.random() < 0.5:
            base = BASE_MAP[base]  # Choose synonym randomly
        # Avoid homopolymer runs
        if base == u1 == u2 == u3:
            base = BASE_MAP[base]
        # Shift last three bases
        u1, u2, u3 = u2, u3, base
        dna_sequence += base
    return dna_sequence

def encode_to_dna(input_file, output_file):
    """Encode a binary file into DNA sequences."""
    random.seed(SEED)  # Set random seed for reproducibility
    n = 0  # Byte counter

    with open(input_file, "rb") as infile, open(output_file, "w") as outfile:
        # Write initial primer and first address
        outfile.write(FORWARD_PRIMER + add_to_base_pair(0))
        while chunk := infile.read(65536):  # Read 64 KB chunks
            for byte in chunk:
                dna_segment = bytes_to_base_pair(byte)
                outfile.write(dna_segment)
                n += 1
                # Add primers and address after every N bytes
                if n % N == 0:
                    outfile.write(REVERSE_PRIMER + "\n" + FORWARD_PRIMER + add_to_base_pair(n // N))
        # Pad last segment to N bytes
        remaining = N - (n % N)
        for _ in range(remaining):
            dna_segment = bytes_to_base_pair(random.randint(0, 255))
            outfile.write(dna_segment)
        # Write final reverse primer
        outfile.write(REVERSE_PRIMER)

# Example usage
input_file = "in.html"  # Replace with your input binary file
output_file = "out.txt"  # Replace with your desired output file
encode_to_dna(input_file, output_file)
