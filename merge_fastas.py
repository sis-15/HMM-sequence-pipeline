import argparse
import os

def read_fasta(file_path: str) -> list[tuple[str, str]]:
    """Reads a FASTA file and returns a list of (header, sequence) tuples."""
    records = []
    current_header = None
    current_seq = []

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                # Standardize sequence to uppercase and strip gaps if unaligned sequence wanted
                current_seq.append(line.upper().replace("-", "").replace(".", ""))

        if current_header is not None:
            records.append((current_header, "".join(current_seq)))

    return records


def merge_and_deduplicate(
    input_files: list[str], 
    output_file: str, 
    query_fasta: str = None, 
    dedup_by: str = "sequence"
) -> None:
    """
    Concatenates and deduplicates multiple FASTA files into a master list.
    
    :param input_files: List of paths to input FASTA files.
    :param output_file: Destination file path.
    :param query_fasta: Optional path to query FASTA to ensure it stays at Index 0.
    :param dedup_by: 'sequence' (deduplicate exact amino acid strings) or 
                     'header' (deduplicate by sequence header ID).
    """
    seen_keys = set()
    master_records = []

    # 1. Optionally place Query sequence as the first record (Row 1)
    if query_fasta and os.path.exists(query_fasta):
        query_records = read_fasta(query_fasta)
        if query_records:
            q_head, q_seq = query_records[0]
            master_records.append((q_head, q_seq))
            key = q_seq if dedup_by == "sequence" else q_head
            seen_keys.add(key)
            print(f"--> Query sequence loaded and fixed at Index 0: {q_head}")

    # 2. Iterate through input files and aggregate hits
    total_raw_records = 0
    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"--> Warning: File '{file_path}' not found. Skipping...")
            continue

        records = read_fasta(file_path)
        total_raw_records += len(records)

        for header, seq in records:
            key = seq if dedup_by == "sequence" else header

            # Ignore empty sequences or previously seen sequences/headers
            if not seq or key in seen_keys:
                continue

            seen_keys.add(key)
            master_records.append((header, seq))

    # 3. Write out the Master FASTA File
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_file, "w", encoding="utf-8") as f:
        for header, seq in master_records:
            f.write(f"{header}\n{seq}\n")

    print("\n" + "=" * 50)
    print(f"Merge Complete!")
    print(f"Total raw sequences processed : {total_raw_records}")
    print(f"Total unique master sequences : {len(master_records)}")
    print(f"Master file saved to          : {output_file}")
    print("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate and deduplicate multiple FASTA files into a single master sequence list."
    )
    parser.add_argument(
        "--inputs", 
        nargs="+", 
        required=True, 
        help="Space-separated paths to input FASTA files (e.g. jackhmmer.fasta hhblits.fasta mgnify.fasta)"
    )
    parser.add_argument(
        "--output", 
        type=str, 
        required=True, 
        help="Destination path for master FASTA file"
    )
    parser.add_argument(
        "--query", 
        type=str, 
        default=None, 
        help="Optional path to query FASTA file to guarantee it is placed at line 1"
    )
    parser.add_argument(
        "--dedup_by", 
        choices=["sequence", "header"], 
        default="sequence", 
        help="Deduplication criterion: 'sequence' (identical amino acids) or 'header' (identical header IDs)"
    )

    args = parser.parse_args()

    merge_and_deduplicate(
        input_files=args.inputs, 
        output_file=args.output, 
        query_fasta=args.query, 
        dedup_by=args.dedup_by
    )