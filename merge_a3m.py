import os
import argparse
from typing import List, Tuple

def parse_a3m(file_path: str) -> Tuple[str, list[tuple[str, str]]]:
    """
    Parses an A3M file into (query_header, query_seq) and a list of 
    aligned records: [(header, sequence), ...].
    """
    records = []
    current_header = None
    current_seq_parts = []

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq_parts)))
                current_header = line
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        if current_header is not None:
            records.append((current_header, "".join(current_seq_parts)))

    if not records:
        return "", []

    query_record = records[0]
    aligned_records = records[1:]
    return query_record, aligned_records


def combine_a3m_files(
    input_files: List[str], 
    output_path: str, 
    filter_unique_seqs: bool = True
) -> str:
    """
    Merges multiple A3M files into one master alignment.
    
    1. Keeps the first sequence as the primary query sequence (Row 1).
    2. Collects all aligned hits across input files.
    3. Deduplicates hit sequences while preserving order.
    4. Writes the unified alignment to `output_path`.
    """
    valid_files = [f for f in input_files if os.path.exists(f) and os.path.getsize(f) > 0]
    
    if not valid_files:
        raise FileNotFoundError("None of the specified input A3M files exist or are non-empty.")

    print(f"--> Merging {len(valid_files)} A3M alignment file(s)...")

    # Read primary query from the first valid file
    primary_query, primary_hits = parse_a3m(valid_files[0])
    
    all_hits = list(primary_hits)
    
    # Process remaining files
    for file_path in valid_files[1:]:
        _, hits = parse_a3m(file_path)
        all_hits.extend(hits)

    # Deduplicate entries based on sequence
    seen_sequences = set()
    deduplicated_hits = []

    if filter_unique_seqs:
        for header, seq in all_hits:
            # Normalize sequence representation for comparison
            normalized_seq = seq.upper()
            if normalized_seq not in seen_sequences:
                seen_sequences.add(normalized_seq)
                deduplicated_hits.append((header, seq))
    else:
        deduplicated_hits = all_hits

    # Format master A3M
    output_lines = [f"{primary_query[0]}\n{primary_query[1]}"]
    for header, seq in deduplicated_hits:
        output_lines.append(f"{header}\n{seq}")

    master_a3m_content = "\n".join(output_lines) + "\n"

    # Save master file
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(master_a3m_content)

    print(f"--> Master alignment created successfully at '{output_path}'")
    print(f"    - Total aligned hits retained: {len(deduplicated_hits)}")

    return master_a3m_content


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple A3M files into a master AlphaFold alignment.")
    parser.add_argument(
        "-i", "--input", 
        nargs="+", 
        required=True, 
        help="Space-separated paths to input A3M files (e.g. jack.a3m hh.a3m mmseqs.a3m)"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path where the merged master A3M file will be saved"
    )
    parser.add_argument(
        "--keep-duplicates", 
        action="store_true", 
        help="Do not deduplicate identical aligned sequences across files"
    )

    args = parser.parse_args()

    combine_a3m_files(
        input_files=args.input, 
        output_path=args.output, 
        filter_unique_seqs=not args.keep_duplicates
    )