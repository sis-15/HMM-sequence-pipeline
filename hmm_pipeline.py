import argparse
import os
import sys
import shlex

from run_jackhmmer import run_jackhmmer_api
from run_hhblits import run_hhblits_api
from Fetch_full import fetch_full_sequences_from_hits
from merge_a3m import combine_a3m_files


def read_sequence_from_fasta(fasta_path: str) -> tuple[str, str]:
    """Reads header and amino acid sequence string from a FASTA file."""
    header, seq_parts = "query", []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip().lstrip(">")
            else:
                seq_parts.append(line.strip())
    return header, "".join(seq_parts)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Zero-Database Cloud HMM Sequence Pipeline"
    )
    parser.add_argument(
        "--tool",
        nargs="+",
        choices=["jackhmmer", "hhblits", "mmseqs", "nhmmer", "all"],
        required=True,
        help="One or more tools to run sequentially, or 'all' to run all available tools.",
    )
    parser.add_argument(
        "--query", type=str, required=True, help="Path to input FASTA file"
    )
    parser.add_argument("--job_name", type=str, default="api_run")
    parser.add_argument("--iterations", type=int, default=3)

    parser.add_argument(
        "--database",
        type=str,
        default="/mnt/d/bio_databases/uniref90.fasta",
        help="Path to target FASTA/HMM database file",
    )
    

    args = parser.parse_args()

    os.makedirs(args.job_name, exist_ok=True)
    query_header, query_seq = read_sequence_from_fasta(args.query)
    bash_script_path = os.path.join(args.job_name, "run.sh")
    reconstructed_cmd = " ".join(shlex.quote(arg) for arg in sys.argv)
    
    with open(bash_script_path, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("# Auto-generated run command\n\n")
        f.write(f"python3 {reconstructed_cmd}\n")
        
    os.chmod(bash_script_path, 0o755)

    tools_to_run = args.tool
    if "all" in tools_to_run:
        tools_to_run = ["jackhmmer", "hhblits", "mmseqs", "nhmmer"]

    all_hit_ids = []
    generated_a3m_files = []

    # Step 1: Run Tools
    for tool in tools_to_run:
        print(f"\n================ Running Tool: {tool.upper()} ================")
        a3m_output_path = os.path.join(args.job_name, f"{tool}.a3m")

        if tool == "jackhmmer":
            hits, _ = run_jackhmmer_api(
                sequence=query_seq,
                database=args.database,
                iterations=args.iterations,
                query_file=args.query,
                output_a3m_path=a3m_output_path,
            )
            all_hit_ids.extend(hits)

        elif tool == "hhblits":
            hits, _ = run_hhblits_api(
                sequence=query_seq,
                database=args.database,
                iterations=args.iterations,
                query_file=args.query,
                output_a3m_path=a3m_output_path,
            )
            all_hit_ids.extend(hits)

        elif tool == "mmseqs":
            from mmseqs2 import extract_accessions_from_a3m, get_mmseqs_hits

            a3m_output_path = os.path.join(args.job_name, "mmseqs.a3m")
            get_mmseqs_hits(query_seq, job_name=args.job_name, output_a3m_path=a3m_output_path)
            hits = extract_accessions_from_a3m(a3m_output_path)
            all_hit_ids.extend(hits)

        if os.path.exists(a3m_output_path) and os.path.getsize(a3m_output_path) > 0:
            generated_a3m_files.append(a3m_output_path)

    # Deduplicate IDs
    seen_ids = set()
    unique_hit_ids = [
        h for h in all_hit_ids if not (h in seen_ids or seen_ids.add(h))
    ]
    print(f"\nRetrieved {len(unique_hit_ids)} unique hit accessions.")

    # Step 2: Fetch Full Native Sequences from UniProt
    output_fasta = os.path.join(
        args.job_name, f"{args.job_name}_full_native_sequences.fasta"
    )
    if unique_hit_ids and not (len(tools_to_run) == 1 and tools_to_run[0] == "nhmmer"):
        print(f"\n================ Fetching Full Sequences ================")
        fetch_full_sequences_from_hits(unique_hit_ids, output_fasta)

    # Step 3: Merge or Build Master A3M File
    master_a3m_path = os.path.join(args.job_name, f"{args.job_name}_master.a3m")
    
    if generated_a3m_files:
        print(f"\n================ Merging A3M Alignments ================")
        combine_a3m_files(
            input_files=generated_a3m_files,
            output_path=master_a3m_path,
            filter_unique_seqs=True,
        )
    elif os.path.exists(output_fasta):
        # Fallback: Construct master .a3m from query + downloaded FASTA hits
        print(f"\n================ Constructing Master A3M from FASTA ================")
        with open(output_fasta, "r", encoding="utf-8") as f_in, open(master_a3m_path, "w", encoding="utf-8") as f_out:
            f_out.write(f">{query_header}\n{query_seq}\n")
            f_out.write(f_in.read())
        print(f"--> Master A3M created at '{master_a3m_path}'")