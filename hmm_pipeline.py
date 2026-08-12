import argparse
import os
import sys
import shlex
import yaml

from run_jackhmmer import run_jackhmmer_api
from run_hhblits import run_hhblits_api
from Fetch_full import fetch_full_sequences_from_hits
from merge_a3m import combine_a3m_files


def load_config(config_path):
    """Loads YAML configuration file safely."""
    if not os.path.exists(config_path):
        sys.exit(f"Error: Config file '{config_path}' not found.")
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="AlphaFold-style automated sequence search and MSA generation pipeline."
    )
    parser.add_argument(
        "--tool",
        nargs="+",
        choices=["jackhmmer", "hhblits", "mmseqs", "nhmmer", "all"],
        required=True,
        help="One or more tools to run sequentially, or 'all' to run all available tools.",
    )
    parser.add_argument("--query", type=str, required=True, help="Path to input FASTA file")
    parser.add_argument("--job_name", type=str, default="api_run", help="Output directory prefix")
    parser.add_argument("--iterations", type=int, default=3, help="Search iterations")

    # Cluster / Config flags
    parser.add_argument(
        "--config",
        default=None,
        help="Path to YAML cluster configuration file (e.g., config/cluster.yaml)",
    )
    parser.add_argument(
        "--database",
        default=None,
        help="Single custom database file (overrides YAML config mapping)",
    )
    
    parser.add_argument(
        "--stream",
        action="store_true",
        help="Stream target database from disk sequentially to save RAM (recommended for large DBs like UniRef90)",
    )

    return parser.parse_args()


def main():
    
    args = parse_args()

    # Create job directory
    os.makedirs(args.job_name, exist_ok=True)
    query_header, query_seq = read_sequence_from_fasta(args.query)

    # 1. Generate executable run.sh for reproducibility
    bash_script_path = os.path.join(args.job_name, "run.sh")
    reconstructed_cmd = " ".join(shlex.quote(arg) for arg in sys.argv)
    with open(bash_script_path, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("# Auto-generated run command\n\n")
        f.write(f"python3 {reconstructed_cmd}\n")
    os.chmod(bash_script_path, 0o755)

    # 2. Load YAML Configuration if provided
    config = None
    if args.config:
        print(f"[INFO] Loading cluster configuration from: {args.config}")
        config = load_config(args.config)

    # 3. Select tools to run
    tools_to_run = args.tool
    if "all" in tools_to_run:
        tools_to_run = ["jackhmmer", "hhblits", "mmseqs", "nhmmer"]

    all_hit_ids = []
    generated_a3m_files = []

    # 4. Execute Tools
    for tool in tools_to_run:
        # Determine databases to run against for this tool
        target_dbs = []
        if args.database:
            target_dbs.append(args.database)
        elif config and "databases" in config and tool in config["databases"]:
            dbs = config["databases"][tool]
            if isinstance(dbs, dict):
                target_dbs = list(dbs.values())
            elif isinstance(dbs, list):
                target_dbs = dbs
        else:
            # Fallback default if no config or database flag provided
            default_db = "/mnt/d/bio_databases/uniref90.fasta"
            if tool in ["jackhmmer", "hhblits"]:
                print(f"[WARNING] No database specified for {tool}. Using fallback: {default_db}")
            target_dbs.append(default_db)

        # Loop through each mapped database for the tool
        for db_path in target_dbs:
            db_name = os.path.splitext(os.path.basename(db_path))[0]
            print(f"\n================ Running Tool: {tool.upper()} against {db_name} ================")
            
            a3m_output_path = os.path.join(args.job_name, f"{tool}_{db_name}.a3m")
            use_streaming = args.stream or (config and config.get("execution", {}).get("stream", False))
            
            if tool == "jackhmmer":
                hits, _ = run_jackhmmer_api(
                    sequence=query_seq,
                    database=db_path,
                    iterations=args.iterations,
                    query_file=args.query,
                    output_a3m_path=a3m_output_path,
                    stream=use_streaming,
                )
                all_hit_ids.extend(hits)

            elif tool == "hhblits":
                hits, _ = run_hhblits_api(
                    sequence=query_seq,
                    database=db_path,
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

    # 5. Deduplicate Hit Accession IDs
    seen_ids = set()
    unique_hit_ids = [h for h in all_hit_ids if not (h in seen_ids or seen_ids.add(h))]
    print(f"\nRetrieved {len(unique_hit_ids)} unique hit accessions across all searches.")

    # 6. Fetch Full Native Sequences from UniProt
    output_fasta = os.path.join(args.job_name, f"{args.job_name}_full_native_sequences.fasta")
    if unique_hit_ids and not (len(tools_to_run) == 1 and tools_to_run[0] == "nhmmer"):
        print(f"\n================ Fetching Full Sequences ================")
        fetch_full_sequences_from_hits(unique_hit_ids, output_fasta)

    # 7. Merge or Build Master A3M File
    master_a3m_path = os.path.join(args.job_name, f"{args.job_name}_master.a3m")

    if generated_a3m_files:
        print(f"\n================ Merging A3M Alignments ================")
        combine_a3m_files(
            input_files=generated_a3m_files,
            output_path=master_a3m_path,
            filter_unique_seqs=True,
        )
    elif os.path.exists(output_fasta):
        print(f"\n================ Constructing Master A3M from FASTA ================")
        with open(output_fasta, "r", encoding="utf-8") as f_in, open(
            master_a3m_path, "w", encoding="utf-8"
        ) as f_out:
            f_out.write(f">{query_header}\n{query_seq}\n")
            f_out.write(f_in.read())
        print(f"--> Master A3M created at '{master_a3m_path}'")


if __name__ == "__main__":
    main()