import os
import re
import subprocess


def extract_acc_from_jackhmmer(target_name: str, target_accession: str) -> str:
    """Extracts a clean accession from jackhmmer target_name or target_accession columns."""
    # 1. Prefer target_accession if it's a valid UniProt ID (e.g., Q9Y6X5)
    if target_accession and target_accession != "-":
        acc = target_accession.split(".")[0].split("/")[0]
        if len(acc) >= 6:
            return acc

    # 2. Parse target_name (e.g., sp|P12345|NAME or UniRef90_P12345)
    clean_name = re.sub(r"^UniRef\d+_", "", target_name)
    parts = clean_name.split("|")
    
    # Check middle or last item if pipe-separated (e.g. tr|A0A024RBG1|...)
    for p in parts:
        token = p.split(".")[0].split("/")[0]
        if len(token) >= 6 and not token.isdigit():
            return token

    return target_name.split(".")[0]


def run_jackhmmer_api(
    sequence: str = "",
    database: str = "",
    iterations: int = 3,
    query_file: str = "",
    output_tbl_path: str = "jackhmmer_out.tbl",
    cpus: int = 16,  # Matches cpu_threads in config/cluster.yaml
    extra_flags: str = "",
) -> tuple[list[str], str]:
    """
    Runs jackhmmer via subprocess against a sequence database FASTA.
    Returns a tuple of (unique_accession_ids, raw_tblout_string).
    """
    # 1. Verify Target Database Path
    if not database or not os.path.exists(database):
        raise FileNotFoundError(
            f"Target database file not found at: {database}"
        )

    # 2. Resolve Query Sequence File
    if query_file and os.path.exists(query_file):
        query_path = query_file
    else:
        query_path = "_temp_jackhmmer_query.fasta"
        clean_seq = "".join(
            [
                l.strip()
                for l in sequence.strip().splitlines()
                if not l.startswith(">")
            ]
        )
        with open(query_path, "w") as f:
            f.write(f">query\n{clean_seq}\n")

    tbl_out = output_tbl_path if output_tbl_path else "jackhmmer_out.tbl"
    os.makedirs(os.path.dirname(os.path.abspath(tbl_out)), exist_ok=True)

    # 3. Build Command
    cmd = [
        "jackhmmer",
        "-N",
        str(iterations),
        "--cpu",
        str(cpus),
        "--tblout",
        tbl_out,
    ]

    # Append any extra flags from config (e.g. --F1 0.05)
    if extra_flags:
        cmd.extend(extra_flags.split())

    cmd.extend([query_path, database])

    print(f"--> Executing Jackhmmer CLI: {' '.join(cmd)}")
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        raise RuntimeError(
            "jackhmmer binary not found in PATH! Ensure HMMER is loaded or in environment."
        )
    except subprocess.CalledProcessError as e:
        print(f"jackhmmer STDOUT:\n{e.stdout}")
        print(f"jackhmmer STDERR:\n{e.stderr}")
        raise e

    # 4. Parse Generated `--tblout` File
    accessions = []
    tbl_string = ""

    if os.path.exists(tbl_out):
        with open(tbl_out, "r", encoding="utf-8") as f:
            tbl_string = f.read()

        for line in tbl_string.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                target_name = parts[0]
                target_acc = parts[1]
                acc = extract_acc_from_jackhmmer(target_name, target_acc)
                if acc:
                    accessions.append(acc)

    # Clean temporary query file if created
    if query_path == "_temp_jackhmmer_query.fasta" and os.path.exists("_temp_jackhmmer_query.fasta"):
        os.remove("_temp_jackhmmer_query.fasta")

    # Deduplicate keeping order
    seen = set()
    unique_accs = [a for a in accessions if not (a in seen or seen.add(a))]
    print(f"[Jackhmmer CLI] Retrieved {len(unique_accs)} unique hit IDs.")

    return unique_accs, tbl_string
