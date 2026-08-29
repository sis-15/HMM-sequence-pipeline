import os
import re
import subprocess

RE_PREFIX = re.compile(r"^UniRef\d+_")


def extract_clean_acc(token: str) -> str:
    """Extracts a valid UniProt accession from hit header tokens in the A3M file."""
    token = RE_PREFIX.sub("", token.lstrip(">").strip())
    parts = token.split("|")
    acc = parts[-1] if len(parts) > 1 else parts[0]
    acc = acc.split("/")[0].split(".")[0]
    return acc if len(acc) >= 6 and not acc.isdigit() else ""


def run_hhblits_api(
    sequence: str = "",
    database: str = "",
    iterations: int = 3,
    query_file: str = "1csm.fasta",
    output_a3m_path: str = "hhblits_out.a3m",
    cpus: int = 16,
) -> tuple[list[str], str]:
    """
    Runs actual HHblits via subprocess against HH-suite binary databases (.ffdata / .ffindex).
    Returns a tuple of (unique_accession_ids, a3m_alignment_string).
    """
    # 1. Resolve Query Sequence File
    if os.path.exists(query_file):
        query_path = query_file
    else:
        query_path = "_temp_query.fasta"
        clean_seq = "".join(
            [
                l.strip()
                for l in sequence.strip().splitlines()
                if not l.startswith(">")
            ]
        )
        with open(query_path, "w") as f:
            f.write(f">query\n{clean_seq}\n")

    # Determine destination for output A3M
    a3m_out = output_a3m_path if output_a3m_path else "hhblits_out.a3m"
    os.makedirs(os.path.dirname(os.path.abspath(a3m_out)), exist_ok=True)

    # 2. Build and execute HHblits CLI command
    cmd = [
        "hhblits",
        "-i",
        query_path,
        "-d",
        database,
        "-n",
        str(iterations),
        "-cpu",
        str(cpus),
        "-oa3m",
        a3m_out,
    ]

    print(f"--> Executing HHblits CLI: {' '.join(cmd)}")
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        raise RuntimeError(
            "HHblits binary not found in PATH! Ensure 'module load hhsuite' is in your submit_job.sh."
        )
    except subprocess.CalledProcessError as e:
        print(f"HHblits STDOUT:\n{e.stdout}")
        print(f"HHblits STDERR:\n{e.stderr}")
        raise e

    # 3. Parse generated A3M file to extract accessions and raw alignment string
    accessions = []
    a3m_string = ""

    if os.path.exists(a3m_out):
        with open(a3m_out, "r", encoding="utf-8") as f:
            a3m_string = f.read()

        for line in a3m_string.splitlines():
            if line.startswith(">"):
                acc = extract_clean_acc(line)
                if acc:
                    accessions.append(acc)

    # Clean up temporary query file if created
    if query_path == "_temp_query.fasta" and os.path.exists("_temp_query.fasta"):
        os.remove("_temp_query.fasta")

    # Deduplicate keeping order
    seen = set()
    unique_accs = [a for a in accessions if not (a in seen or seen.add(a))]
    print(f"[HHblits CLI] Retrieved {len(unique_accs)} unique hit IDs.")

    return unique_accs, a3m_string