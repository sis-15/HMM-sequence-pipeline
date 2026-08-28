import os
import re
import subprocess
import tempfile

RE_PREFIX = re.compile(r"^UniRef\d+_")


def extract_clean_acc(token: str) -> str:
    """Extracts a valid UniProt accession from hit header tokens in the A3M file."""
    token = RE_PREFIX.sub("", token.lstrip(">").strip())
    parts = token.split("|")
    acc = parts[-1] if len(parts) > 1 else parts[0]
    acc = acc.split("/")[0].split(".")[0]
    return acc if len(acc) >= 6 and not acc.isdigit() else ""


def parse_fasta_records(fasta_path: str) -> list[tuple[str, str]]:
    """Parses a multi-FASTA file into a list of (header, sequence) tuples."""
    records = []
    current_header = None
    current_seq = []

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    records.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            records.append((current_header, "".join(current_seq)))

    return records


def run_hhblits_api(
    sequence: str = "",
    database: str = "",
    iterations: int = 3,
    query_file: str = "1csm.fasta",
    output_a3m_path: str = "hhblits_out.a3m",
    cpus: int = 8,
) -> tuple[list[str], str]:
    """
    Runs HHblits via subprocess against HH-suite binary databases (.ffdata / .ffindex).
    Handles multi-sequence FASTA files by looping per sequence.
    Returns a tuple of (unique_accession_ids, aggregated_a3m_alignment_string).
    """
    # 1. Extract sequence records
    if os.path.exists(query_file):
        records = parse_fasta_records(query_file)
    else:
        clean_seq = "".join(
            [
                l.strip()
                for l in sequence.strip().splitlines()
                if not l.startswith(">")
            ]
        )
        records = [(">query", clean_seq)]

    a3m_out = output_a3m_path if output_a3m_path else "hhblits_out.a3m"
    os.makedirs(os.path.dirname(os.path.abspath(a3m_out)), exist_ok=True)

    # Clear target output file if it already exists
    if os.path.exists(a3m_out):
        os.remove(a3m_out)

    all_accessions = []
    aggregated_a3m_lines = []

    # 2. Loop through each sequence individually
    for idx, (header, seq) in enumerate(records, start=1):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_f:
            tmp_f.write(f"{header}\n{seq}\n")
            tmp_query_path = tmp_f.name

        tmp_a3m_out = f"{a3m_out}.seq_{idx}.tmp"

        cmd = [
            "hhblits",
            "-i",
            tmp_query_path,
            "-d",
            database,
            "-n",
            str(iterations),
            "-cpu",
            str(cpus),
            "-oa3m",
            tmp_a3m_out,
        ]

        print(f"--> [HHblits {idx}/{len(records)}] Executing CLI for {header.split()[0]}...")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except FileNotFoundError:
            if os.path.exists(tmp_query_path):
                os.remove(tmp_query_path)
            raise RuntimeError(
                "HHblits binary not found in PATH! Ensure PATH contains HH-suite."
            )
        except subprocess.CalledProcessError as e:
            print(f"HHblits STDOUT:\n{e.stdout}")
            print(f"HHblits STDERR:\n{e.stderr}")
            if os.path.exists(tmp_query_path):
                os.remove(tmp_query_path)
            raise e
        finally:
            if os.path.exists(tmp_query_path):
                os.remove(tmp_query_path)

        # 3. Parse individual result A3M
        if os.path.exists(tmp_a3m_out):
            with open(tmp_a3m_out, "r", encoding="utf-8") as f:
                content = f.read()
                aggregated_a3m_lines.append(content)
                for line in content.splitlines():
                    if line.startswith(">"):
                        acc = extract_clean_acc(line)
                        if acc:
                            all_accessions.append(acc)
            os.remove(tmp_a3m_out)

    # 4. Save combined A3M file
    full_a3m_string = "\n".join(aggregated_a3m_lines)
    with open(a3m_out, "w", encoding="utf-8") as f:
        f.write(full_a3m_string)

    # Deduplicate keeping order
    seen = set()
    unique_accs = [a for a in all_accessions if not (a in seen or seen.add(a))]
    print(f"[HHblits CLI] Finished {len(records)} sequences. Retrieved {len(unique_accs)} unique hit IDs.")

    return unique_accs, full_a3m_string