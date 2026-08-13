import os
import re
import pyhmmer

RE_PREFIX = re.compile(r"^UniRef\d+_")

def extract_clean_acc(token: str) -> str:
    """Extracts a valid accession/family name from hit header tokens."""
    token = RE_PREFIX.sub("", token.lstrip(">").strip())
    parts = token.split("|")
    acc = parts[-1] if len(parts) > 1 else parts[0]
    acc = acc.split("/")[0].split(".")[0]
    return acc if len(acc) >= 3 else token

def run_hhblits_api(
    sequence: str = "", 
    database: str = "dbCAN_HMMdb_V10", 
    iterations: int = 1,
    query_file: str = "1csm.fasta",
    output_a3m_path: str = None
) -> tuple[list[str], str]:
    """
    Searches a local HMM profile database (e.g. dbCAN) against the query sequence 
    using PyHMMER hmmscan.
    Returns a tuple of (unique_domain_hits, a3m_alignment_string).
    """
    # 1. Resolve local database file path
    target_db = database
    if not os.path.exists(target_db):
        for ext in [".txt", ".hmm", ".txt.hmm"]:
            if os.path.exists(f"{database}{ext}"):
                target_db = f"{database}{ext}"
                break

    if not os.path.exists(target_db):
        raise FileNotFoundError(f"Local HMM database file not found at: '{database}'")

    # 2. Load Query Sequence
    if os.path.exists(query_file):
        query_path = query_file
    else:
        query_path = "_temp_query.fasta"
        clean_seq = "".join([l.strip() for l in sequence.strip().splitlines() if not l.startswith(">")])
        with open(query_path, "w") as f:
            f.write(f">query\n{clean_seq}\n")

    print(f"--> Reading query sequence from '{query_path}'...")
    with pyhmmer.easel.SequenceFile(query_path, digital=True) as seq_file:
        query_seq = seq_file.read()

    # Read original query header/sequence string
    if hasattr(query_seq, "name") and query_seq.name:
        query_name = query_seq.name.decode("utf-8") if isinstance(query_seq.name, bytes) else str(query_seq.name)
    else:
        iquery_name = "query"
    query_text_seq = query_seq.sequence

    # 3. Read target HMM database and execute hmmscan
    print(f"--> Scanning query against local HMM database '{target_db}'...")
    
    hits = []
    a3m_records = [f">{query_name}\n{query_text_seq}"]

    with pyhmmer.plan7.HMMFile(target_db) as hmm_file:
        for top_hits in pyhmmer.hmmer.hmmscan([query_seq], hmm_file, cpus=4):
            for hit in top_hits:
                if getattr(hit, "included", True):
                    raw_name = hit.name
                    hmm_name = raw_name.decode("utf-8") if isinstance(raw_name, bytes) else str(raw_name)
                    cleaned = extract_clean_acc(hmm_name)
                    if cleaned:
                        hits.append(cleaned)

                    # Extract aligned target segments for A3M
                    for domain in hit.domains:
                        aligned_segment = domain.alignment.target_sequence
                        a3m_records.append(f">{hmm_name}\n{aligned_segment}")

    # Build A3M string
    a3m_string = "\n".join(a3m_records) + "\n"

    # Save A3M file if output path is requested
    if output_a3m_path and len(a3m_records) > 1:
        os.makedirs(os.path.dirname(output_a3m_path) or ".", exist_ok=True)
        with open(output_a3m_path, "w") as f:
            f.write(a3m_string)
        print(f"--> Saved HMM Scan A3M alignment to '{output_a3m_path}'")

    # Clean up temp file
    if query_path == "_temp_query.fasta" and os.path.exists("_temp_query.fasta"):
        os.remove("_temp_query.fasta")

    # Deduplicate keeping order
    seen = set()
    unique_hits = [h for h in hits if not (h in seen or seen.add(h))]
    print(f"[HMM Scan Local] Retrieved {len(unique_hits)} matching domain/family hits from {target_db}.")
    
    return unique_hits, a3m_string