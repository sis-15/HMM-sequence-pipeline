import os
import re
import pyhmmer

RE_PREFIX = re.compile(r"^UniRef\d+_")

def extract_clean_acc(token: str) -> str:
    """Extracts a valid UniProt accession from hit header tokens."""
    token = RE_PREFIX.sub("", token.lstrip(">").strip())
    parts = token.split("|")
    acc = parts[-1] if len(parts) > 1 else parts[0]
    acc = acc.split("/")[0].split(".")[0]
    return acc if len(acc) >= 6 and not acc.isdigit() else ""

def run_jackhmmer_api(
    sequence: str = "", 
    database: str = "/mnt/d/bio_databases/uniref90.fasta", 
    iterations: int = 3,
    query_file: str = "1csm.fasta",
    output_a3m_path: str = None
) -> tuple[list[str], str]:
    """
    Runs local PyHMMER Jackhmmer directly against a target FASTA file.
    Returns a tuple of (unique_accession_ids, a3m_alignment_string).
    """
    # 1. Resolve target database path
    target_db = database
    if not os.path.exists(target_db):
        if os.path.exists("uniprot_sprot.fasta"):
            target_db = "uniprot_sprot.fasta"
            print(f"--> Target '{database}' not found. Falling back to local 'uniprot_sprot.fasta'...")
        else:
            raise FileNotFoundError(f"Target database file not found at: {database}")

    # 2. Prepare Query Sequence
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

    # 3. Stream Database and Run PyHMMER Jackhmmer
    print(f"--> Executing local Jackhmmer against '{target_db}' ({iterations} iterations)...")
    
    with pyhmmer.easel.SequenceFile(target_db, digital=True) as db_file:
        jackhmmer_result = pyhmmer.hmmer.jackhmmer(
            query_seq, 
            db_file, 
            cpus=4, 
            max_iterations=iterations
        )

    # 4. Extract Hits and MSA Alignment
    accessions = []
    last_iteration = None

    if hasattr(jackhmmer_result, "__iter__"):
        all_iterations = list(jackhmmer_result)
        if all_iterations:
            last_iteration = all_iterations[-1]
            hits_list = last_iteration.hits
    else:
        last_iteration = jackhmmer_result
        hits_list = getattr(jackhmmer_result, "hits", None)

    if hits_list:
        for hit in hits_list:
            if getattr(hit, "included", True):
                raw_name = hit.name
                target_name = raw_name.decode("utf-8") if isinstance(raw_name, bytes) else str(raw_name)
                cleaned = extract_clean_acc(target_name)
                if cleaned:
                    accessions.append(cleaned)

    # 5. Extract A3M alignment string from final MSA iteration
    a3m_string = ""
    if last_iteration:
        # Check if MSA result exists
        msa = getattr(last_iteration, "msa", None)
        if msa is None and hasattr(last_iteration, "to_msa"):
            msa = last_iteration.to_msa()

        if msa is not None:
            a3m_lines = []
            # Extract sequences from PyHMMER MSA object
            for seq in msa.sequences:
                seq_name = seq.name.decode("utf-8") if isinstance(seq.name, bytes) else str(seq.name)
                # Convert digital sequence to string if necessary
                if hasattr(seq, "text_sequence"):
                    seq_str = seq.text_sequence()
                else:
                    seq_str = str(seq)
                
                # A3M formatting: replace gap characters '.' with '-'
                seq_str = seq_str.replace(".", "-")
                a3m_lines.append(f">{seq_name}\n{seq_str}")
                
            a3m_string = "\n".join(a3m_lines) + "\n"

    # 6. Save A3M file if output path is requested
    if output_a3m_path and a3m_string:
        os.makedirs(os.path.dirname(os.path.abspath(output_a3m_path)), exist_ok=True)
        with open(output_a3m_path, "w", encoding="utf-8") as f:
            f.write(a3m_string)

    # Clean up temporary query file if generated
    if query_path == "_temp_query.fasta" and os.path.exists("_temp_query.fasta"):
        os.remove("_temp_query.fasta")

    seen = set()
    unique_accs = [a for a in accessions if not (a in seen or seen.add(a))]
    print(f"[Jackhmmer Local] Retrieved {len(unique_accs)} unique hit IDs.")
    
    return unique_accs, a3m_string