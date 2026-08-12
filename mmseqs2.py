import os
import argparse
from colabfold.colabfold import run_mmseqs2
import re

RE_PREFIX = re.compile(r"^UniRef\d+_")

def get_mmseqs_hits(
    sequence: str,
    job_name: str = "query",
    output_a3m_path: str = None,
    use_filter: bool = True,
    use_env: bool = True
) -> str:
    os.makedirs(job_name, exist_ok=True)
    
    # Use output_a3m_path if provided by pipeline, else fallback to default
    a3m_path = output_a3m_path if output_a3m_path else os.path.join(job_name, f"{job_name}.a3m")
    
    print(f"--> Querying MMseqs2 API for '{job_name}'...")
    a3m_lines = run_mmseqs2(
        x=[sequence],
        prefix=job_name,
        use_env=use_env,
        use_filter=use_filter,
        use_templates=False,
        user_agent="ColabFold-Custom-MSA/1.0"
    )
    
    with open(a3m_path, "w", encoding="utf-8") as f:
        f.write(a3m_lines[0])
        
    print(f"--> Done! Saved file to: {a3m_path}\n")
    return a3m_path

def extract_accessions_from_a3m(a3m_path: str) -> list[str]:
    """
    Parses clean UniProt / UniRef accession IDs from an A3M alignment file header.
    """
    accessions = []
    
    with open(a3m_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                # Grab the first token (e.g., >UniRef90_P01234/1-150 -> UniRef90_P01234/1-150)
                first_token = line.split()[0].lstrip(">")
                
                # Clean prefix and strip position suffixes like /1-150 or |sp|
                clean_id = RE_PREFIX.sub("", first_token).split("/")[0]
                if "|" in clean_id:
                    clean_id = clean_id.split("|")[-1]
                
                # Check for valid UniProt accession format
                clean_id = clean_id.split(".")[0]
                if len(clean_id) >= 6 and not clean_id.isdigit():
                    accessions.append(clean_id)

    # Deduplicate while preserving sequence order
    seen = set()
    unique_ids = [acc for acc in accessions if not (acc in seen or seen.add(acc))]
    return unique_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch MMseqs2 hits for a target protein sequence.")
    parser.add_argument("--job_name", type=str, default="1csm", help="Name of output folder/file prefix.")
    parser.add_argument("--no_env", action="store_true", help="Disable metagenomic environmental database search.")
    parser.add_argument("--no_filter", action="store_true", help="Disable redundancy filtering.")
    
    args = parser.parse_args()

    # Target sequence (1CSM)
    target_seq = (
        "MDFTKPETVLNLQNIRDELVRMEDSIIFKFIERSHFATCPSVYEANHPGLEIPNFKGSFLDWALSNLEIAHSRIRRFES"
        "PDETPFFPDKIQKSFLPSINYPQILAPYAPEVNYNDKIKKVYIEKIIPLISKRDGDDKNNFGSVATRDIECLQSLSRRI"
        "HFGKFVAEAKFQSDIPLYTKLIKSKDVEGIMKNITNSAVEEKILERLTKKAEVYGVDPTNESGERRITPEYLVKIYKEI"
        "VIPITKEVEVEYLLRRLEE"
    )

    get_mmseqs_hits(
        sequence=target_seq,
        job_name=args.job_name,
        use_filter=not args.no_filter,
        use_env=not args.no_env
    )