import argparse
import os
import pyhmmer

def run_nhmmer(query_path: str, database_path: str, job_name: str, cpus: int = 4):
    """Runs nhmmer for DNA/RNA nucleotide sequence searches using pyhmmer."""
    os.makedirs(job_name, exist_ok=True)
    out_tsv = os.path.join(job_name, f"{job_name}_hits.tsv")
    out_accs = os.path.join(job_name, f"{job_name}_accessions.txt")

    # Load nucleotide query sequences in digital DNA mode
    with pyhmmer.easel.SequenceFile(query_path, digital=True, alphabet=pyhmmer.easel.Alphabet.dna()) as q_file:
        queries = list(q_file)

    # Load target database sequences
    with pyhmmer.easel.SequenceFile(database_path, digital=True, alphabet=pyhmmer.easel.Alphabet.dna()) as db_file:
        target_seqs = list(db_file)

    print(f"Running nhmmer on {len(queries)} query against {len(target_seqs)} target sequences...")
    
    # Run pyhmmer's nhmmer pipeline
    results = pyhmmer.hmmer.nhmmer(queries, target_seqs, cpus=cpus)

    accessions = set()
    with open(out_tsv, "w") as tsv_f:
        tsv_f.write("query_name\ttarget_name\tevalue\tbitscore\tstart\tend\n")
        
        for top_hits in results:
            for hit in top_hits:
                if hit.included:
                    target_name = hit.name.decode()
                    accessions.add(target_name)
                    
                    for domain in hit.domains:
                        tsv_f.write(
                            f"{top_hits.query_name.decode()}\t"
                            f"{target_name}\t"
                            f"{hit.evalue:.2e}\t"
                            f"{hit.score:.1f}\t"
                            f"{domain.alignment.target_from}\t"
                            f"{domain.alignment.target_to}\n"
                        )

    # Write accessions for downstream processing
    with open(out_accs, "w") as acc_f:
        for acc in sorted(accessions):
            acc_f.write(f"{acc}\n")

    print(f"nhmmer complete. Found {len(accessions)} unique hits.")
    print(f"Results saved to: {out_tsv}")
    print(f"Accessions saved to: {out_accs}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run nhmmer nucleotide search using PyHMMER")
    parser.add_argument("--query", required=True, help="Path to input nucleotide FASTA")
    parser.add_argument("--database", required=True, help="Path to target nucleotide database FASTA")
    parser.add_argument("--job_name", default="nhmmer_run", help="Output directory and run ID")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPU threads")
    args = parser.parse_args()

    run_nhmmer(args.query, args.database, args.job_name, args.cpus)