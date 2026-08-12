import concurrent.futures
import os
import re
import threading
import requests

# Precompiled regex for faster pattern matching
RE_PREFIX = re.compile(r"^UniRef\d+_")
RE_GAPS = re.compile(r"[-.]")

# Thread-local storage to reuse requests.Session per worker thread
thread_local = threading.local()

def get_session() -> requests.Session:
    """Returns a persistent requests.Session per thread for reusing connections."""
    if not hasattr(thread_local, "session"):
        session = requests.Session()
        session.headers.update({"User-Agent": "Full-Sequence-Fetcher/1.0"})
        thread_local.session = session
    return thread_local.session

# Extract IDs from A3M / Hits from mmseqs2
def extract_accessions_from_a3m(a3m_path: str) -> list[str]:
  """Parses accession IDs directly from an existing A3M file."""
  accessions = []
  with open(a3m_path, "r", encoding="utf-8") as f:
    for line in f:
      if line.startswith(">"):
        # Take first token (e.g., >UniRef100_A0A5B7BN98/10-150 -> A0A5B7BN98)
        first_token = line.split()[0].lstrip(">")
        clean_id = RE_PREFIX.sub("", first_token).split("/")[0]
        # Ignore numerical IDs that won't exist on UniProt
        if not clean_id.isdigit():
          accessions.append(clean_id)

  # Preserve order while removing duplicate hits
  seen = set()
  unique_accs = [
      acc for acc in accessions if not (acc in seen or seen.add(acc))
  ]
  print(
      f"[Part 1 Complete] Extracted {len(unique_accs)} unique UniProt/UniRef"
      " hit IDs."
  )
  return unique_accs


# Retrieve Full Native Sequences
def fetch_single_sequence(acc: str) -> tuple[str, str | None]:
    session = get_session()
    urls = [
        f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
        f"https://rest.uniprot.org/uniref/UniRef100_{acc}.fasta",
    ]

    for url in urls:
        try:
            resp = session.get(url, timeout=5)
            if resp.status_code == 200 and resp.text:
                lines = resp.text.splitlines()
                seq = "".join(l.strip() for l in lines if not l.startswith(">"))
                if seq:
                    return acc, seq
        except Exception:
            continue
    return acc, None


def fetch_full_sequences_from_hits(
    accessions: list[str], output_fasta: str, max_workers: int = 25
):
  """Downloads full-length sequences in parallel for all IDs."""
  print(
      f"[Part 2 Starting] Fetching {len(accessions)} full-length sequences"
      " from UniProt..."
  )

  results = {}
  completed = 0
  total = len(accessions)

  with concurrent.futures.ThreadPoolExecutor(
      max_workers=max_workers
  ) as executor:
    future_to_acc = {
        executor.submit(fetch_single_sequence, acc): acc for acc in accessions
    }

    for future in concurrent.futures.as_completed(future_to_acc):
      acc, seq = future.result()
      completed += 1
      if seq:
        results[acc] = seq

      if completed % 200 == 0 or completed == total:
        print(
            f"  Progress: {completed}/{total} retrieved ({len(results)}"
            " successful)"
        )

  # Write full-length sequences to Multi-FASTA
  with open(output_fasta, "w", encoding="utf-8") as f:
    for acc in accessions:
      if acc in results:
        f.write(f">{acc}\n{results[acc]}\n")

  print(
      f"\nDone! Saved {len(results)} full-length native sequences to"
      f" '{output_fasta}'"
  )


if __name__ == "__main__":
  # Path to your existing hits file
  a3m_file = "1csm/1csm.a3m"
  output_file = "1csm/1csm_full_native_sequences.fasta"

  # Part 1: Extract hits
  hit_ids = extract_accessions_from_a3m(a3m_file)

  # Part 2: Fetch sequences directly
  fetch_full_sequences_from_hits(hit_ids, output_file)