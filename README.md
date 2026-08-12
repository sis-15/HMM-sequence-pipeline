# HMM-sequence-pipeline
An automated AlphaFold-inspired sequence search and MSA pipeline. Generates A3M alignments and native protein FASTA datasets from query sequences across Jackhmmer, HHblits, and MMseqs2.

### The contents of this markdown documented were generated with the help of Gemini

# HMM & Sequence Search Pipeline
A lightweight, local Python pipeline for bioinformatic sequence searches, native sequence retrieval, and multi-database sequence consolidation inspired by AlphaFold's sequence retrieval architecture.

# 1. Overview & Architecture
The pipeline consists of six main modules working together:
* **`hmm_pipeline.py`**: The main CLI driver script. Coordinates tool execution, writes reproducible `run.sh` scripts, manages output paths, and controls the alignment/sequence workflow.
* **`run_jackhmmer.py`**: Uses `pyhmmer` to run iterative Jackhmmer searches against unaligned FASTA databases (e.g., UniRef90, Swiss-Prot, MGnify) and extracts A3M alignments directly from MSA objects.
* **`run_hhblits.py`**: Uses `pyhmmer` / HHblits backends to scan query sequences against profile databases.
* **`mmseqs2.py`**: Queries the ColabFold MMseqs2 API to rapidly retrieve MSAs in A3M format.
* **`Fetch_full.py`**: Accepts matching accession IDs from search steps and downloads full native protein sequences directly from the UniProt API into standard FASTA format.
* **`merge_a3m.py`**: Combines multiple `.a3m` alignment files, preserves the primary query sequence as Row 1, deduplicates identical aligned sequences, and builds the Master A3M alignment for AlphaFold.
* **`merge_fastas.py`**: A fast, linear-time utility that merges and deduplicates unaligned sequences from multiple search streams into a single Master FASTA file while preserving the query sequence at Row 1.

1. **Query Input:** `1csm.fasta`
2. **Parallel Searches:**
   ├── `Jackhmmer` vs. UniRef90 / Swiss-Prot / MGnify
   ├── `HHblits` vs. dbCAN / HMM Database
   └── `MMseqs2 API Search`
3. **Dual Processing Stream**
   ├── Alignment Stream:
   │   └── Generated .a3m files -> merge_a3m.py -> <job_name>_master.a3m (For AlphaFold)
   └── Accession Stream:
       └── Extract Hit Accessions -> Fetch_full.py -> <job_name>_full_native_sequences.fasta
4. Reproducibility:
   └── Auto-generates executable <job_name>/run.sh

# 2. Requirements & Installation
### Requirements
*   Python 3.9+
*   `pyhmmer`
*   `requests`
*   `colabfold`

### Install Dependencies
```bash
pip install pyhmmer requests colabfold
```

# 3. Command Line Parameters
### Driver Pipeline (hmm_pipeline.py)
| Parameter | Type | Required | Default | Description |
| :--- | :--- | :---: | :--- | :--- |
| `--tool` | `string` | **Yes** | — | Search engine backend: `jackhmmer`, `hhblits` (HMM scan), or `mmseqs`. |
| `--query` | `string` | **Yes** | — | Path to the input FASTA sequence file. |
| `--database` | `string` | No | `/mnt/d/bio_databases/uniref90.fasta` | Path to the local target database file (`.fasta`, `.txt`, `.hmm`). |
| `--job_name` | `string` | No | `api_run` | Output directory name and prefix for results. |
| `--iterations` | `int` | No | `3` | Number of search iterations (used by Jackhmmer). |

### Master Sequence Merger (merge_fastas.py)
| Parameter | Type | Required | Default | Description |
| :--- | :--- | :---: | :--- | :--- |
| `--inputs` | `string list` | **Yes** | — | Space-separated paths to input FASTA files to merge. |
| `--output` | `string` | **Yes** | — | Destination path for the master FASTA output file. |
| `--query` | `string` | No | `None` | Path to original query FASTA (guarantees query is placed at line 1). |
| `--dedup_by` | `string` | No | `sequence` | Deduplication mode: `sequence` (identical amino acids) or `header` (identical IDs). |

### Master MSA Merger (merge_a3m.py)
| Parameter | Type | Required | Default | Description |
| :--- | :--- | :---: | :--- | :--- |
| `-i, --input` | `string list` | **Yes** | — | Space-separated paths to input .a3m files to merge. |
| `-o, --output` | `string` | **Yes** | — | Destination path for the master `.a3m` alignment file. |
| `--keep-duplicates` | `flag` | No | `False` | Disables deduplication of identical hit sequences across inputs. |

# 4. Example Commands
### Runnning Jackhmmer Against Swiss-Prot
```bash
python3 hmm_pipeline.py --tool jackhmmer --query 1csm.fasta --database uniprot_sprot.fasta --iterations 3 --job_name jackhmmer_swissprot_test
```
### Running HHBlits HMM Profile Scan Against dbCAN
```bash
python3 hmm_pipeline.py --tools hhblits --query 1csm.fasta --database dbCAN_HMMdb_V10.txt --job_name dbcan_results
```

### Running MMSeqs2
```bash
python3 hmm_pipeline.py --tools mmseqs --query 1csm.fasta --job_name mmseqs2_test
```

### Merging And Deduplicating Multi-Database Hits
```bash
python3 merge_fastas.py --inputs mmseqs_test/mmseqs_test_full_native_sequences.fasta jackhmmer_test/jackhmmer_test_full_native_sequences.fasta 1csm/1csm_full_native_sequences.fasta --output master_run/master_sequences.fasta --query 1csm.fasta --dedup_by sequence
```

### Running All Tools in a Single Multi-Search Pipeline
```bash
python3 hmm_pipeline.py --tool all --query 1csm.fasta --job_name multi_tool_run
```

### Rerunning a Saved Job via Bash Script
```bash
./test_run/run.sh
```

### Standalone Merging of A3M alignments
```bash
python3 merge_a3m.py -i test_run/jackhmmer.a3m test_run/mmseqs.a3m -o master_run/master_alignment.a3m
```

### Merging and Deduplicating Multi-Database FASTA Hits
```bash
python3 merge_fastas.py --inputs mmseqs_test/mmseqs_test_full_native_sequences.fasta jackhmmer_test/jackhmmer_test_full_native_sequences.fasta --output master_run/master_sequences.fasta --query 1csm.fasta --dedup_by sequence
```

### Quick Commands for Sequence Auditing
```bash
# Count total aligned sequences in an A3M or FASTA file
grep -c "^>" filename.fasta
grep -c "^>" filename.a3m

# Inspect first few aligned headers
grep "^>" filename.a3m | head -n 10
```
