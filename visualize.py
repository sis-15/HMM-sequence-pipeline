import argparse
import os
import re
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

# Single-character tokens for Logomaker compatibility
GROUPS_CONFIG = {
    "H": ("Hydrophobic (LVIMC)", ["L", "V", "I", "M", "C"], "black"),
    "Ω": ("Aromatic (FWY)", ["F", "W", "Y"], "purple"),
    "-": ("Acidic/Polar (EDNQ)", ["E", "D", "N", "Q"], "red"),
    "+": ("Basic (KRH)", ["K", "R", "H"], "blue"),
    "o": ("Small/Polar (STAG)", ["S", "T", "A", "G"], "green"),
    "P": ("Proline (P)", ["P"], "orange"),
}

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def stream_and_process_a3m(a3m_path):
  """Streams an A3M file line-by-line to handle giant files without memory issues.

  Calculates amino acid counts directly without storing all sequences in RAM.
  """
  aa_to_idx = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

  counts_matrix = None
  seq_count = 0
  align_length = 0

  print(f"--> Streaming and parsing A3M file: {a3m_path}")

  # Stream fasta records one by one
  for record in SeqIO.parse(a3m_path, "fasta"):
    # Remove insertion states (lowercase) and deletion dots
    clean_seq = re.sub(r"[a-z.]", "", str(record.seq))

    # Initialize count matrix based on the first sequence length
    if counts_matrix is None:
      align_length = len(clean_seq)
      counts_matrix = np.zeros((align_length, len(AMINO_ACIDS)), dtype=int)
      print(
          "--> Alignment length determined by query:"
          f" {align_length} positions"
      )

    # Skip sequences that don't match expected match-state length
    if len(clean_seq) != align_length:
      continue

    # Increment counts in matrix
    for pos, char in enumerate(clean_seq):
      if char in aa_to_idx:
        counts_matrix[pos, aa_to_idx[char]] += 1

    seq_count += 1
    if seq_count % 10000 == 0:
      print(f"    Processed {seq_count} sequences...")

  print(f"--> Done! Total aligned sequences processed: {seq_count}")

  counts_df = pd.DataFrame(counts_matrix, columns=AMINO_ACIDS)
  total_per_col = counts_df.sum(axis=1).replace(0, 1)
  prob_df = counts_df.div(total_per_col, axis=0)

  # Build Grouped Matrix
  grouped_df = pd.DataFrame(index=prob_df.index)
  for char_code, (label, aa_list, color) in GROUPS_CONFIG.items():
    grouped_df[char_code] = prob_df[aa_list].sum(axis=1)

  return counts_df, prob_df, grouped_df, align_length


def generate_logomaker_plot(
    raw_prob_df, grouped_prob_df, align_length, args
):
  """Generates wrapped sequence logo plots using Logomaker."""
  import logomaker

  print("--> Generating Logomaker Sequence Logos...")

  num_chunks = int(np.ceil(align_length / args.chunk_size))

  fig_width = max(12, args.chunk_size * args.width_per_pos)
  fig_height = num_chunks * 2 * args.row_height

  fig, axes = plt.subplots(
      nrows=num_chunks * 2,
      ncols=1,
      figsize=(fig_width, fig_height),
      sharey=True,
  )

  if num_chunks * 2 == 2:
    axes = np.array([axes[0], axes[1]])

  group_colors = {code: config[2] for code, config in GROUPS_CONFIG.items()}

  for i in range(num_chunks):
    start_pos = i * args.chunk_size
    end_pos = min((i + 1) * args.chunk_size, align_length)

    ax_top = axes[i * 2]
    ax_bottom = axes[i * 2 + 1]

    chunk_raw = raw_prob_df.iloc[start_pos:end_pos]
    chunk_grouped = grouped_prob_df.iloc[start_pos:end_pos]

    logomaker.Logo(chunk_raw, ax=ax_top, color_scheme="chemistry", vpad=0.05)
    ax_top.set_title(
        f"Positions {start_pos} - {end_pos} (Standard AAs)",
        fontsize=11,
        fontweight="bold",
    )
    ax_top.set_ylabel("Frequency")

    logomaker.Logo(
        chunk_grouped, ax=ax_bottom, color_scheme=group_colors, vpad=0.05
    )
    ax_bottom.set_title(
        f"Positions {start_pos} - {end_pos} (Grouped Properties)",
        fontsize=11,
        fontweight="bold",
    )
    ax_bottom.set_ylabel("Frequency")

    for ax in [ax_top, ax_bottom]:
      ax.spines["top"].set_visible(False)
      ax.spines["right"].set_visible(False)
      ax.set_xlim([start_pos - 0.5, end_pos - 0.5])

  legend_patches = [
      mpatches.Patch(color=config[2], label=f"'{code}': {config[0]}")
      for code, config in GROUPS_CONFIG.items()
  ]
  axes[1].legend(
      handles=legend_patches,
      bbox_to_anchor=(1.01, 1),
      loc="upper left",
      borderaxespad=0.0,
  )

  plt.tight_layout()
  plt.savefig(args.output, dpi=300, bbox_inches="tight")
  plt.close()
  print(f"--> Saved wide output logo to: {args.output}")


def generate_residue_histograms(counts_df, seq_id, output_dir):
  """Generates and saves individual amino acid frequency histograms for each residue position.

  Filenames are formatted as: seq_<seq_id>_residue_<residue_num>.png
  """
  os.makedirs(output_dir, exist_ok=True)
  align_length = len(counts_df)

  print(
      f"--> Generating residue histograms for {align_length} positions into:"
      f" '{output_dir}'..."
  )

  plt.ioff()  # Turn off interactive plotting for speed

  for pos in range(align_length):
    counts = counts_df.iloc[pos]
    total_seqs = counts.sum()

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Plot Bar Chart
    bars = ax.bar(AMINO_ACIDS, counts.values, color="skyblue", edgecolor="navy")

    res_num = pos + 1  # 1-based biological indexing
    ax.set_title(
        f"Sequence {seq_id} | Residue Position {res_num} (N={total_seqs:,})",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xlabel("Amino Acid", fontsize=10, fontweight="bold")
    ax.set_ylabel("Count", fontsize=10, fontweight="bold")
    ax.set_ylim(0, max(counts.max() * 1.15, 1))

    # Annotate bars with non-zero counts
    for bar in bars:
      height = bar.get_height()
      if height > 0:
        ax.annotate(
            f"{int(height)}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    # File naming format: seq_1391_residue_1.png
    filename = f"seq_{seq_id}_residue_{res_num}.png"
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=150)
    plt.close(fig)

    if (res_num) % 50 == 0 or res_num == align_length:
      print(f"    Saved histograms up to residue position {res_num}/{align_length}")

  print(f"--> Done! All residue histograms saved to '{output_dir}'.")


def main():
  parser = argparse.ArgumentParser(
      description="Process giant A3M files into logos and per-residue histograms."
  )
  parser.add_argument(
      "--input", "-i", required=True, help="Path to input A3M file"
  )
  parser.add_argument(
      "--output",
      "-o",
      default="msa_logo.png",
      help="Output image filename for Logomaker",
  )
  parser.add_argument(
      "--seq-id",
      default="X",
      help="Sequence ID tag for file naming (e.g., 1391)",
  )

  # Flags to select output mode
  parser.add_argument(
      "--logo", action="store_true", help="Generate Logomaker sequence logo"
  )
  parser.add_argument(
      "--histogram",
      action="store_true",
      help="Generate individual per-residue histograms",
  )
  parser.add_argument(
      "--hist-dir",
      default="histograms",
      help="Output directory for residue histograms (default: ./histograms)",
  )

  # Custom sizing controls for Logomaker
  parser.add_argument(
      "--chunk-size",
      "-c",
      type=int,
      default=60,
      help="Positions per row (default: 60)",
  )
  parser.add_argument(
      "--width-per-pos",
      type=float,
      default=0.35,
      help="Width in inches per alignment position (default: 0.35)",
  )
  parser.add_argument(
      "--row-height",
      type=float,
      default=2.5,
      help="Height in inches per logo row (default: 2.5)",
  )

  args = parser.parse_args()

  # If neither flag is explicitly set, default to running BOTH
  if not args.logo and not args.histogram:
    args.logo = True
    args.histogram = True

  counts_df, raw_prob_df, grouped_prob_df, align_length = stream_and_process_a3m(
      args.input
  )

  if args.logo:
    generate_logomaker_plot(raw_prob_df, grouped_prob_df, align_length, args)

  if args.histogram:
    generate_residue_histograms(counts_df, args.seq_id, args.hist_dir)


if __name__ == "__main__":
  main()