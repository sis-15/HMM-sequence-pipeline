import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logomaker
from Bio import SeqIO

# Single-character tokens for Logomaker compatibility
GROUPS_CONFIG = {
    # We can change these as desired - this is just a quick example
    'H': ('Hydrophobic (LVIMC)', ['L', 'V', 'I', 'M', 'C'], 'black'),
    'Ω': ('Aromatic (FWY)',      ['F', 'W', 'Y'],           'purple'),
    '-': ('Acidic/Polar (EDNQ)', ['E', 'D', 'N', 'Q'],      'red'),
    '+': ('Basic (KRH)',         ['K', 'R', 'H'],           'blue'),
    'o': ('Small/Polar (STAG)',  ['S', 'T', 'A', 'G'],      'green'),
    'P': ('Proline (P)',         ['P'],                     'orange')
}

def stream_and_process_a3m(a3m_path):
    """
    Streams an A3M file line-by-line to handle giant files without memory issues.
    Calculates amino acid counts directly without storing all sequences in RAM.
    """
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
    
    counts_matrix = None
    seq_count = 0
    align_length = 0

    print(f"--> Streaming and parsing A3M file: {a3m_path}")
    
    # Stream fasta records one by one
    for record in SeqIO.parse(a3m_path, "fasta"):
        # Remove insertion states (lowercase) and deletion dots
        clean_seq = re.sub(r'[a-z.]', '', str(record.seq))
        
        # Initialize count matrix based on the first sequence length
        if counts_matrix is None:
            align_length = len(clean_seq)
            counts_matrix = np.zeros((align_length, len(amino_acids)), dtype=int)
            print(f"--> Alignment length determined by query: {align_length} positions")

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

    # Build DataFrames
    counts_df = pd.DataFrame(counts_matrix, columns=amino_acids)
    total_per_col = counts_df.sum(axis=1).replace(0, 1)
    prob_df = counts_df.div(total_per_col, axis=0)

    # Build Grouped Matrix
    grouped_df = pd.DataFrame(index=prob_df.index)
    for char_code, (label, aa_list, color) in GROUPS_CONFIG.items():
        grouped_df[char_code] = prob_df[aa_list].sum(axis=1)

    return prob_df, grouped_df, align_length

def main():
    parser = argparse.ArgumentParser(description="Process giant A3M files into multi-row wrapped sequence logos.")
    parser.add_argument("--input", "-i", required=True, help="Path to input A3M file")
    parser.add_argument("--output", "-o", default="msa_logo.png", help="Output image filename")
    parser.add_argument("--chunk-size", "-c", type=int, default=60, help="Positions per row (default: 60)")
    # ADDED: Custom sizing controls
    parser.add_argument("--width-per-pos", type=float, default=0.35, help="Width in inches per alignment position (default: 0.35)")
    parser.add_argument("--row-height", type=float, default=2.5, help="Height in inches per logo row (default: 2.5)")
    args = parser.parse_args()

    raw_prob_df, grouped_prob_df, align_length = stream_and_process_a3m(args.input)

    # Calculate number of wrapped rows needed for long proteins
    num_chunks = int(np.ceil(align_length / args.chunk_size))
    
    # DYNAMIC SIZING CALCULATIONS:
    # Scale total width by the number of positions per chunk so it never smushes
    fig_width = max(12, args.chunk_size * args.width_per_pos)
    fig_height = num_chunks * 2 * args.row_height

    # Create subplots with scaled dimensions
    fig, axes = plt.subplots(
        nrows=num_chunks * 2, 
        ncols=1, 
        figsize=(fig_width, fig_height), 
        sharey=True
    )
    
    # Ensure axes is an array even if num_chunks == 1
    if num_chunks * 2 == 2:
        axes = np.array([axes[0], axes[1]])

    group_colors = {code: config[2] for code, config in GROUPS_CONFIG.items()}

    for i in range(num_chunks):
        start_pos = i * args.chunk_size
        end_pos = min((i + 1) * args.chunk_size, align_length)
        
        ax_top = axes[i * 2]
        ax_bottom = axes[i * 2 + 1]

        # Slice data for current window
        chunk_raw = raw_prob_df.iloc[start_pos:end_pos]
        chunk_grouped = grouped_prob_df.iloc[start_pos:end_pos]

        # Plot Raw AA Logo (vpad=0.05 ensures characters take up full width)
        logomaker.Logo(chunk_raw, ax=ax_top, color_scheme='chemistry', vpad=0.05)
        ax_top.set_title(f"Positions {start_pos} - {end_pos} (Standard AAs)", fontsize=11, fontweight='bold')
        ax_top.set_ylabel("Frequency")

        # Plot Grouped Logo
        logomaker.Logo(chunk_grouped, ax=ax_bottom, color_scheme=group_colors, vpad=0.05)
        ax_bottom.set_title(f"Positions {start_pos} - {end_pos} (Grouped Properties)", fontsize=11, fontweight='bold')
        ax_bottom.set_ylabel("Frequency")

        for ax in [ax_top, ax_bottom]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # Ensure x-axis ticks align tightly to integer position boundaries
            ax.set_xlim([start_pos - 0.5, end_pos - 0.5])

    # Add legend to the top-right corner outside the plot area
    legend_patches = [
        mpatches.Patch(color=config[2], label=f"'{code}': {config[0]}")
        for code, config in GROUPS_CONFIG.items()
    ]
    axes[1].legend(handles=legend_patches, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"--> Saved wide output logo to: {args.output}")

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()