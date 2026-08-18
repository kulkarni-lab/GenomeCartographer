#!/usr/bin/env python3

"""
RUN as: conda activate hoxpull
python3 1.Hox-multi-extract.py -p peps/ -r ref_hox_peps.fasta -i /pfs/home/siddharth/miniconda3/envs/hoxpull/bin/iqtree3   -o HoxTree -t 29

find_hox_iqtree_multi.py

Identify putative Hox genes from MULTIPLE peptide FASTA files in ./peps/,
combine all candidates with reference Hox proteins, align sequences,
remove gap-only taxa, and reconstruct a phylogenetic tree using IQ-TREE2/3.

Dependencies:
    - BLAST+
    - MAFFT
    - IQ-TREE2 or IQ-TREE3

USAGE:
    python find_hox_iqtree_multi.py \
        -p ./peps \
        -r ref_hox_peps.fasta \
        -i /path/to/iqtree2 \
        -o output_dir



INPUT:
    ./peps/*.fa
    ./peps/*.faa
    ./peps/*.fasta

OUTPUT:
    output_dir/
        all_candidate_hox.fa
        hox_alignment.cleaned.fa
        iqtree_run/hox.treefile
"""

import argparse
import glob
import os
import subprocess
import sys
from collections import OrderedDict
from shutil import which


###############################################################################
# Helpers
###############################################################################

def run_cmd(cmd):

    print(f"\n[CMD] {' '.join(cmd)}\n")

    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"\n[ERROR] Command failed:\n{' '.join(cmd)}")
        sys.exit(1)


def check_program(prog):

    if which(prog) is None:
        print(f"[ERROR] Missing dependency: {prog}")
        sys.exit(1)


###############################################################################
# FASTA parsing
###############################################################################

def read_fasta(fasta_file):

    seqs = OrderedDict()

    current = None

    with open(fasta_file) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                current = line[1:].split()[0]
                seqs[current] = []

            else:
                seqs[current].append(line)

    return {
        k: "".join(v)
        for k, v in seqs.items()
    }


###############################################################################
# Extract matching sequences
###############################################################################

def extract_matching_sequences(query_fasta,
                               hit_ids,
                               output_handle,
                               species_prefix):

    seqs = read_fasta(query_fasta)

    count = 0

    for seqid, seq in seqs.items():

        if seqid in hit_ids:

            new_id = f"{species_prefix}|{seqid}"

            output_handle.write(f">{new_id}\n")
            output_handle.write(seq + "\n")

            count += 1

    return count


###############################################################################
# Remove gap-only sequences
###############################################################################

def remove_gap_only_sequences(alignment_file, cleaned_file):

    seqs = OrderedDict()

    current = None

    with open(alignment_file) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                current = line
                seqs[current] = []

            else:
                seqs[current].append(line)

    removed = []

    with open(cleaned_file, "w") as out:

        for header, seq_lines in seqs.items():

            seq = "".join(seq_lines)

            nongap = (
                seq.replace("-", "")
                   .replace("?", "")
                   .replace("N", "")
                   .replace("X", "")
                   .replace("*", "")
                   .strip()
            )

            if len(nongap) == 0:

                removed.append(header[1:])
                continue

            out.write(header + "\n")
            out.write(seq + "\n")

    if removed:

        print("\n[WARNING] Removed gap-only sequences:")

        for r in removed:
            print(f"  - {r}")

    return removed


###############################################################################
# Main
###############################################################################

def main():

    parser = argparse.ArgumentParser(
        description="Find Hox genes from multiple peptide FASTAs."
    )

    parser.add_argument(
        "-p",
        "--pepdir",
        required=True,
        help="Directory containing peptide FASTA files"
    )

    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Reference Hox peptide FASTA"
    )

    parser.add_argument(
        "-i",
        "--iqtree",
        required=True,
        help="Path to IQ-TREE executable"
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Threads"
    )

    parser.add_argument(
        "-e",
        "--evalue",
        default="1e-5",
        help="BLAST e-value"
    )

    parser.add_argument(
        "-m",
        "--max_targets",
        type=int,
        default=5,
        help="Max target sequences"
    )

    args = parser.parse_args()

    ###########################################################################
    # Dependency checks
    ###########################################################################

    check_program("makeblastdb")
    check_program("blastp")
    check_program("mafft")

    if not os.path.isfile(args.iqtree):

        print(f"[ERROR] IQ-TREE executable not found: {args.iqtree}")
        sys.exit(1)

    ###########################################################################
    # Output directories
    ###########################################################################

    os.makedirs(args.outdir, exist_ok=True)

    tmpdir = os.path.join(args.outdir, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    iqdir = os.path.join(args.outdir, "iqtree_run")
    os.makedirs(iqdir, exist_ok=True)

    ###########################################################################
    # Collect peptide FASTAs
    ###########################################################################

    fasta_files = []

    extensions = ["*.fa", "*.faa", "*.fasta"]

    for ext in extensions:
        fasta_files.extend(
            glob.glob(os.path.join(args.pepdir, ext))
        )

    fasta_files = sorted(fasta_files)

    if len(fasta_files) == 0:

        print("[ERROR] No FASTA files found.")
        sys.exit(1)

    print("\n[INFO] Found peptide FASTAs:")

    for f in fasta_files:
        print(f"  - {f}")

    ###########################################################################
    # Combined candidate FASTA
    ###########################################################################

    combined_fa = os.path.join(
        args.outdir,
        "all_candidate_hox.fa"
    )

    with open(combined_fa, "w") as combined_out:

        #######################################################################
        # Write reference sequences first
        #######################################################################

        with open(args.reference) as ref:
            combined_out.write(ref.read())

        #######################################################################
        # Process each peptide FASTA
        #######################################################################

        for query_fasta in fasta_files:

            species = os.path.basename(query_fasta).split(".")[0]

            print(f"\n================================================")
            print(f"[INFO] Processing: {species}")
            print("================================================")

            db_prefix = os.path.join(
                tmpdir,
                f"{species}_db"
            )

            ###################################################################
            # Build BLAST DB
            ###################################################################

            run_cmd([
                "makeblastdb",
                "-in", query_fasta,
                "-dbtype", "prot",
                "-out", db_prefix
            ])

            ###################################################################
            # BLAST search
            ###################################################################

            blast_out = os.path.join(
                tmpdir,
                f"{species}.blast.tsv"
            )

            run_cmd([
                "blastp",
                "-query", args.reference,
                "-db", db_prefix,
                "-evalue", str(args.evalue),
                "-max_target_seqs", str(args.max_targets),
                "-outfmt",
                "6 qseqid sseqid pident length evalue bitscore",
                "-num_threads", str(args.threads),
                "-out", blast_out
            ])

            ###################################################################
            # Collect hit IDs
            ###################################################################

            hit_ids = set()

            with open(blast_out) as f:

                for line in f:

                    fields = line.strip().split("\t")

                    if len(fields) >= 2:
                        hit_ids.add(fields[1])

            print(f"[INFO] Candidate hits found: {len(hit_ids)}")

            ###################################################################
            # Extract candidate sequences
            ###################################################################

            n = extract_matching_sequences(
                query_fasta,
                hit_ids,
                combined_out,
                species
            )

            print(f"[INFO] Sequences added: {n}")

    ###########################################################################
    # MAFFT alignment
    ###########################################################################

    print("\n================================================")
    print("[INFO] Running MAFFT alignment")
    print("================================================")

    aln_file = os.path.join(
        args.outdir,
        "hox_alignment.fa"
    )

    with open(aln_file, "w") as out:

        mafft_cmd = [
            "mafft",
            "--auto",
            "--thread", str(args.threads),
            combined_fa
        ]

        subprocess.run(
            mafft_cmd,
            stdout=out,
            check=True
        )

    ###########################################################################
    # Remove gap-only sequences
    ###########################################################################

    print("\n================================================")
    print("[INFO] Removing gap-only sequences")
    print("================================================")

    clean_aln = os.path.join(
        args.outdir,
        "hox_alignment.cleaned.fa"
    )

    remove_gap_only_sequences(
        aln_file,
        clean_aln
    )

    ###########################################################################
    # IQ-TREE
    ###########################################################################

    print("\n================================================")
    print("[INFO] Running IQ-TREE")
    print("================================================")

    prefix = os.path.join(iqdir, "hox")

    run_cmd([
        args.iqtree,
        "-s", clean_aln,
        "-m", "Q.INSECT+F+G4",
#        "-bb", "1000",
#        "-alrt", "1000",
        "-nt", str(args.threads),
        "-pre", prefix
    ])

    ###########################################################################
    # Done
    ###########################################################################

    treefile = prefix + ".treefile"

    print("\n================================================")
    print("[DONE]")
    print("================================================")
    print(f"Combined sequences : {combined_fa}")
    print(f"Alignment          : {clean_aln}")
    print(f"Tree               : {treefile}")
    print("================================================\n")


###############################################################################
# Run
###############################################################################

if __name__ == "__main__":
    main()
