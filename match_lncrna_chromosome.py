"""
Usage:
python3 chr_lnc_match.py <chromosome_key_to_match> <path_to_nucleotide_file> \
    <path_to_lncRNA_mapping_file>

ex:
python3 chrX DMS-MaPseq_hg38_Nature_Methods.bed export.bed
"""

from argparse import ArgumentParser
from pathlib import Path
import re
import sys

fp = Path(
    "DMS-MaPseq_hg38_Nature_Methods_2017_DMS_invivo_"
    "score_both_score_DMSMapSeq_2017.score.bed"
)
chr_map = Path("chromosome_lncrna_match/export.bed")


NUC_RE = re.compile(r"(?P<chr>chr\S+)\s+(?P<nuc1>\d+)\s+(?P<nuc2>\d+)")
CHR_RE = re.compile(
    r"(?P<chr>chr\S+)\s+(?P<start>\d+)\s+(?P<end>\d+)\s+(?P<name>\S+)\s+"
)


def get_nuc1s(fp, chr):
    with open(fp, "r") as f:
        counter = 0
        nuc1s = []
        for line in f:
            if line.startswith(chr):
                counter += 1
                match = NUC_RE.match(line)
                nuc1s.append(int(match.groupdict()["nuc1"]))

    if len(nuc1s) == 0:
        print(f"No matching chromosome label: {chr}")
    return nuc1s


def get_chr_dict(filepath, chr):
    with open(filepath, "r") as f:
        lnc_dict = dict()
        for line in f:
            if line.startswith(chr):
                match = CHR_RE.match(line).groupdict()
                lnc_dict[(int(match["start"]), int(match["end"]))] = match["name"]

    if len(lnc_dict) == 0:
        print(f"No matching chromosome label {chr}")
    return lnc_dict


def main(chr, nuc_fp, chr_fp):
    nuc1s = get_nuc1s(nuc_fp, chr)
    chr_map = get_chr_dict(chr_fp, chr)

    lncrnas = set()
    for nuc in nuc1s:
        # print(nuc)
        for key in chr_map:
            # print(key)
            if key[0] < nuc < key[1]:
                lncrnas.add(chr_map[key])
    outfile = Path(f"{chr}_lncrna_matches.txt")
    while outfile.exists():
        outfile = Path(f"{outfile}_cp")
    with open(outfile, "w") as f:
        f.writelines([x + "\n" for x in lncrnas])
    print(lncrnas)
    print(f"lncRNA names saved to {outfile}")


if __name__ == "__main__":
    argv = sys.argv[1:]
    parser = ArgumentParser()
    parser.add_argument(
        "chromosome",
        type=str,
        help="Chromosome code to be matched, e.g., 'chrX'",
    )
    parser.add_argument(
        "nuc_fp",
        type=str,
        help="Filepath for nucleotide file (line looks like 'chrY	57130217"
        "	57130218	NA	0.006	+')",
    )
    parser.add_argument(
        "chr_fp",
        type=str,
        help="Filepath for chromosome-to-lncRNA mapping file (export.bed)",
    )
    args = parser.parse_args(argv)
    chromosome = args.chromosome
    nuc_fp = args.nuc_fp
    chr_fp = args.chr_fp
    main(chromosome, nuc_fp, chr_fp)
