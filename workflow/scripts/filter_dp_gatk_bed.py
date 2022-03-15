## This script writes all positions with AF lower than threshold to a BED file

import argparse
import re


def main():
    global args
    # index_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    with open(args.vcffile, "r") as inf, open(args.outfile, "w") as out:
        for line in inf:
            if not line.startswith("#"):
                fields = line.split("\t")
                genome = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                pos_filter = fields[6]
                info = fields[7]
                gt = fields[9]
                sec_pos = pos + (len(ref) - len(alt))
                dp = 0
                dp_sum = 0
                m = re.search(r"DP=(\d+);", info)
                if m is not None:
                    dp = float(m.group(1))
                if dp >= args.dp:
                    m = re.match(r"\d/\d:(\d+(,\d+)+):", gt)
                    if m is not None:
                        counts = [int(x) for x in m.group(1).split(",")]
                        dp_sum = sum(counts)
                if dp < args.dp or dp_sum < args.dp or pos_filter == "LowQual":
                    out.write("\t".join([genome, str(pos-1), str(sec_pos)]) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="vcffile", help="Input VCF file.", required=True)
    parser.add_argument("-o", "--output", dest="outfile", help="Output BED file.", required=True)
    parser.add_argument("-d", "--dp", dest="dp", help="Depth threshold, ", required=True, type=float)
    args = parser.parse_args()
    main()
