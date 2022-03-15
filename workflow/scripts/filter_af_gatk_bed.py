import argparse
import re
# This script writes all positions with AF lower than threshold to a BED file

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
                info = fields[7]
                gt = fields[9]
                sec_pos = pos + (len(ref) - len(alt))
                af = 0
                m = re.match(r"\d/\d:(\d+(,\d+)+):", gt)
                if m is not None:
                    # GATK does weird stuff - DP = 0 but GT is 1/1 --> filter!
                    alts = [int(x) for x in m.group(2).split(",")[1:]]
                    counts = [int(x) for x in m.group(1).split(",")]
                    count_a = max(alts)
                    count_sum = sum(counts)
                    if count_sum > 0:
                        # GATK does weird stuff - DP = 0 but GT is 1/1 --> filter!
                        af = count_a / count_sum
                if af < args.af:
                    out.write("\t".join([genome, str(pos - 1), str(sec_pos)]) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="vcffile", help="Input VCF file.", required=True)
    parser.add_argument("-o", "--output", dest="outfile", help="Output BED file.", required=True)
    parser.add_argument("-a", "--af", dest="af", help="Allele Frequency threshold (inclusive) ",
                        required=True, type=float)

    args = parser.parse_args()
    main()
