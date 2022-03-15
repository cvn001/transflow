import argparse
import re
# This script writes all positions with GT ./.   (uncalled by GATK)  to a BED file


def main():
    global args
    # index_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    with open(args.vcffile, "r") as inf, open(args.vcfoutfile, "w") as outvcf, open(args.qualoutfile, "w") as outqual:
        for line in inf:
            if not line.startswith("#"):
                fields = line.split("\t")
                genome = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                gt = fields[9]
                length = (len(ref) - len(alt))
                sec_pos = pos + length if length > 0 else pos
                m = re.match(r"(\d/\d)", gt)
                if m is not None:
                    gt = m.group(1)
                    if gt == "1/1" or gt == "0/1" or gt == "1/0":
                        if alt != "*":
                            outvcf.write(line)
                    elif gt == "1/2" or gt == "1/3":
                        outqual.write("\t".join([genome, str(pos - 1), str(pos)]) + "\n")
                    elif gt != "0/0":
                        print("weird GT found!!!!!!!:" + line)
                else:
                    m = re.match(r"(\./\.)", gt)
                    if m is not None:
                        outqual.write("\t".join([genome, str(pos - 1), str(sec_pos)]) + "\n")
                    else:
                        print("NO GT found!!!!!!!:" + line)
            else:
                outvcf.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="vcffile", help="Input VCF file.", required=True)
    parser.add_argument("--output_vcf", dest="vcfoutfile", help="Output VCF file.", required=True)
    parser.add_argument("--output_bed", dest="qualoutfile", help="Output BED file with uncalled regions.", required=True)
    # parser.add_argument("-g", "--gq", dest="gq", help="Genotype Quality threshold ", required=True, type=float)
    args = parser.parse_args()
    main()
