########################################################################
# Extract the percentage of reads mapping to MTBC from a kraken report #
########################################################################

import sys

kraken_report = sys.argv[1]
output = sys.argv[2]
sample_name = sys.argv[3]

result_lines = ''
with open(kraken_report, 'r') as f1:
    all_lines = f1.readlines()
    species = 'None MTBC'
    percentage = 0
    for each_line in all_lines:
        if 'Mycobacterium tuberculosis complex' in each_line:
            m_list = each_line.strip().split('\t')
            species = 'MTBC'
            percentage = m_list[0].strip()
    result_lines += '{}\t{}\t{}\n'.format(sample_name, percentage, species)
with open(output, 'w') as f2:
    f2.write(result_lines)
