#############################################################
# Automatic transmission detection for all genomic clusters #
#############################################################
import os
import sys
import subprocess as sub
import argparse
from collections import defaultdict


def run_cmd(command, debug):
    """
    Linux system bash runner
    :param command: command line
    :param debug: full output
    :return: running status
    """
    if not debug:
        ret = sub.run(command, shell=True, stderr=sub.PIPE, stdout=sub.PIPE)
    else:
        ret = sub.run(command, shell=True)
    if ret.returncode == 0:
        failure = False
    else:
        failure = True
    return failure


def fetch_cluster(cluster_file, out_dir, method):
    """
    Fetching all clusters and reorder the clusters by their members.
    :param method: transmission clustering method.
    :param cluster_file: a cluster file from transcluster export.
    :param out_dir: the output directory for all transmission analyses.
    :return: all_cluster_dict: a dict with all cluster ids and members.
    :export: transmission_clusters.txt
    """
    all_cluster_dict = defaultdict(list)
    trans_result_file = os.path.join(out_dir, 'transmission_clusters.txt')
    result_lines = 'sample\ttransmission_cluster\tclustering_type\n'
    tmp_dict = defaultdict(int)
    if method == 'SNP':
        group_dict = defaultdict(list)
        with open(cluster_file, 'r') as f1:
            i = 0
            for each_line in f1.readlines()[1:]:
                c_list = each_line.strip('\r|\n|,').split(',')
                sample = c_list[0]
                tmp_group = c_list[1]
                if tmp_group not in ['ungrouped', 'proxy']:
                    group_dict[tmp_group].append(sample)
                else:
                    i += 1
                    group_dict['ungrouped' + str(i)].append(sample)
        for c, c_list in group_dict.items():
            tmp_dict[tuple(c_list)] = len(c_list)
    else:
        with open(cluster_file, 'r') as f1:
            for each_line in f1.readlines()[1:]:
                c_list = each_line.strip('\r|\n|,').split(',')
                tmp_dict[tuple(c_list)] = len(c_list)
    sorted_list = sorted(tmp_dict.items(), key=lambda d: d[1], reverse=True)
    i = 1
    clustered_sample_num = 0
    all_sample_num = 0
    cluster_num = 0
    for each_item in sorted_list:
        s_list = each_item[0]
        if len(s_list) > 1:
            c = 'Clustered'
            cluster_num += 1
        else:
            c = 'Unique'
        for s in s_list:
            all_sample_num += 1
            if c == 'Clustered':
                clustered_sample_num += 1
                result_lines += '{}\t{}\t{}\n'.format(s, i, c)
            else:
                result_lines += '{}\t\t{}\n'.format(s, c)
        all_cluster_dict[i] = s_list
        i += 1
    cluster_rate = '{:.1%}'.format(clustered_sample_num / all_sample_num)
    with open(trans_result_file, 'w') as f2:
        header = '# {} ({}) of the {} MTB samples were grouped into {} ' \
                 'transmission clusters.'.format(str(clustered_sample_num), 
                                                 str(cluster_rate),
                                                 str(all_sample_num), 
                                                 str(cluster_num))
        print(header)
        f2.write(header + '\n' + result_lines)
    if cluster_num == 0:
        print('[*] No transmission cluster found.')
    return all_cluster_dict, cluster_num


def check_files(output_dir):
    """
    Checking the output files from seqtrack_analysis.R
    :param output_dir: the directory where seqtrack_result.txt, transmission_network.[pdf/png] located.
    :return: all files exists or not.
    """
    res_file = os.path.join(output_dir, 'seqtrack_result.txt')
    image_pdf = os.path.join(output_dir, 'transmission_network.pdf')
    image_png = os.path.join(output_dir, 'transmission_network.png')
    export_ok = True
    for each_file in [res_file, image_pdf, image_png]:
        if not os.path.exists(each_file):
            export_ok = False
            break
    return export_ok


def run_seqtrack(dist_file, sample_file, date_file, output_dir, use_coord):
    """
    Running seqtrack_analysis.R (locates in the same directory of this script) to do transmission detection.
    :param use_coord: using coordinate data.
    :param dist_file: the pairwise distance file from PANPASCO pipeline.
    :param sample_file: the tsv format file containing all samples of a cluster.
    :param date_file: the tsv format file containing all samples collection dates.
    :param output_dir: the output directory for each transmission cluster with more than 3 samples..
    :return: None
    """
    src_dir = os.path.split(os.path.realpath(__file__))[0]
    r_script = os.path.join(src_dir, "seqtrack_analysis.R")

    if not os.path.exists(r_script):
        print("Error: the seqtrack R script [{}] is not exists.".format(
            r_script))
        sys.exit(1)

    if use_coord in ['TRUE', 'True', 'true', 'T', 't']:
        print('Using longitude and latitude information data.')
    log_file = os.path.join(output_dir, "seqtrack.log")
    r_cmd = 'Rscript {r} -i {i} -s {s} -m {m} -c {c} -o {o} 2> {l}'.format(
        r=r_script,
        i=dist_file,
        s=sample_file,
        m=date_file,
        c=use_coord,
        o=output_dir,
        l=log_file)

    failure = run_cmd(r_cmd, debug=True)
    if failure:
        print("[*] Error: run seqtrack R script failed.")
        sys.exit(1)


def transmission_detection(all_cluster_dict, cluster_num, out_dir, dist_file,
                           network_inference, date_file, use_coord, xmin):
    """
    Running transmission detection for every valid clusters (with more than 3 samples).
    :param cluster_num: detected clusters number.
    :param use_coord: using coordinate data.
    :param all_cluster_dict: the dict with all clusters information.
    :param out_dir: the output directory for all transmission analyses.
    :param dist_file: the pairwise distance file from PANPASCO pipeline.
    :param date_file: the tsv format file containing all samples collection dates.
    :param xmin: the minimum number of samples in a cluster to be considered valid cluster
    :return: None
    """
    check_lines = ''
    if not network_inference:
        print("[*] Skip transmission network inference.")
        check_lines = 'No cluster network inference.\n'
    elif cluster_num > 0:
        print(
            '=> Using SeqTrack to infer transmission events for all clusters with at least 4 samples.'
        )
        i = 0
        for cluster, samples in all_cluster_dict.items():
            if len(samples) >= xmin:
                i += 1
                print('==> Cluster {} ... '.format(cluster), end='')
                cluster_dir = os.path.join(out_dir,
                                           'cluster_{}'.format(cluster))
                os.makedirs(cluster_dir, exist_ok=True)
                sample_file = os.path.join(cluster_dir, 'samples.txt')
                with open(sample_file, 'w') as out:
                    for each_sample in samples:
                        out.write(each_sample + '\n')
                run_seqtrack(dist_file, sample_file, date_file, cluster_dir,
                             use_coord)
                check_ok = check_files(cluster_dir)
                if check_ok:
                    result = 'cluster_{}\tOK\n'.format(cluster)
                else:
                    result = 'cluster_{}\tFailed\n'.format(cluster)
                check_lines += result
                print('Done')
        if i == 0:
            check_lines += 'No cluster with at least 3 samples.\n'
            print(
                "Warning: no cluster with more than 3 samples for transmission detection."
            )
    else:
        check_lines = 'No cluster found in this study.\n'
    result_file = os.path.join(out_dir, 'transmission_detection_run_info.txt')
    with open(result_file, 'w') as out:
        out.write(check_lines)


def main():
    """
    Main function for transmission detection.
    :return: None
    """
    dist_file = args.distance
    cluster_file = args.cluster
    date_file = args.date
    out_dir = args.out_dir
    use_coord = args.coord
    method = args.method
    network = args.network
    xmin = args.xmin
    all_cluster_dict, cluster_num = fetch_cluster(cluster_file, 
                                                  out_dir,
                                                  method)
    network_inference = False
    if network in ['TRUE', 'True', 'true', 'T', 't']:
        network_inference = True
    transmission_detection(all_cluster_dict, cluster_num, out_dir, dist_file,
                           network_inference, date_file, use_coord, xmin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",
                        "--cluster",
                        dest="cluster",
                        help="Input transcluster result file.",
                        required=True)
    parser.add_argument("-i",
                        "--distance",
                        dest="distance",
                        help="Input pairwise distance file.",
                        required=True)
    parser.add_argument("-d",
                        "--date",
                        dest="date",
                        help="Input sample collection date file.",
                        required=True)
    parser.add_argument("-o",
                        "--output",
                        dest="out_dir",
                        help="Transmission cluster output directory.",
                        required=True)
    parser.add_argument("-l",
                        "--coord",
                        dest="coord",
                        help="Using coordinate information.",
                        required=True)
    parser.add_argument("-m",
                        "--method",
                        dest="method",
                        help="Transmission clustering method.",
                        required=True)
    parser.add_argument("-n",
                        "--network",
                        dest="network",
                        help="Network inference.",
                        required=True)
    parser.add_argument("-x",
                        "--xmin",
                        dest="xmin",
                        help="Minimum cluster size.",
                        required=True)
    args = parser.parse_args()
    main()
