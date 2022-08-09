from ast import arg
import os
import sys
import yaml
import argparse
from collections import defaultdict


def network_detection(info_file):
    cluster_info = ''
    if os.path.exists(info_file):
        with open(info_file, 'r') as f:
            all_lines = f.readlines()
            header = all_lines[0]
            if 'No cluster' in header:
                cluster_info = 'False'
            else:
                cluster_info = 'True'
    else:     
        print('[{}] not exist!'.format(info_file))
        sys.exit(1)
    return cluster_info


def file_exist(file_list):
    for each_file in file_list:
        if not os.path.exists(each_file):
            print('[{}] not exist!'.format(each_file))
            sys.exit(1)
        else:
            pass


def generate_params():
    files = []
    kraken_dir = os.path.join(input_dir, '2.MTBC_identification')
    kraken = os.path.join(kraken_dir, 'kraken.result')
    p_dict['kraken'] = kraken
    files.append(kraken)
    qc_dir = os.path.join(input_dir, '1.Quality_control')
    seq_qc_dir = os.path.join(qc_dir, 'Sequence_QC')
    before_seq_qc_dir = os.path.join(seq_qc_dir, 'Before_trimming')
    after_seq_qc_dir = os.path.join(seq_qc_dir, 'After_trimming')
    seq_before_multiqc = '../../1.Quality_control/Sequence_QC/Before_trimming/MultiQC_report.html'
    seq_after_multiqc = '../../1.Quality_control/Sequence_QC/After_trimming/MultiQC_report.html'
    p_dict['seq_before_multiqc'] = seq_before_multiqc
    p_dict['seq_after_multiqc'] = seq_after_multiqc
    files.append(os.path.join(before_seq_qc_dir, 'MultiQC_report.html'))
    files.append(os.path.join(after_seq_qc_dir, 'MultiQC_report.html'))
    bam_qc_dir = os.path.join(qc_dir, 'Alignment_QC')
    bam_multiqc = '../../1.Quality_control/Alignment_QC/MultiQC_report.html'
    p_dict['bam_multiqc'] = bam_multiqc
    files.append(os.path.join(bam_qc_dir, 'MultiQC_report.html'))
    dist_dir = os.path.join(input_dir, '4.SNP_distance')
    snp_diff_matrix = os.path.join(dist_dir, '{}_pairwise_distance_matrix.txt'.format(prefix))
    p_dict['snp_diff_matrix'] = snp_diff_matrix
    files.append(snp_diff_matrix)
    snp_diff_heatmap = os.path.join(dist_dir, '{}_pairwise_distance_heatmap.png'.format(prefix))
    p_dict['snp_diff_heatmap'] = snp_diff_heatmap
    files.append(snp_diff_heatmap)
    snp_diff_histogram = os.path.join(dist_dir, '{}_pairwise_distance_distribution.png'.format(prefix))
    p_dict['snp_diff_histogram'] = snp_diff_histogram
    files.append(snp_diff_histogram)
    result_dir = os.path.join(input_dir, '5.Transmission_cluster')
    cluster_out_dir = os.path.join(result_dir, cluster_dir)
    cluster_table = os.path.join(cluster_out_dir, 'transmission_clusters.txt')
    files.append(cluster_table)
    p_dict['cluster_table'] = cluster_table
    clustering_visualization = os.path.join(cluster_out_dir, 'clustering_visualization.png')
    p_dict['clustering_visualization'] = clustering_visualization
    files.append(clustering_visualization)
    cluster_member_visualization = os.path.join(cluster_out_dir, 'cluster_member_visualization.png')
    p_dict['cluster_member_visualization'] = cluster_member_visualization
    files.append(cluster_member_visualization)
    network_info = os.path.join(cluster_out_dir, 'transmission_detection_run_info.txt')
    show_network = network_detection(network_info)
    p_dict['show_network'] = show_network
    example_dir = os.path.join(cluster_out_dir, 'cluster_1')
    network_image = os.path.join(example_dir, 'transmission_network.png')
    p_dict['network_image'] = network_image
    if show_network == 'True':
        files.append(network_image)
    risk_dir = os.path.join(input_dir, '6.Risk_factor')
    risk_out_dir = os.path.join(risk_dir, cluster_dir)
    risk_file = os.path.join(risk_out_dir, 'characteristic_data.txt')
    p_dict['risk_file'] = risk_file
    show_risk = 'False'
    if os.path.exists(risk_file):
        show_risk = 'True'
    p_dict['show_risk'] = show_risk
    file_exist(files)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", dest="a", help="result directory", required=True)
    parser.add_argument("-b", dest="b", help="kraken cutoff", required=True)
    parser.add_argument("-m", dest="m", help="non-MTBC reads filtering", required=True)
    parser.add_argument("-c", dest="c", help="mapping quality threshold", required=True)
    parser.add_argument("-d", dest="d", help="allele frequency threshold", required=True)
    parser.add_argument("-e", dest="e", help="depth threshold", required=True)
    parser.add_argument("-f", dest="f", help="clustering method", required=True)
    parser.add_argument("-g", dest="g", help="SNP threshold", required=True)
    parser.add_argument("-i", dest="i", help="transmission threshold", required=True)
    parser.add_argument("-j", dest="j", help="clustering directory", required=True)
    parser.add_argument("-p", dest="p", help="prefix of output file", required=True)
    args = parser.parse_args()
    # print('Generating parameters for R markdown report.')
    p_dict = defaultdict()
    input_dir = args.a
    kraken_cutoff = args.b
    p_dict['kraken_cutoff'] = kraken_cutoff
    mtbc_only = args.m
    p_dict['mtbc_only'] = mtbc_only
    mapping_quality = args.c
    p_dict['mapping_quality'] = mapping_quality
    allele_frequency = args.d
    p_dict['allele_frequency'] = allele_frequency
    mapping_depth = args.e
    p_dict['mapping_depth'] = mapping_depth
    method = args.f
    p_dict['cluster_method'] = method
    cluster_dir = args.j
    snp_cutoff = args.g
    p_dict['snp_cutoff'] = snp_cutoff
    trans_cutoff = args.i
    p_dict['trans_cutoff'] = trans_cutoff
    prefix = args.p
    p_dict['prefix'] = prefix
    generate_params()
    output_dir = os.path.join(input_dir, '7.Summary_report')
    os.makedirs(output_dir, exist_ok=True)
    output_method_dir = os.path.join(output_dir, cluster_dir)
    os.makedirs(output_method_dir, exist_ok=True)
    result_yaml = os.path.join(output_method_dir, 'report.yaml')
    with open(result_yaml, 'w') as outfile:
        yaml.dump(p_dict, outfile, default_flow_style=False)
