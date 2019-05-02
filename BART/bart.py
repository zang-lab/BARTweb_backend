import os,sys
import argparse
import yaml # requirements: pyyaml

import rpregress_adaptive_lasso
import enhancerIdentifier_adaptive_lasso
import AUCcalc_with_json


def marge_bart_pipeline(args, conf_data):
    if args.species == 'hg38':
        sys.stdout.write("Start prediction on hg38...\n")
    elif args.species == 'mm10':
        sys.stdout.write("Start prediction on mm10...\n")

    # load config data
    rp_file = conf_data[args.species]['H3K27ac_rp']
    rpkm_file = conf_data[args.species]['H3K27ac_RPKM']
    sym_file = conf_data[args.species]['TSS_bed']
    annotation_file = conf_data[args.species]['annotation_file']
    tf_json = conf_data[args.species]['TF_json']
    overlap_json = conf_data[args.species]['overlap_json']
    UAC_info = conf_data[args.species]['UAC_info']

    if args.input_type == 'genelist':
        sys.stdout.write("Do adaptive lasso to select H3K27ac samples...\n")

        rpregress_adaptive_lasso.main(args.species, rp_file, args.input_file, sym_file, args.output_name, 'target', 20, False, True, 'Gene_Only', True, False, annotation_file)
        regression_info = args.output_name + '_adaptive_lasso_Info.txt'
        if not os.path.exists(regression_info):
            sys.stderr.write("adaptive lasso information does not exist! \n")
            sys.exit(0)

        sys.stdout.write("generate cis-regulatory profile...\n")
        enhancerIdentifier_adaptive_lasso.main(regression_info, args.input_file, args.output_name, args.species, UAC_info, rpkm_file)
        enhancer_info = args.output_name + '_enhancer_prediction_lasso.txt'
        if not os.path.exists(enhancer_info):
            sys.stderr.write("enhancer information does not exist! \n")
            sys.exit(0)

        sys.stdout.write("fast BART!...\n")
        norm_file = conf_data[args.species]['MSigDB_norm']
        auc_file = args.output_name + '_auc.txt'
        result_file = args.output_name + '_bart_results.txt'
        AUCcalc_with_json.cal_auc(enhancer_info, tf_json, overlap_json, auc_file, result_file, norm_file)

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.realpath(__file__))
    conf_path = os.path.join(script_dir, 'config.yaml')
    if not os.path.exists(conf_path):
        sys.stderr.write("config.yaml does not exist..\n")
        sys.exit(0)

    sys.stdout.write("loading configuration file...\n")
    with open(conf_path, 'r') as fconf:
        try:
            conf_data = yaml.safe_load(fconf)
        except yaml.YAMLError as e:
            sys.stderr.write("config.yaml does not exist.. \n")
            sys.stderr.write(e)
            sys.exit(0)

    try:
        parser = argparse.ArgumentParser(description="""Adaptive lasso + json storage on TF prediction.""")
        parser.add_argument( '-s','--species', dest='species', required = True, type = str, choices=['hg38','mm10'], help='genome species')
        parser.add_argument( '-t','--input_type', dest='input_type', required = True, type = str, choices=['genelist', 'ChIP-seq'], help='input file type')
        parser.add_argument( '-i','--input_file', dest='input_file', required = True, type = str, help='input file')
        parser.add_argument( '-o','--output_name', dest='output_name', required = True, type = str, help='output prefix name')

        args = parser.parse_args()
        marge_bart_pipeline(args, conf_data)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted by user... bye ;) \n")
        sys.exit(0)




