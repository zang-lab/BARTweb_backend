import os,sys

# import from package
from bart2 import OptValidator
from bart2 import ReadCount
from bart2 import RPRegress, EnhancerIdentifier
from bart2 import AUCcalc, StatTest

script_dir = os.path.dirname(os.path.realpath(__file__))
ADAPTIVE_LASSO_MAXSAMPLES = 20 # TODO: should we fix it?

def bart(options):
    args = OptValidator.opt_validate(options)
 
    # create output directory
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except:
        sys.stderr.write('Output directory: {} could not be created. \n'.format(args.outdir))
        sys.exit(0)
    sys.stdout.write("Output directory will be {} \n".format(args.outdir))
    sys.stdout.write("Output file prefix will be {} \n".format(args.ofilename))

    if args.species == 'hg38':
        sys.stdout.write("Start prediction on hg38...\n")
    elif args.species == 'mm10':
        sys.stdout.write("Start prediction on mm10...\n")

    # bart geneset [-h] <-i genelist.txt> [--refseq] <-s species> [-t target] [-p processes] [--outdir] [options]
    if args.subcommand_name == 'geneset':
        '''
        Use adaptive lasso regression to select H3K27ac samples.

        RPRegress parameters:
        species, rp matrix, gene file, refseq TSS, output directory, target gene method, adaptive lasso max sapmle numbers, log transform, square transform, gene symbol or refseqID, gene symbol or not, separate by chrom, sample description file
        '''
        sys.stdout.write("Do adaptive lasso to select informative H3K27ac samples...\n")
        if options.refseq:
            RPRegress.main(args.species, args.rp, args.infile, args.tss, args.ofilename, 'target', ADAPTIVE_LASSO_MAXSAMPLES, False, True, 'Gene_Response', False, False, args.desc)
        else:
            RPRegress.main(args.species, args.rp, args.infile, args.tss, args.ofilename, 'target', ADAPTIVE_LASSO_MAXSAMPLES, False, True, 'Gene_Only', True, False, args.desc)
        regression_info = args.ofilename + '_adaptive_lasso_Info.txt'
        if not os.path.exists(regression_info):
            sys.stderr.write("Error: selecting samples from H3K27ac compendium! \n")
            sys.exit(0)

        '''
        Generate cis-regulatory profile based on adaptive lasso model weights multiply H3K27ac samples RPKM signal.

        EnhancerIdentifier parameters:
        selected samples file, gene file, output directory, species, UDHS, rpkm matrix
        '''
        sys.stdout.write("Generate cis-regulatory profile...\n")
        EnhancerIdentifier.main(regression_info, args.ofilename, args.dhsfile, args.rpkm)
        enhancer_profile = args.ofilename + '_enhancer_prediction_lasso.txt'
        if not os.path.exists(enhancer_profile):
            sys.stderr.write("Error: generating enhancer profile! \n")
            sys.exit(0)

        # get ranked score UDHS positions from enhancer profile
        positions = AUCcalc.get_position_list(enhancer_profile)

    # bart profile [-h] <-i ChIP-seq profile> <-f format> <-s species> [-t target] [-p processes] [--outdir] [options]
    elif args.subcommand_name == 'profile':
        sys.stdout.write('Start mapping the {} file...\n'.format(args.format.upper()))
        counting = ReadCount.read_count_on_DHS(args)
        # get ranked score UDHS positions from read count
        positions = sorted(counting.keys(),key=counting.get,reverse=True)


    '''
    Start using revised BART on calculating the AUC score for each TF ChIP-seq dataset
    '''
    sys.stdout.write("revised BART!...\n")
    sys.stdout.write('Prediction starts...\n\nRank all DHS...\n')
    positions = AUCcalc.get_position_list(enhancer_profile)
    tf_aucs, tf_index = AUCcalc.cal_auc(args, positions)

    stat_file = args.ofilename + '_bart_results.txt'
    StatTest.stat_test(tf_aucs, tf_index, stat_file, args.normfile)
    sys.stdout.write("BART job finished successfully!\n")

