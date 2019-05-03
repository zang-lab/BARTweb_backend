import os,sys

from revised_bart import OptValidator
from revised_bart import ReadCount
from revised_bart import RPRegress, EnhancerIdentifier
from revised_bart import AUCcalc, StatTest

script_dir = os.path.dirname(os.path.realpath(__file__))
ADAPTIVE_LASSO_MAXSAMPLES = 20

def bart(options):
    args = OptValidator.opt_validate(options)
 
    # create output directory
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except:
        sys.exit('Output directory: {} could not be created.'.format(args.outdir))
    
    print("Output directory will be {}".format(args.ofilename))

    if args.species == 'hg38':
        sys.stdout.write("Start prediction on hg38...\n")
    elif args.species == 'mm10':
        sys.stdout.write("Start prediction on mm10...\n")

    # bart geneset [-h] <-i genelist.txt> [--refseq] <-s species> [-t target] [-p processes] [--outdir] [options]
    if args.subcommand_name == 'geneset':
        '''
        Use adaptive lasso regression to select H3K27ac samples.
        '''
        sys.stdout.write("Do adaptive lasso to select H3K27ac samples...\n")
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
        '''
        sys.stdout.write("Generate cis-regulatory profile...\n")
        EnhancerIdentifier.main(regression_info, args.infile, args.ofilename, args.species, args.dhsfile, args.rpkm)
        enhancer_profile = args.ofilename + '_enhancer_prediction_lasso.txt'
        if not os.path.exists(enhancer_profile):
            sys.stderr.write("Error: generating enhancer profile! \n")
            sys.exit(0)

        '''
        Start using revised BART on calculating the AUC score for each TF ChIP-seq dataset
        '''
        sys.stdout.write("revised BART!...\n")
        sys.stdout.write('Prediction starts...\n\nRank all DHS...\n')
        AUCcalc.cal_auc(enhancer_profile, args)

    # bart profile [-h] <-i ChIP-seq profile> <-f format> <-s species> [-t target] [-p processes] [--outdir] [options]
    elif args.subcommand_name == 'profile':
        sys.stdout.write('Start mapping the {} file...\n'.format(args.format.upper()))
        counting = ReadCount.read_count_on_DHS(args)
        positions = sorted(counting.keys(),key=counting.get,reverse=True)

