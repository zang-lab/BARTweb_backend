import argparse,math,numpy,os,sys,tables
import time,re
import numpy as np

# read in list
def read_list(fname):
    fp = open(fname)
    l = []
    for line in fp.readlines():
        if not ( line[0] == '#' ):
            f = line.strip().split()
            if not f:
                continue
            else:
                sym = f[0]
                l  += [sym]
    fp.close()
    return l

# read in H3K27ac sample IDs
def read_sample_list(fname):
    fp = open(fname)
    sl,scores = [],[]
    for line in fp:
        line  = line.strip().split('\t')
        #print(line)
        if re.match('[0-9]',line[0]):
            sl += [line[0]]
            scores += [float(line[1])]
    fp.close()
    return sl,scores

# read in H3K27ac signal on union DHS sites for sample
def read_hdf5( h5file, sample_names ):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 

    #h5file = tables.open_file( file_name, driver="H5FD_CORE")    
    X = None

    for elem in sample_names:
        a = h5file.get_node("/", elem )
        m = a.read()
        if X is None:
            X = m
        else:
            X = numpy.vstack((X,m))

    X = X.transpose()
    return X

# rank numpy vector 
def rank(a):
    """ Rank from low to high """
    temp = a.argsort()
    rank = numpy.empty(len(a), int)
    rank[temp] = numpy.arange(len(a))
    return rank

def sqrttansform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount
    med = np.median( xt,0 )
    x = np.sqrt(xt) - np.sqrt(med)
    return x


def median_norm(X):
    col_median = numpy.median(X,0)
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - col_median[i]
    return X

def quantile_norm(X):
    """ input: numpy array data_points x features """
    Y = numpy.copy(X)
    Y.sort(0)
    row_mu = Y.mean(1)
    for i in range(Y.shape[1]):
        r = rank(X[:,i])
        X[r,i] = row_mu
    return X

def rowcenter(X):
    X = (X.transpose()-X.mean(1)).transpose()
    return X

def read_table(tffile,col=1,ctype=int):
    fp = open(tffile)
    l = []
    for line in fp.readlines():
        f = line.split()
        l += [ctype(f[col])]
    return l

def get_H3K27ac_rpkm(sample_name, UAC_H3K27ac, idx=1):
    RPKM_file = os.path.join(UAC_H3K27ac, sample_name + '_Strength.txt')
    if not os.path.exists(RPKM_file):
        sys.stderr.write('{} RPKM sample does not exist!'.format(sample_name))
        return None

    dhs = []
    with open(RPKM_file, 'r') as fopen:
        for line in fopen:
            dhs.append(float(line.strip().split('\t')[idx]))

    return dhs


def main(samplefile,genefile,casename,genome,UAC_info, UAC_H3K27ac,tffile=None):  # UAC_info to extract chrom, start, end, ID; UAC_H3K27ac is a directory

    fpo = open(casename+'_enhancer_prediction_lasso.txt','w')

    sample_names,sample_weights = read_sample_list(samplefile)
    sample_names = [i for i in sample_names if re.match('[0-9]',i)] # incase 'start/end/chr' in sample names
    # sample_names = [i for i in sample_names if re.match('GSM',i)] # incase 'start/end/chr' in sample names
    sample_weights = np.array(sample_weights)
    
    DHS_sample_names = [ elem+'_Strength' for elem in sample_names ]
    print(DHS_sample_names)
    #print(DHS_sample_names)
    # udhs_h5file = tables.open_file( UDHS_H3K27ac_HDF5s, driver="H5FD_CORE")    
    udhs_h5file = tables.open_file( UAC_H3K27ac, driver="H5FD_CORE")    
    DHS = None
    for dhs_samplename in DHS_sample_names:
        dhs = read_hdf5( udhs_h5file, [dhs_samplename] )
        # dhs = get_H3K27ac_rpkm(dhs_samplename, UAC_H3K27ac)
        dhs = numpy.array(dhs)
        dhs = dhs.transpose()
        if DHS is None:
            DHS = dhs
        else:
            DHS =  numpy.vstack((DHS,dhs))

    udhs_h5file.close()
    DHS = DHS.transpose()
    # DHS = sqrttansform(DHS)

    DHS = numpy.sqrt(DHS)
    DHS = median_norm(DHS)
    DHS = rowcenter(DHS)
        
    T = np.dot(DHS,sample_weights)#;print(T,T.shape)

    # #h5file = UDHS_H3K27ac_HDF5s[0]
    # udhs_h5file = tables.open_file( UDHS_H3K27ac_HDF5s, driver="H5FD_CORE")    
    # chrom = read_hdf5( udhs_h5file, ["chrom"] )
    # start = read_hdf5( udhs_h5file, ["start"] )
    # end = read_hdf5( udhs_h5file, ["end"] )
    # # ID = read_hdf5( udhs_h5file, ["end"] )
    # ID = read_hdf5( udhs_h5file, ["ID"] )
    # udhs_h5file.close()

    chrom = []
    start = []
    end = []
    with open(UAC_info, 'r') as fopen:
        for line in fopen:
            res = line.strip().split('\t')
            chrom.append(res[0])
            start.append(int(res[1]))
            end.append(int(res[2]))

    out_res = []
    for i,elem in enumerate(T):
        out_res.append((chrom[i],start[i],end[i],str(i+1),elem))
    sorted_out_res = sorted(out_res, key=lambda x: float(x[4]),reverse=True) 

    print('%s\t%s\t%s\t%s\t%s' % ("chromosom",'start','end','UDHSID',"Score"), file=fpo)
    for line in sorted_out_res:
        print('%s\t%d\t%d\t%s\t%3.6f' % (line[0],line[1],line[2],line[3],line[4]), file=fpo)
    fpo.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Predict enhancer elements from gene list and H3K27ac in union DNase-seq peaks.""")
    parser.add_argument( '-s', dest='samplefile', type=str, required=True, help='File that lists informative H3K27ac samples. One sample ID per line.' )
    parser.add_argument( '-i', dest='inputfile', type=str, required=True, help='Input gene list file. Gene symbols of refseq IDs.' )
    parser.add_argument( '-r', dest='refgeneAnn', type=str, required=True, help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>' )
    parser.add_argument( '--rp', dest='rphdf5', type=str, required=True, help='Path for the the hdf5 format rp files.' )
    parser.add_argument( '--k27ac', dest='k27achdf5', type=str, required=True, help='Path for the hdf5 format H3K27ac reads counts in UDHS regions.' )
    parser.add_argument( '-n', dest='name', type=str, required=True, help='Name of study, for output file naming.' )
    parser.add_argument( '-g', dest='genome', default='hg38', choices=['mm9','mm10','hg19','hg38'], required=False, help='genome' )
    parser.add_argument( '-t', dest='tffile', default=None, required=False, help='Indicators of TF binding at each UDHS sites 0 or 1. For performance evaluation.' )
    
    args = parser.parse_args()
    # main( args.samplefile, args.inputfile, args.name, args.genome, args.rphdf5, args.k27achdf5, args.refgeneAnn, tffile=args.tffile )
    # main( args.samplefile, args.inputfile, args.name, args.genome, args.k27achdf5, tffile=args.tffile )
    UAC_info = '/nv/vol190/zanglab/wm9tr/data/marge_data/RPKM/hg38_splitted_udhs_equally.bed'
    main( args.samplefile, args.inputfile, args.name, args.genome, UAC_info, args.k27achdf5, tffile=args.tffile )
