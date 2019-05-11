import os,sys,argparse
import tables
import time,re
import numpy as np

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

# rank np vector 
def rank(a):
    """ Rank from low to high """
    temp = a.argsort()
    rank = np.empty(len(a), int)
    rank[temp] = np.arange(len(a))
    return rank

def median_norm(X):
    col_median = np.median(X,0)
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - col_median[i]
    return X

def quantile_norm(X):
    """ input: numpy array data_points x features """
    Y = np.copy(X)
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

def main(samplefile, output_name, UAC_info, UAC_H3K27ac, tffile=None):
    '''generate cis-regulatory profile from selected H3K27ac samples

    Import arguments:
    - samplefile:     selected samples in *_adaptive_lasso_Info.txt
    - output_name:       output_dir/output_prefix
    - UAC_info:       UDHS or UAC bed file
    - UAC_H3K27ac:    RPKM on H3K27ac signals based on UDHS or UAC
    - tffile:         for background correction, default: None

    Output file:
    Enhancer profile based on UDHS or UAC.
    '''
    sample_names,sample_weights = read_sample_list(samplefile)
    sample_names = [i for i in sample_names if re.match('[0-9]',i)] # incase 'start/end/chr' in sample names
    sample_weights = np.array(sample_weights)
    
    DHS_sample_names = [ elem+'_Strength' for elem in sample_names ]
    print(DHS_sample_names)
    # read data from RPKM H3K27ac hdf5 file
    udhs_h5file = tables.open_file( UAC_H3K27ac, driver="H5FD_CORE")    
    chrom = read_hdf5(udhs_h5file, ["chrom"])
    start = read_hdf5(udhs_h5file, ["start"])
    end = read_hdf5(udhs_h5file, ["end"])
    ID = read_hdf5(udhs_h5file, ["ID"])
    DHS = None
    for dhs_samplename in DHS_sample_names:
        dhs = read_hdf5(udhs_h5file, [dhs_samplename])
        dhs = np.array(dhs)
        dhs = dhs.transpose()
        if DHS is None:
            DHS = dhs
        else:
            DHS =  np.vstack((DHS,dhs))
    udhs_h5file.close()
    DHS = DHS.transpose()

    DHS = np.sqrt(DHS)
    DHS = median_norm(DHS)
    DHS = rowcenter(DHS)
        
    T = np.dot(DHS,sample_weights) # coef * H3K27ac RPKM signals

    out_res = []
    for i,elem in enumerate(T):
        out_res.append((chrom[i].decode('utf-8'),start[i],end[i],str(i+1),elem))
    sorted_out_res = sorted(out_res, key=lambda x: float(x[4]),reverse=True)
    
    fpo = open(output_name+'_enhancer_prediction_lasso.txt','w')
    print('%s\t%s\t%s\t%s\t%s' % ("chromosom",'start','end','UDHSID',"Score"), file=fpo)
    for line in sorted_out_res:
        print('%s\t%d\t%d\t%s\t%3.2f' % (line[0],line[1],line[2],line[3],line[4]), file=fpo)
    fpo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Predict enhancer elements from gene list and H3K27ac in union DNase-seq peaks.""")
    parser.add_argument( '-s', dest='samplefile', type=str, required=True, help='File that lists informative H3K27ac samples. One sample ID per line.' )
    parser.add_argument( '--k27ac', dest='k27achdf5', type=str, required=True, help='Path for the hdf5 format H3K27ac reads counts in UDHS regions.' )
    parser.add_argument( '-n', dest='name', type=str, required=True, help='Name of study, for output file naming.' )
    parser.add_argument( '-t', dest='tffile', default=None, required=False, help='Indicators of TF binding at each UDHS sites 0 or 1. For performance evaluation.' )
    
    args = parser.parse_args()
    main( args.samplefile, args.name, args.genome, args.k27achdf5, tffile=args.tffile )
