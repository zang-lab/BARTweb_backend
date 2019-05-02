import zlib
import argparse,math,os,sys,tables
import numpy as np
from operator import itemgetter, attrgetter, methodcaller
from sklearn import linear_model, decomposition, datasets
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn import cross_validation
from sklearn import metrics
#import regions
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


UNCHANGED,DOWN,UP,TARGET = 0,1,2,3
statdict = { UNCHANGED:'.', DOWN:'DOWN', UP:'UP', TARGET:'TARGET' }

class dataset(object):
    def __init__(self):
        self.index = None
        self.info  = None # include chrom and TSS and gene symbol
        self.x     = None
        self.y     = None
        self.bwdir = None 
        self.bwfiles = None
        self.rpfiles = None 

def gene_sym(symfile):
    """
    One representative of each gene symbol is chosen. 
    """
    # return {"refseq":["symbol","chr"]}
    fp = open(symfile)
    symdict = {} # {"gene_symbol": ["refseq_ID", "chr"]}  the first refseq ID in the TSS file
    for line in fp:
        sline = line.strip().split('\t')
        symbol = sline[3]
        chr_info = sline[0]
        symdict[symbol] = chr_info
        
    # for line in fp:
    #     f = line.strip().split('\t')
    #     #print f
    #     g = f[3]
    #     IDs = g.split(':')
    #     #print IDs
    #     if IDs[1] not in symdict:
    #         symdict[IDs[1]] = [IDs[0], f[0]]
 
    # rsymdict = {}
    # for elem in symdict:
    #     rsymdict[symdict[elem][0]] = [elem,symdict[elem][1]]
    fp.close()

    return symdict

def scalartransform(x):
    pcount = 1
    x = np.log2(x+pcount)
    return x

def transform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount 
    med = np.median( xt,0 )
    x = np.log2(xt) - np.log2(med)
    return x

def sqrttansform(x):
    xt = np.array(x)
    pcount = 1
    xt += pcount
    med = np.median( xt,0 )
    x = np.sqrt(xt) - np.sqrt(med)
    return x
#====================================================================

def read_hdf5_ori( h5file, sample_names ):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 
    X = None
    for elem in sample_names:
        a = h5file.get_node("/", elem )
        m = a.read()
        if X is None:
            X = m
        else:
            X = numpy.vstack((X,m))

    X1 = X.transpose()
    return X

def read_hdf5( h5file, sample_name):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 
    a = h5file.get_node("/", sample_name )
    m = a.read()
    return m


def getSampleNames_hdf5(h5file):
    samplenames = []
    for array in h5file.walk_nodes("/","Array"):
        if array.name not in samplenames:
            samplenames.append(array.name)
        else:
            continue
    #samplenames = samplenames[340:]
    return samplenames

def readregpotfiles(sym,genome,samplenames,h5file):

    # make list of regulatory potential file
    # read in regulatory potential from files in directory 

    index = None
    x = None
    nsamples = len(samplenames)

    # refseqID = read_hdf5( h5file, 'refseqID' )
    # print(refseqID)
    symID = read_hdf5(h5file, 'symbol')
    print(symID)
    chrom = read_hdf5(h5file, 'chr')
    #print(chrom)
    start = read_hdf5(h5file, 'start')
    #print(start)
    for k,name in enumerate(samplenames):
        if index == None:
            print(name)
            index = {}
            info = {}
            i = 0
            # for j,geneid in enumerate(refseqID):
            for j,geneid in enumerate(symID):
                geneid = geneid.decode("utf-8")  # gene symbol
                if geneid in sym:
                    # symid = sym[geneid][0]
                    # if symid not in index:
                    if geneid not in index:
                        # index[symid] = i
                        # info[symid] = [chrom[j].decode("utf-8"),start[j]]
                        index[geneid] = i
                        info[geneid] = [chrom[j].decode("utf-8"),start[j]]
                        i += 1
            ngenes = len(index)
            x = np.zeros((ngenes,nsamples))
            print(np.shape(x))
            
        RP = read_hdf5( h5file, name )
        # for i,geneid in enumerate(refseqID):
        for i,geneid in enumerate(symID):
            geneid = geneid.decode("utf-8")
            if geneid in sym:
                symid = geneid
                # symid = sym[geneid][0]
                # symid = sym[geneid][0]
                rp = RP[i]
                try:
                    x[index[symid],k] = rp  ### float num was ignored here, e.g., 'chr', 'refseqID', 'start', 'symbol'
                except:
                    pass
    
    z         = dataset()  
    z.rpfiles = samplenames 
    z.x = x # x.shape = ngenes,nsamples  # x is RP not relative RP, change to relative RP
    print(np.median(x, axis=1))

    z.index = index # {symbol:'start position'}
    z.info  = info # {'symbol':[chr,start]}
    return z


def read_genelistOnly(sym, fname, index, gname2, sepby_chrom=True):
    
    status = np.zeros( len(index) )
    genenames = np.ndarray(shape=(len(index)),dtype=object)
    print(list(index.keys())[0:20])

    train_chroms = ['chr1','chr3','chr5','chr7','chr9','chr11','chr13','chr15','chr17','chr19','chr21']
    test_chroms = ['chr2','chr4','chr6','chr8','chr10','chr12','chr14','chr16','chr18','chr20','chr22']
    train_index = []
    test_index = []

    fp = open(fname).readlines()
    genes = [g.strip() for g in fp]

    allgenes = list(sym.keys())
    print(allgenes[0:20])

    for ag in allgenes:
        if gname2:
            try:
                # added by @marvinquiet
                # i = index[sym[ag][0]]
                # if sym[ag][1] in train_chroms:
                #     train_index.append(i)
                # elif sym[ag][1] in test_chroms:
                #     test_index.append(i)
                i = index[ag]
                if sym[ag] in train_chroms:
                    train_index.append(i)
                elif sym[ag] in test_chroms:
                    test_index.append(i)
                else:
                    pass
                #print i
                # if sym[ag][0] in genes:
                #     #print sym[ag][0]
                #     status[i] = TARGET
                # else:
                #     status[i] = UNCHANGED
                if ag in genes:
                    #print sym[ag][0]
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
                genenames[i] = ag
            except:
                continue
        else:
            try:
                i = index[sym[ag][0]]
                if sym[ag][1] in train_chroms:
                    train_index.append(i)
                elif sym[ag][1] in test_chroms:
                    test_index.append(i)
                else:
                    pass
                if ag in genes:
                    status[i] = TARGET
                else:
                    status[i] = UNCHANGED
                genenames[i] = ag
            except:
                continue

    print('file: %s\ttarget: %d\tunchanged: %d\n' % ( fname, sum( status == TARGET ), sum( status == UNCHANGED ) ))
    print(genenames[0:20])
    return (genenames, status,train_index,test_index)  

def dataset_annotation(annotationf):
    #get the cell annotation for each datasetID
    # change original description to 
    inf = open(annotationf,'rU')
    ann = {}
    for line in inf:
        if line.startswith('datasetID'):
            pass
        else:
            # line = line.strip().split('\t')
            # ID = line[0] # dataset id -> GSM id
            # info = [line[4],line[5],line[7]] # CellLineName, Tissue/Organ, DetailedTissue
            # try:
            #     ann[ID] = info
            # except:
            #     ann[ID] = 'NA'
            line = line.strip().split(',')
            ID = line[1] # dataset id -> GSM id
            info = [line[2],line[3],line[4]] # CellLineName, Tissue/Organ, DetailedTissue
            try:
                ann[ID] = info
            except:
                ann[ID] = 'NA'
    return ann


def regress(x, y, train_index, test_index, samplenames, name, change, maxsamples,sepby_chrom, ann, genenames):
    print(x)
    fout = open(name+'_auc_record.txt','w')
    fout2 = open(name + '_predited_scores.txt','w')
    print(x.shape)
    print(y.shape)

    cross_validate_flag = True
    record   = []
    cvrecord = []

    col_list   = samplenames
    maxsamples = min( len(col_list), maxsamples )

    # forward step-wise regression - include the samples that produce the best cross validation results
    for j in range(maxsamples):
        print ('first j = %s'%str(j))
        if cross_validate_flag:
            cvscore = []

            for i,elem in enumerate(col_list):
                if i not in record:
                    #print ('second i = %s'%str(i))
                    trial = record + [i]
                    rand_idx = list(range(x.shape[0]))
                    if not sepby_chrom:
                        # select random subset for trial and testing
                        random.shuffle( rand_idx )
                        xt = x[:,trial]
                        xt = xt[rand_idx,:]
                        yt = y[rand_idx] 
                        xt.reshape( (x.shape[0],len(trial))) #???
                        X_train, X_test, y_train, y_test = cross_validation.train_test_split( xt, yt )
                        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
                        #LR_l1.fit(X_train,y_train)
                        #cvscore += [ np.mean( cross_validation.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ) ) ]
                        cvscore += [ np.mean( cross_validation.cross_val_score( LR_l1, xt, yt, scoring='roc_auc', cv=5 ) ) ]
                    else:
                        xt = x[:,trial]
                        yt = y
                        X_train = xt[train_index,:]
                        y_train = yt[train_index]
                        
                        X_train.reshape((X_train.shape[0],len(trial)))
                        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
                        LR_l1.fit(X_train,y_train)
                        cvscore += [ np.mean( cross_validation.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ) ) ]
                else:
                    cvscore += [0]
                    continue

            #record the auc values when do the first regression       
            if j == 0:
                for m in range(len(col_list)):
                    fout.write(col_list[m] + '\t' + str(cvscore[m]) + '\n')
            fout.close()

            k = np.argmax(cvscore) 
            cvrecord += [max(cvscore)]
            record += [k]
            print (j,'record',record)

    if sepby_chrom:
        xt = x[:,record]
        X_test = xt[test_index,:] 
        y_test = y[test_index]
        LR_l1.fit(X_test,y_test)
        y_test_hat = LR_l1.predict_log_proba(X_test)
        yhat = LR_l1.predict_log_proba(xt)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_test_hat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        print("==============")
    else:
        xt = x[:,record]
        xt.reshape( (xt.shape[0],len(record)) )
        LR_l1.fit(xt,y)
        yhat = LR_l1.predict_log_proba(xt)
        fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        
        for l in range(len(y)):
            fout2.write(genenames[l] + '\t' + str(y[l]) + '\t'  + str(yhat[l,1]) + '\n')
        fout2.close()

        plt.plot(fpr,tpr)
        plt.xlabel('true positive rate')
        plt.ylabel('false positive rate')
        plt.title('ROC curve for gene prediction')
        plt.savefig('%s_roc.png'%name)
        plt.show()
        print("==============")
    outf = open(name + '_' + change + '_regressionInfo.txt','w')
    for k,elem in enumerate(record):
        dataID = col_list[elem].split('_')[0]
        if dataID in list(ann.keys()):
            annInfo = ann[dataID]
        else:
            annInfo = ['NA','NA','NA']

        print(dataID, cvrecord[k], LR_l1.coef_[0,k], annInfo)
        outf.write(dataID + '\t' + str(cvrecord[k]) + '\t' + str(LR_l1.coef_[0,k]) + '\t' + '\t'.join(annInfo) + '\n')
    outf.write('AUC = %s'%(str(round(auc,3))))
    outf.close()   

    return LR_l1.coef_, [ col_list[i]  for i in record ], yhat, record

def lasso_test(x,y):
    
    # given x,y, return auc score
    LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
    LR_l1.fit(x,y)
    #np.mean( cross_validation.cross_val_score( LR_l1, X_train, y_train, scoring='roc_auc', cv=5 ))
    yhat = LR_l1.predict_log_proba(x)
    fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
    auc = metrics.auc(fpr,tpr)
    selected_features = len([i for i in LR_l1.coef_[0] if i !=0])
    return auc,selected_features
    
def lasso_test_best_alpha(x,y,prename):
    # given x,y, return alpha used for adaptive lasso
    alphas = [i for i in np.logspace(-2,1,10)]
    alpha_cvs = []
    plt.figure()
    for alpha in alphas:
        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01,fit_intercept=True,C=alpha);print(alpha)
        cvs_scores = cross_validation.cross_val_score( LR_l1, x, y, scoring='roc_auc', cv=5 )
        alpha_cvs.append(cvs_scores)
        
        LR_l1.fit(x,y)
        yhat = LR_l1.predict_log_proba(x)
        fpr, tpr, thresholds = metrics.roc_curve(y, yhat[:,1], pos_label=1)
        auc = metrics.auc(fpr,tpr)
        selected_features = len([i for i in estimator.coef_[0] if i !=0])
        print(alpha,np.mean(cvs_scores),auc,selected_features)
        # plot the auc figs
        y_mean = np.mean(cvs_scores)
        y_err  = np.std(cvs_scores)
        plt.errorbar(alpha,y_mean,y_err,color='r',ecolor='grey',fmt='o',capsize=4)
    plt.ylim([0,1])
    plt.xscale('log')
    plt.savefig(prename+'_alpha_auc.png',bbox_inches='tight',pad_inches=0.1,transparent=True)
    plt.close()
    #alpha_cvs_mean = [i.mean() for i in alpha_cvs]
    #best_alpha = alphas[alpha_cvs_mean.index(max(alpha_cvs_mean))]
    return alphas,alpha_cvs


def best_alpha(x,y):
    # given x,y, return alpha used for adaptive lasso
    alphas = [i for i in np.logspace(-2,1,10)]
    #alphas = [0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2]
    alpha_cvs = []
    
    for alpha in alphas:
        LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01,fit_intercept=True,C=alpha)
        cvs_scores = cross_validation.cross_val_score( LR_l1, x, y, scoring='roc_auc', cv=5 )
        alpha_cvs.append(cvs_scores) 
        print('  =best-alpha= ',alpha, '==mean-cvs==',np.mean(cvs_scores))
    alpha_cvs_mean = [i.mean() for i in alpha_cvs]
    best_alpha = alphas[alpha_cvs_mean.index(max(alpha_cvs_mean))]
   
    return best_alpha,max(alpha_cvs_mean)


def adaptive_lasso(x,y,samplefiles,name,maxsamples,ann,genenames):
    # test of adaptive lasso
 
    g = lambda w:np.sqrt(np.abs(w))
    gprime = lambda w: 1.0/(2.*np.sqrt(np.abs(w))+np.finfo(float).eps)
    n_samples,n_features = x.shape
    
    n_lasso_iterations = 10
    weights = np.ones(n_features)
    selected_features = n_features
    for k in range(n_lasso_iterations):
        if selected_features >maxsamples: 
            alpha=0.02
        else:
            alpha=0.2
        X_w = x / weights
        #alpha,best_cvs = best_alpha(X_w,y) # if you need to select best alpha for each step later
        #alpha = 0.1
        estimator = linear_model.LogisticRegression(penalty='l1', tol=0.01,fit_intercept=True,C=alpha)  # TODO: set fixed seed
        estimator.fit(X_w,y)
        coef_ = estimator.coef_/weights#;print(coef_)
        weights = gprime(coef_)
        selected_features = len([i for i in coef_[0] if i !=0])
        print(k,alpha,selected_features)

    rand_idx = list(range(x.shape[0]))
    random.shuffle( rand_idx )
    # xt = np.multiply(x,coef_);print(xt.shape)
    xt,yt = X_w[rand_idx,:], y[rand_idx]
    cvs_scores = cross_validation.cross_val_score(estimator ,xt,yt,  scoring='roc_auc', cv=5 )
    best_cvs = np.mean(cvs_scores)
    yhat = estimator.predict_log_proba(xt)
    fpr, tpr, thresholds = metrics.roc_curve(yt, yhat[:,1], pos_label=1)
    auc = metrics.auc(fpr,tpr)
    

#         print(k,'alpha',alpha)
#         print(k,'best_cvs',best_cvs)
#         print(k,'auc',auc)
#         print(k,'selected_features',selected_features)
# 
    outf = open('{}_adaptive_lasso_Info.txt'.format(name),'w')
    # for coef in [ i for i in coef_[0] if i!=0]:
    for coef in sorted([ i for i in coef_[0] if i!=0], key=abs, reverse=True):
        samplefile = samplefiles[list(coef_[0]).index(coef)]
        dataID = samplefile.split('_')[0]
        if dataID in list(ann.keys()):
            annInfo =  ann[dataID]
        else:
            annInfo = ['NA','NA','NA']

        outf.write('{}\t{}\t{}\n'.format(dataID, coef, '\t'.join(annInfo)))
    outf.write('AUC = {}\n'.format(auc))
    outf.write('best_cvs = {}\n'.format(best_cvs))
    outf.write('selected_features = {}\n'.format(selected_features))
    return auc,selected_features

def main( genome, rp_hdf5, gxfile, symfile, name, change, maxsamples, logtrans, sqrttrans, exptype,gname2, sepby_chrom, annotation):
    #symfile is the gene annotation file, change to gene symbol file
    basedir = os.path.dirname(name);print(basedir)
    os.makedirs(basedir,exist_ok=True)

    sym = gene_sym(symfile) # {"resfseq": {"gene_symbol", "chr"}} -> {"gene_symbol": "chr"}

    h5file = tables.open_file( rp_hdf5, driver="H5FD_CORE")
    #TODO: For testing the old GSM ids
    # samplenames  = getSampleNames_hdf5(h5file)

    old_samplenames = '/nv/vol190/zanglab/wm9tr/data/marge_data/RelativeRP/test_marge/marge_bart_pipeline/human_366_only_one_rep_samples_GSMs.txt' # 352 unique GSM ids
    # old_samplenames = '/nv/vol190/zanglab/wm9tr/data/marge_data/RelativeRP/test_marge/marge_bart_pipeline/human_366_samples_GSMs.txt' # 445 GSM ids
    samplenames = []
    with open(old_samplenames, 'r') as fopen:
        for line in fopen:
            samplenames.append(line.strip())

    z   = readregpotfiles(sym,genome,samplenames,h5file)
    h5file.close()
    
    if logtrans:
        #z.x = scalartransform(z.x)
        z.x = transform(z.x)
    if sqrttrans:
        z.x = sqrttansform(z.x)

    (genenames,z.y,train_index,test_index) = read_genelistOnly(sym, gxfile, z.index, gname2, sepby_chrom)

    if change == 'down':
        print('Do regrssion with DOWN genes')
        y = 1*( z.y == DOWN )
    elif change == 'up':
        print('Do regrssion with UP genes')
        y = 1*( z.y == UP )
    elif change == 'target':
        print('Do regrssion with TARGET genes')
        y = 1*( z.y == TARGET )
    else:
        print("Please choose the specfic direction, UP or DOWN or TARGET.")
        sys.exit()
    
    x = z.x[:,:-5] # remove the last few columns: refseq, start, chr, etc...
    ann = dataset_annotation( annotation )
    #coef,rprecord,yhat,record = regress( x, y, train_index, test_index, z.rpfiles, name, change, maxsamples, sepby_chrom, ann, genenames )


    #x = z.x#[:,:20]
    #z.rpfiles = z.rpfiles[:20]

    auc,selected_features = adaptive_lasso(x,y,z.rpfiles,name,maxsamples,ann,genenames)
    print(auc)
    print(selected_features)
    
    
if __name__ == "__main__":

    try:
        parser = argparse.ArgumentParser(description="""Regression of regulatory potential to gene expression changes.""")
        parser.add_argument( '-e','--expr', dest='expr', required = True, type = str, help = 'The related differential expression file')
        parser.add_argument( '--exptype', dest='exptype',required = True, choices=['Gene_Response','Gene_Only'], type = str, help = 'Gene_Response includes 2 columns, one is the geneID, and the other\
                                             is 1/0, 1 for target and 0 for un-target; Gene_Only includes 1 column, only the gene list of the targets. Only official gene symbol or \
                                             refseqID are allowd for the geneID.')
        parser.add_argument( '--gname2', dest='gname2', default=False, action='store_true', help = 'If this switch is on, gene or transcript IDs in files given through -e will be considered as official gene symbols, otherwise, it is RefSeqID, DEFAULT=FALSE')
        parser.add_argument( '-r','--historicalRP', dest='histRP', required = True, type = str, help = 'The file with hdf5 format which contain the H3K27ac RP information')
        parser.add_argument( '-n','--name', dest='name',required = True, type = str, help = 'The prefix of the output names')
        parser.add_argument( '-g','--genome', dest="genome", type=str, default='hg38', choices=['mm9','hg19','hg38','mm10'], required=False, help='genome')
        parser.add_argument( '--maxsamples', dest='maxsamples',  type=int, default=20, required=False, help='Maximum number of samples to include in regression model.' )
        parser.add_argument( '-l', '--logtrans', dest='logtrans',  default=False,  action='store_true', help='Use log transformation on regulatory potential (default is do nothing), this is complementary with -s, --sqrttrans.' )
        parser.add_argument( '-s', '--sqrttrans', dest='sqrttrans',  default=False,  action='store_true', help='Use sqrt transformation on regulatory potential (default is do nothing). this is complementary with -l, --logtrans.' )
        # parser.add_argument( '-m', dest='sym', type=str, required=True, help='refseqTSS is six columns: <chromosome name> <TSS> <TSS+1> <refseq:genesymbok> <score> <strand>')
        parser.add_argument( '-m', dest='sym', type=str, required=True, help='genesymbolTSS is six columns: <chromosome name> <TSS-1> <TSS+1> <genesymbol> <score> <strand>')
        parser.add_argument( '--sc', dest='sepby_chrom',  default=False,  action='store_true', help='Whether to seperate chroms to do cross validation, default = False' )
        parser.add_argument( '-a', dest='annotation', required=True,  type=str, help='The annotation file for each dataset' )

        args = parser.parse_args()

 
        main( args.genome, args.histRP, args.expr, args.sym, args.name, 'target', args.maxsamples, args.logtrans, args.sqrttrans, args.exptype, args.gname2, args.sepby_chrom, args.annotation )

    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)


