from labtools.adtools.finder import pull_AD, pull_barcode
from labtools.adtools.seqlib import read_fastq, read_fastq_big
import pandas as pd

def seq_counter(fastq, design_to_use = None, barcoded = False, only_bcs = False, **kwargs):
    """Counts occurences of ADs or AD-barcode pairs in a fastq file.
    
    Parameters
    ----------
    fastq : str 
        Path to fastq or fastq.gz file.
    design_to_use : str, default None
        Path to csv file containing ArrayDNA column.
    barcoded : bool, default False
        Whether to count ADs with different barcodes separately.
    only_bcs : default False
        True, False or the barcode map to use. If True, no map is used.
    **kwargs : dict
        Add additional arguments to pass to pull_AD or pull_barcode.
    
    Returns
    ----------
    counts : pandas.core.series.Series
        Pandas series where indices are AD or AD/barcode sequences and values are counts.
    
    Examples
    ----------
    >>> seq_counter("../exampledata/mini.fastq")
    GGTTCTTCTAAATTGAGATGTGATAATAATGCTGCTGCTCATGTTAAATTGGATTCATTTCCAGCTGGTGTTAGATTTGATACATCTGATGAAGAATTGTTGGAACATTTGGCTGCTAAA    1
    GAAGAATTGTTTTTACATTTGTCTGCTAAGATTGGTAGATCTTCTAGGAAACCACATCCATTCTTGGATGAATTTATTCATACTTTGGTTGAAGAAGATGGTATTTGTAGAACTCATCCA    3
    dtype: int64
    """
    seqCounts = {}
    # LT: Extract loss_table from kwargs and make design_loss_table to hold temporary loss_table calculations
    loss_table = kwargs.pop('loss_table', None) 
    design_loss_table = {'total': 0, 'filtered': 0}

    # merge lists of fastq files pertaining to the same sample
    if type(fastq) == list:
        # compute by only counting barcodes
        if only_bcs != False and design_to_use == None:
            for file in fastq:
                for line in read_fastq_big(file, **kwargs):
                    if loss_table is not None:
                        bc,loss_reason = pull_barcode(line[1], loss_table=loss_table, **kwargs)
                        if loss_reason is not None:
                            loss_table[loss_reason] += 1
                    else:
                        bc,loss_reason = pull_barcode(line[1], **kwargs)
                    
                    if bc not in seqCounts and bc != None: 
                        seqCounts[bc] = 1
                        # LT: increment total read counts for design_loss_table
                        design_loss_table['total'] += 1
                    elif bc != None: 
                        seqCounts[bc] += 1
                        # LT: increment total read counts for design_loss_table
                        design_loss_table['total'] += 1
                counts = pd.Series(seqCounts)
        # deny barcode search if a design file is provided
        elif only_bcs != False and design_to_use != None:
            raise TypeError("Using a design file is not compatible with only barcodes.")
        # compute by counting ADs or ADs + barcodes
        else:        
            for file in fastq:
                for line in read_fastq_big(file, **kwargs):
                    if loss_table is not None:
                        AD,bc,loss_reason = pull_AD(line[1], barcoded, loss_table=loss_table, **kwargs)
                        if loss_reason is not None:
                            loss_table[loss_reason] += 1
                    else:
                        AD,bc,loss_reason = pull_AD(line[1], barcoded, **kwargs)
                    
                    if barcoded and AD != None:
                        pair = (AD, bc)
                    
                        if pair not in seqCounts and pair[0] != None: seqCounts[pair] = 1
                        elif pair[0] != None: seqCounts[pair] += 1
                        # LT: increment total read counts for design_loss_table
                        design_loss_table['total'] += 1

                    elif AD != None and AD not in seqCounts: 
                        seqCounts[AD] = 1
                        # LT: increment total read counts for design_loss_table
                        design_loss_table['total'] += 1
                    elif AD != None: 
                        seqCounts[AD] += 1
                        # LT: increment total read counts for design_loss_table
                        design_loss_table['total'] += 1
            counts = pd.Series(seqCounts)
    # fastq files are not provided in lists (one file per sample)
    else:
        # compute by only counting barcodes
        if only_bcs != False and design_to_use == None:
            for line in read_fastq_big(fastq, **kwargs):
                if loss_table is not None:
                    bc,loss_reason = pull_barcode(line[1], loss_table=loss_table, **kwargs)
                    if loss_reason is not None:
                        loss_table[loss_reason] += 1
                else:
                    bc,loss_reason = pull_barcode(line[1], **kwargs)
                
                if bc not in seqCounts and bc != None: 
                    seqCounts[bc] = 1
                    # LT: increment total read counts for design_loss_table
                    design_loss_table['total'] += 1
                elif bc != None: 
                    seqCounts[bc] += 1
                    # LT: increment total read counts for design_loss_table
                    design_loss_table['total'] += 1
            counts = pd.Series(seqCounts)
        # deny barcode search if a design file is provided
        elif only_bcs != False and design_to_use != None:
            raise TypeError("Using a design file is not compatible with only barcodes.")
        # compute by counting ADs or ADs + barcodes
        else:
            for line in read_fastq_big(fastq, **kwargs):
                if loss_table is not None:
                    AD,bc,loss_reason = pull_AD(line[1], barcoded, loss_table=loss_table, **kwargs)
                    if loss_reason is not None:
                        loss_table[loss_reason] += 1
                else:
                    AD,bc,loss_reason = pull_AD(line[1], barcoded, **kwargs)
                
                if barcoded and AD != None:
                    pair = (AD, bc)
                
                    if pair not in seqCounts and pair[0] != None: seqCounts[pair] = 1
                    elif pair[0] != None: seqCounts[pair] += 1
                    # LT: increment total read counts for design_loss_table
                    design_loss_table['total'] += 1
                    
                elif AD != None and AD not in seqCounts: 
                    seqCounts[AD] = 1
                    # LT: increment total read counts for design_loss_table
                    design_loss_table['total'] += 1
                elif AD != None: 
                    seqCounts[AD] += 1
                    # LT: increment total read counts for design_loss_table
                    design_loss_table['total'] += 1
            counts = pd.Series(seqCounts)
    
    # remove non-perfect matches if required
    if design_to_use:
        design = pd.read_csv(design_to_use)
        if barcoded:
            counts = counts.where(counts.index.droplevel(1).isin(design.ArrayDNA)).dropna()
        else:
            counts = counts.where(counts.index.isin(design.ArrayDNA)).dropna()
    # LT: Calculate the number of reads filtered because of the design file
    design_loss_table['filtered'] = design_loss_table['total'] - counts.sum()
    return counts, design_loss_table

def create_map(ad_bcs, filter = False):
    """Converts output of seq_counter with AD,bc pairs to a dict map.
    
    If the barcode is found with two different ADs, it is not included in 
    the dictionary.

    Parameters
    ----------
    ad_bcs : pd.Series 
        output counts from seq_counter with barcoded = True.
    filter : int, default False
        Number of reads below which to ignore the barcode.

    Returns
    ----------
    bc_dict : dict
        Dictionary with barcodes as keys and 1 AD as value.
    """
    bc_dict = {}
    bad_bcs = 0
    for line in zip(ad_bcs.index, ad_bcs):
        ad = line[0][0]
        bc = line[0][1]
        count = line[1]
        if filter:
            if count < filter:
                pass
            elif count >= filter:
                if bc not in bc_dict:
                    bc_dict[bc] = ad
                elif bc in bc_dict and bc_dict[bc] == ad:
                    pass
                elif bc in bc_dict and bc_dict[bc] != ad:
                    del bc_dict[bc]
                    bad_bcs += 1
        else:
            if bc not in bc_dict:
                bc_dict[bc] = ad
            elif bc in bc_dict and bc_dict[bc] == ad:
                pass
            elif bc in bc_dict and bc_dict[bc] != ad:
                del bc_dict[bc]
                bad_bcs += 1

    return(bc_dict, bad_bcs)

def convert_bcs_from_map(bcs, bc_dict):
    """Takes bc only data and uses a barcode dictionary to return AD counts.
    
    If the barcode is found with two different ADs, it is not included in 
    the dictionary.

    Parameters
    ----------
    bcs : pd.Series 
        output counts from seq_counter with only_bcs = True.
    bc_dict : dict
        Dictionary with barcodes as keys and 1 AD as value from create_map().

    Returns
    ----------
    converted : pd.Series
        Pandas series where indices are AD sequences and values are counts.
    """
    ads = []
    for bc in bcs.index:
        ad = None
        if bc in bc_dict:
            ad = bc_dict[bc]
        ads.append(ad)
    ad_col = pd.Series(ads)
    x = pd.DataFrame(bcs).reset_index()
    x["AD"] = ads
    converted = x[[0, "AD"]].groupby("AD").sum()[0]
    return converted


def sort_normalizer(pair_counts, bin_counts, loss_table=None, thresh = 10):
    """Normalize by reads per sample, reads per tile and reads per bin.
    
    Parameters
    ----------
    pair_counts : list of pandas.core.series.Series 
        List of pandas series where indices are AD or AD/barcode sequences and values are counts.
    bin_counts : list
        List of number of cells per bin in the same order as the pair counts.
    loss_table: dictionary
        Dictionary that stores the value of reads removed during various filtering steps.
    thresh : int, default 10
        Number of reads above which to count the unique sequence.
    
    Returns
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the normalized read counts.
    numreads : pandas.DataFrame
        Total read counts for each unique sequence.
    reads : pandas.DataFrame
        Read counts per bin for each unique sequence.
    
    Examples
    ----------
    >>> sort_normalizer([count1, count2], [1000,1000])
    """
    df = pd.DataFrame(pair_counts)
    df.fillna(0, inplace=True)
    # LT: counting reads that will be removed by thresholding
    loss_table['thresh'] += df.loc[:, (df.sum() <= thresh)].sum().sum()
    # 10 is the read minimum, should make this user defined
    df = df.loc[:, (df.sum() > thresh)]
    df = df.transpose()
    numreads = df.sum(axis = 1)
    reads = df.copy(deep = True)
    #df = df_in.copy(deep=True)
    # row i column j
    for j in df:
        df[j] = (df[j]/df[j].sum())/bin_counts[j]
    for i, pair in enumerate(df.index):
        df.iloc[i] = df.iloc[i]/df.iloc[i].sum()
    
    return df, numreads, reads, loss_table

def calculate_activity(df_in, bin_values, min_max = False):
    """Calculate the activity of a normalized sort df.
    
    Parameters
    ----------
    df_in : pandas.DataFrame
        Dataframe output of sort_normalizer()
    bin_values : list
        List of mean or median fluorescence values per bin in the same order as the pair counts.
    min_max : bool, default False
        Whether to normalize the activity using min 0 max 1.
    
    Returns
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the activity values per sequence or sequence-barcode pair.
    """
    df = df_in.copy(deep=True)
    activities = df_in.dot(bin_values)
    if min_max:
        activities = minmax_scale(activities)
    
    df.loc[:,"Activity"] = activities
    return df

# def seq_counter_parallel(fastq, design_to_use = None, barcoded = False, only_bcs = False, **kwargs):
#     """Counts occurences of ADs or AD-barcode pairs in a fastq file in parallel. THIS IS BETA
    
#     Parameters
#     ----------
#     fastq : str 
#         Path to fastq or fastq.gz file.
#     design_to_use : str, default None
#         Path to csv file containing ArrayDNA column.
#     barcoded : bool, default False
#         Whether to count ADs with different barcodes separately.
#     only_bcs : default False
#         True, False or the barcode map to use. If True, no map is used.
    
#     Returns
#     ----------
#     counts : pandas.core.series.Series
#         Pandas series where indices are AD or AD/barcode sequences and values are counts.
    
#     Examples
#     ----------
#     >>> seq_counter("../exampledata/mini.fastq")
#     GGTTCTTCTAAATTGAGATGTGATAATAATGCTGCTGCTCATGTTAAATTGGATTCATTTCCAGCTGGTGTTAGATTTGATACATCTGATGAAGAATTGTTGGAACATTTGGCTGCTAAA    1
#     GAAGAATTGTTTTTACATTTGTCTGCTAAGATTGGTAGATCTTCTAGGAAACCACATCCATTCTTGGATGAATTTATTCATACTTTGGTTGAAGAAGATGGTATTTGTAGAACTCATCCA    3
#     dtype: int64
#     """
#     seqCounts = {}
    
#     def helper(file, design_to_use, barcoded, only_bcs, seqCounts, **kwargs):
#         for line in read_fastq_big(file, **kwargs):
#             AD,bc = pull_AD(line[1], barcoded, **kwargs)
                
#             if barcoded and AD != None:
#                 AD = (AD, bc)
#             if AD not in seqCounts and AD != None: seqCounts[AD] = 1
#             elif AD != None: seqCounts[AD] += 1
#         return seqCounts

#     # merge lists of fastq files pertaining to the same sample
#     if type(fastq) == list:
#         for file in fastq:
#             for line in read_fastq_big(file, **kwargs):
#                 seqCounts = helper(file, design_to_use, barcoded, only_bcs, seqCounts, **kwargs)
#         counts = pd.Series(seqCounts)

#     # remove non-perfect matches if required
#     if design_to_use:
#         design = pd.read_csv(design_to_use)
#         if barcoded:
#             counts = counts.where(counts.index.droplevel(1).isin(design.ArrayDNA)).dropna()
#         else:
#             counts = counts.where(counts.index.isin(design.ArrayDNA)).dropna()
#     return counts

def main():
    pass

if __name__ == '__main__':
    main
