from labtools.adtools.finder import pull_AD
from labtools.adtools.seqlib import read_fastq, read_fastq_big
import pandas as pd

def seq_counter(fastq, design_to_use = None, barcoded = False, **kwargs):
    """Counts occurences of ADs or AD-barcode pairs in a fastq file.
    
    Parameters
    ----------
    fastq : str 
        Path to fastq or fastq.gz file.
    
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
    for line in read_fastq_big(fastq):
        AD,bc = pull_AD(line[1], barcoded, **kwargs)
        
        if barcoded and AD != None:
            AD = (AD, bc)
        if AD not in seqCounts and AD != None: seqCounts[AD] = 1
        elif AD != None: seqCounts[AD] += 1
    counts = pd.Series(seqCounts)
    
    if design_to_use:
        design = pd.read_csv(design_to_use)
        if barcoded:
            counts = counts.where(counts.index.droplevel(1).isin(design.ArrayDNA)).dropna()
        else:
            counts = counts.where(counts.index.isin(design.ArrayDNA)).dropna()
    return counts

def sort_normalizer(pair_counts, bin_counts):
    """Normalize by reads per sample, reads per tile and reads per bin.
    
    Parameters
    ----------
    pair_counts : list of pandas.core.series.Series 
        List of pandas series where indices are AD or AD/barcode sequences and values are counts.
    bin_counts : list
        List of number of cells per bin in the same order as the pair counts.
    
    Returns
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the normalized read counts.
    numreads : 
        unknown
    reads :
        unknown
    
    Examples
    ----------
    >>> sort_normalizer([count1, count2], [1000,1000])
    """
    df = pd.DataFrame(pair_counts)
    df.fillna(0, inplace=True)
    # 10 is the read minimum, should make this user defined
    df = df.loc[:, (df.sum() > 10)]
    df = df.transpose()
    numreads = df.sum(axis = 1)
    reads = df.copy(deep = True)
    #df = df_in.copy(deep=True)
    # row i column j
    for j in df:
        df[j] = (df[j]/df[j].sum())/bin_counts[j]
    for i, pair in enumerate(df.index):
        df.iloc[i] = df.iloc[i]/df.iloc[i].sum()
    
    return df, numreads, reads

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

def main():
    pass

if __name__ == '__main__':
    main