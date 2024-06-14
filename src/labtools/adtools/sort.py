from labtools.adtools.counter import *
from labtools.adtools.seqlib import get_numreads, read_bc_dict
from labtools.adtools.counter import seq_counter
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

class Sort():
    """A sort object.
    
    Attributes
    ----------
    data_files : list
        List of str paths to fastq files in order of bins.
    bin_counts : list
        List of number of cells per bin in order of data files.
    bin_values : list
        List of mean or median fluorescence values per bin in order of data files.
    design_file : str or None
        Path to design csv file containing "ArrayDNA" header with DNA sequences to search
    bc_dict : str or False, default False
        Path to bc_dict csv file (see seqlib.write_bc_dict)
    """

    def __init__(self, data_files: list, bin_counts: list, bin_values: list, 
                 design_file: str=None, bc_dict: str=False):
        self.data_files = data_files
        self.bin_counts = bin_counts
        self.bin_values = bin_values
        self.design_file = design_file
        self.bc_dict = bc_dict

    def process(self, csv=False, **kwargs):
        """Calculate the activity for each tile.
        
        Parameters
        ----------
        csv : bool, default False
            Use a csv input for predetermined AD counts (for example, bowtie output).
        **kwargs : dict, optional
            Extra arguments to pull_AD (for example changing the anchor sequences).
            
        Returns
        ----------
        processed_sort : pandas.DataFrame
        numreads : 
        reads :
       
        Examples
        ----------
        >>>Sort.process()
        """

        sort_list = []
        # LT: Initialize loss_table
        # Note: bc_flanks and bc_length are not actually removed. 
        # Those reads become bc_is_na if they pass the design file and thresholding.
        loss_table = {'ad_preceder': 0, 'design_file': 0, 'thresh': 0, # This line shows truly removed reads. Sum of this line + sum of included reads should = total reads.
                      'bc_is_na': 0, 'bc_flanks': 0, 'bc_length': 0}
        
        # super weird and needs to be double checked
        # converts output from bowtie bam to csv into same format
        # as the seq_counter() fxn output for downstream processing
        if csv == True: 
            for sample in self.data_files:
                parsed_sample = pd.read_csv(sample, index_col=1, header = None).squeeze('columns')
                sort_list.append(parsed_sample)
        else: 
            for sample in self.data_files:
                parsed_sample, design_loss_table = seq_counter(sample, loss_table=loss_table, design_to_use = self.design_file, 
                                            only_bcs = self.bc_dict, **kwargs)
                # LT: Assign calculations from temporary design_loss_table to loss_table for each file.
                loss_table['design_file'] += design_loss_table['filtered']
                if self.bc_dict != False and self.bc_dict != True:
                    bd = read_bc_dict(self.bc_dict)
                    parsed_sample = convert_bcs_from_map(parsed_sample, bd)
                sort_list.append(parsed_sample)
        
        normed_sort, numreads, reads, loss_table = sort_normalizer(sort_list, self.bin_counts, 
                                                                   thresh = kwargs.get("thresh", 10), 
                                                                  loss_table=loss_table)
        # LT: Count reads where barcode is None and record in the loss_table
        normnobcs = reads[reads.index.get_level_values(1).isna()]
        loss_table['bc_is_na'] += normnobcs.sum(axis=1).sum()

        processed_sort = calculate_activity(normed_sort, self.bin_values)

        return processed_sort, numreads, reads, loss_table



    def downsample(self, subset_size: int) -> pd.DataFrame:
        """Perform downsampling on raw reads, then error compared to using all reads.

        Parameters
        ----------
        subset_size : int
            Size in # of reads to subset each file by. 
            
        Returns
        ----------
        errors : pandas.DataFrame

        Examples
        ----------
        >>>Sort.downsample(100000)
        """
        
        sort_list = []
        for sample in self.data_files:
            numreads = get_numreads(sample)
            subsample_list = []
            
            for num in range(0, numreads, subset_size):
                if num == 0:
                    num = None
                subsample = seq_counter(sample, design_to_use = self.design_file, subset = num)
                subsample_list.append(subsample)
            sort_list.append(subsample_list)
        
        all_subsamples = []

        min_reads = min([len(sort) for sort in sort_list])
        for subsample in range(0, min_reads):
            group = [sort_list[i][subsample] for i in range(0, len(sort_list))]
            df, _, _ = sort_normalizer(group, self.bin_counts)
            final = calculate_activity(df, self.bin_values)
            all_subsamples.append(final)
        

        for i in range(1, len(all_subsamples)):
            all_subsamples[0][f"{i}"] = all_subsamples[i].Activity

        downsampling_df = all_subsamples[0].iloc[:,len(sort_list):]

        RMSE_list = []
        for j in downsampling_df:
            if j == "Activity":
                RMSE = downsampling_df.Activity
            else:
                RMSE = ((downsampling_df.Activity - downsampling_df[j]) ** 2 ) ** 0.5
            RMSE.rename(j, inplace = True)
            RMSE_list.append(RMSE)

        errors = pd.DataFrame(RMSE_list)
        plot_error = errors.transpose().iloc[:,1:]

        fig, axes = plt.subplots(1,1, figsize=(18, 10), sharex = True)
        axes.set_title("Subsampling RMSE", size=14)
        #axes.set_ylim(0,75)
        axes.set_ylabel("RMSE")
        axes.set_xlabel(f"Number of {subset_size} reads")
        sns.boxplot(data = plot_error)

        return errors.transpose()
