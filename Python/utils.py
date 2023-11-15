import os
import sys

from pathlib import Path
from textwrap import wrap

import pandas as pd

# Demultiplex samples based on the first barcode ligated
def demultiplex_samples(alevin_fry_dir, barcode_mapping_dir, mapping_file, filt_other=False):
    """ 
    This function takes in the alevin-fry quantification directory, the barcode mapping directory and the name of the file that
    stores in the oligoDt to random primes mapping and the associated information of sample identity to be able to demultiplex
    different experiments. 
    Optionally, it can filter out barcodes that aren't identified by any of the sample names, (i.e. 1st barcodes isn't associated
    with any known well and it's the result of random nucleotides sequencing.
    It outputs in the barcode mapping directory, a .txt file that contains the sample info for each barcode that is later added as metadata
    to the sce/anndata object.
    """
    # BC
    # alevin-fry mtx directory, where matrices are stored
    bc = Path(os.path.join(alevin_fry_dir,'alevin','quants_mat_rows.txt'))
    # Check if file exists
    if bc.exists() is not True:
        raise FileNotFoundError("Barcode files doesn't exist or the input dir isn't an alevin-fry mtx directory. Quitting")

    # Barcode-samples mapping file
    # This is a csv file created by the experimentalists. Samples to demultiplex are in 'demultiplexing' column.
    # The order of barcodes is always the same so i can match by position
    # This is basically a copy of the oligoDt-randomPrimer file mapping
    mapping = Path(os.path.join(barcode_mapping_dir, mapping_file))
    # Check if the file exists
    if mapping.exists() is not True:
        raise FileNotFoundError("Demultiplexing file doesn't exist. Quitting")
    
    # Both files exist, read them in
    barcodes = pd.read_table(bc, header=None, names=['barcodes'])

    # this file should be read in with already sample info for demultiplexing on the right most column, needs to be standardized
    samples = pd.read_table(mapping, skiprows=1, names=['oligoDt','randomP', 'samples'])
    
    # Check if the sample files has the correct annotation and order
    if not (len(samples.columns)==3) and (samples.columns[2]=='samples'):
        raise Exception('Check demultiplexing files. Should have 3 columns')
        
    
    # split cells into single barcodes keeping only the last 8 bases (aka the first one ligated)
    barcodes['to_keep'] = barcodes['barcodes'].str[-8:]
    # Merge the two DFs
    barcodes_merged = barcodes.merge(samples, left_on="to_keep", right_on="oligoDt", how='left') 
    # Change NAs with other and drop unnecessary columns
    barcodes_merged.fillna('other', inplace=True)
    barcodes_merged.drop(['oligoDt', 'randomP'], axis=1, inplace=True)
    
    # Filter out other barcodes or not
    if filt_other is True:
        filt_df = barcodes_merged[~barcodes_merged["samples"].str.contains("other")]
        # write the file to disk to have it 
        filt_df.to_csv(os.path.join(barcode_mapping_dir, 'bc_sample_mapping.txt'), 
                        sep='\t',index=False)
    else:
        barcodes_merged.to_csv(os.path.join(barcode_mapping_dir, 'bc_sample_mapping.txt'), 
                        sep='\t',index=False)
