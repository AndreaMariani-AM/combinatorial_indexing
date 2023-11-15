# import os
# import sys

# from pathlib import Path
# from textwrap import wrap
import argparse
import pandas as pd
# Import custom functions
import utils 

######################
## ARGUMENT PARSING ##
######################
parser = argparse.ArgumentParser(description='Demultiplexing samples based on the first barcode that has been ligated.')
parser.add_argument('--alevin_fry_dir', 
                    help='This is the main alevin-fry direcotry.', required=True)

parser.add_argument('--barcode_mapping_dir', 
                    help='Directory that contains the barcode to samples mapping file.', required=True)

parser.add_argument('--mapping_file', 
                    help='The barcodes to samples mapping file.', required=True)

parser.add_argument('--filter_barcodes', 
                    help="Whether or not to filter barcodes that cannot be associated with a sample. "
                        "Default is not. If you want to filter barcodes, set this flag.", 
                    action="store_true", required=False)

# Reading in inputs
args = parser.parse_args()

alevin_fry_dir      = args.alevin_fry_dir
barcode_mapping_dir = args.barcode_mapping_dir
mapping_file        = args.mapping_file
filter_barcodes     = args.filter_barcodes

# testing the funciton
# alevin_fry_dir="/hpcnfs/scratch/DP/amariani/Dcotugno/SPLiT-seq/data/alevin_fry/counts/alevin"
# barcode_mapping_dir="/hpcnfs/scratch/DP/amariani/Dcotugno/SPLiT-seq/data/alevin_fry/bc_ex_mapping"
# mapping_file="demultiplexing_barcodes.txt"

utils.demultiplex_samples(alevin_fry_dir=alevin_fry_dir,
                          barcode_mapping_dir=barcode_mapping_dir,
                          mapping_file=mapping_file, filt_other=filter_barcodes)
