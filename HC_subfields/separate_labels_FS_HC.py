#!/usr/bin/env python3
"""
Write separate volume for each label found in the input image
Usage
----
separate_labels.py <input label image>
separate_labels.py -h
Example
----
>>> separate_labels.py atlas.nii.gz 
"""

__version__ = '0.2.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Smooth one or more atlas labels')
    parser.add_argument('in_file', help="source atlas labels filename")
    
    args = parser.parse_args()
    
    in_file = args.in_file
    
    # Convert relative to absolute path
    in_file = os.path.abspath(in_file)
    
    # Load the source atlas image
    print('Opening %s' % in_file)
    in_nii = nib.load(in_file)
    
    # Load label image
    src_labels = in_nii.get_fdata()

    # Construct list of unique label values in image
    unique_labels = np.unique(src_labels).astype(int)

    print(unique_labels)

    # loop over each unique label value
    for label in unique_labels:
        
        if label > 0:
        
            # Create mask for current label value
            out_mask = (src_labels == label).astype(int)
            out_mask = np.asanyarray(out_mask)
            
            # Construct output filename. Use zero-padded indexing
            print(f'label:{label}')
            out_file = in_file.strip('.nii.gz') + f'_{label}.nii.gz'
            
            # Save smoothed labels image
            print('Saving label %d to %s' % (label, out_file))
            out_nii = nib.Nifti1Image(out_mask, in_nii.affine)
            out_nii.to_filename(out_file)
    
    print('Done')
    
    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()