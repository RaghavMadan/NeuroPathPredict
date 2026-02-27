#!/usr/bin/env python3

"""
Script 05: Block Visualization & Assignment (NIfTI only)
Part of NeuroPathPredict Pipeline V1.1

This script performs KMeans clustering on coordinate data from a subject's
_uk_data.csv file for specified numbers of splits. It generates NIfTI visualizations
for each split using the input coordinates as voxel indices in the MNI152 2009b 0.5mm template.

MAPPING LOGIC:
- The input coordinates (X, Y, Z) are assumed to be in voxel space for the MNI152 2009b 0.5mm template.
- No affine transformation is applied; coordinates are rounded and used directly as voxel indices.
- The output NIfTI files are written with the standard MNI affine in the header for compatibility.

NOTE FOR FUTURE MODIFICATIONS:
- If your input coordinates are in MNI mm space (not voxel indices), you must apply the affine transformation:
    voxel_indices = np.round(np.dot(np.linalg.inv(affine), np.hstack((coords, np.ones((coords.shape[0], 1)))).T).T[:, :3]).astype(int)
- If your input coordinates are in a different voxel space or resolution, update the dims and affine accordingly.

Usage:
    python3 05_block_visualization.py <uk_data_file_path> <n_splits_str> <output_vis_dir> <subject_id>

Arguments:
    uk_data_file_path : Path to the subject's _uk_data.csv file (e.g., output/subject_id/universal_kriging/subject_id_uk_data.csv)
    n_splits_str      : Comma-separated string of n_splits values (e.g., "2,4,8,16")
    output_vis_dir    : Directory to save NIfTI visualizations (e.g., output/subject_id/universal_kriging/block_visualizations/)
    subject_id        : Subject ID for file naming

Outputs (in <output_vis_dir>):
    - {subject_id}_blocks_{n_splits}splits.nii.gz : NIfTI visualization for each n_split.
    - The input CSV is updated with new columns for each split (e.g., 2_block_splits, 4_block_splits, ...)
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import nibabel as nib
import os
import sys

def create_blocks(coords, n_splits):
    """Create blocks using KMeans clustering."""
    print(f"\nCreating {n_splits} blocks using KMeans clustering...")
    kmeans = KMeans(n_clusters=n_splits, random_state=42, n_init=10)
    labels = kmeans.fit_predict(coords)
    return labels + 1  # 1-based labels

def create_nifti_visualization(coords, labels, output_file):
    """Create NIfTI visualization of the blocks using voxel indices directly (input is already in voxel space)."""
    # MNI152 2009b dimensions (0.5mm isotropic)
    dims = np.array([394, 466, 378])
    volume = np.zeros(dims, dtype=np.int16)
    # Use coordinates directly as voxel indices
    voxel_indices = np.round(coords).astype(int)
    # Debug: print first 5 voxel indices
    print("First 5 voxel indices:")
    for i in range(min(5, voxel_indices.shape[0])):
        print(f"Voxel: {voxel_indices[i]}")
    valid_mask = np.all((voxel_indices >= 0) & (voxel_indices < dims), axis=1)
    valid_voxels = voxel_indices[valid_mask]
    valid_labels = labels[valid_mask]
    for idx, label in zip(valid_voxels, valid_labels):
        volume[idx[0], idx[1], idx[2]] = label
    # Use standard MNI affine for header
    affine = np.array([
        [0.5, 0, 0, -98],
        [0, 0.5, 0, -134],
        [0, 0, 0.5, -72],
        [0, 0, 0, 1]
    ])
    nifti_img = nib.Nifti1Image(volume, affine)
    nifti_img.header.set_zooms((0.5, 0.5, 0.5))
    nifti_img.header.set_sform(affine, code=1)
    nifti_img.header.set_qform(affine, code=1)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    nib.save(nifti_img, output_file)
    print(f"Saved block visualization to: {output_file}")

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 05_block_visualization.py <uk_data_file_path> <n_splits_str> <output_vis_dir> <subject_id>")
        print(__doc__)
        sys.exit(1)
    uk_data_file_path = sys.argv[1]
    n_splits_str = sys.argv[2]
    output_vis_dir = sys.argv[3]
    subject_id = sys.argv[4]
    if not os.path.exists(uk_data_file_path):
        print(f"Error: Input file not found: {uk_data_file_path}")
        sys.exit(1)
    df = pd.read_csv(uk_data_file_path)
    if not all(col in df.columns for col in ['X', 'Y', 'Z']):
        print("Error: CSV file must contain 'X', 'Y', 'Z' columns.")
        sys.exit(1)
    coords = df[['X', 'Y', 'Z']].values
    try:
        n_splits_values = [int(s.strip()) for s in n_splits_str.split(',')]
    except ValueError:
        print(f"Error: n_splits_str must be a comma-separated list of integers (e.g., '2,4,8,16'). Got: {n_splits_str}")
        sys.exit(1)
    for current_n_splits in n_splits_values:
        print(f"--- Processing for {current_n_splits} splits ---")
        labels = create_blocks(coords, current_n_splits)
        # Add block split column to DataFrame
        col_name = f"{current_n_splits}_block_splits"
        df[col_name] = labels
        vis_file = os.path.join(output_vis_dir, f"{subject_id}_blocks_{current_n_splits}splits.nii.gz")
        create_nifti_visualization(coords, labels, vis_file)
    # Save updated DataFrame to CSV
    df.to_csv(uk_data_file_path, index=False)
    print(f"\nBlock NIfTI creation and CSV update completed successfully! Updated CSV: {uk_data_file_path}")

if __name__ == "__main__":
    main() 