#Script to apply MNI mask for extraction
import nilearn
from nilearn.masking import apply_mask

nifti2extract = "/Users/raghav_m/Desktop/NeuroPathPredict/input_files/nifti/mni_icbm152_t1_tal_nlin_sym_09b_hires"
brain_mask = "/Users/raghav_m/Desktop/NeuroPathPredict/input_files/nifti/mni_icbm152_t1_tal_nlin_sym_09a_mask"

masked_data = apply_mask(nifti2extract, brain_mask)

nilearn.plotting.plot_anat(anat_img=masked_data)