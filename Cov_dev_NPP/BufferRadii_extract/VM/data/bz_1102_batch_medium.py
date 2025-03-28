#Bufferzone_niftispheremasker_adapted
#Adapted from https://github.com/nilearn/nilearn/blob/98a3ee060/nilearn/maskers/nifti_spheres_masker.py#L443

"""
Transformer for computing seeds signals
----------------------------------------
Mask nifti images by spherical volumes for seed-region analyses
"""
import numpy as np
import warnings
from sklearn import neighbors
from joblib import Memory
from scipy import sparse
from nilearn import image, masking
from nilearn._utils import CacheMixin, logger, fill_doc
from nilearn.maskers.base_masker import _filter_and_extract, BaseMasker
from nilearn._utils.class_inspect import get_params
from nilearn._utils.niimg import img_data_dtype
from nilearn._utils.niimg_conversions import (_safe_get_data, check_niimg_4d, check_niimg_3d,)
from numpy import genfromtxt
import os
from os.path import join, exists, split
import nilearn
from nilearn.maskers import NiftiMasker
import nibabel as nib
import numpy as np
import datetime
from datetime import datetime
import pandas as pd
from scipy.sparse import isspmatrix_lil, isspmatrix_csr
import shutil
import gc
import logging

def _apply_mask_and_get_affinity(seeds, niimg, radius, allow_overlap,
                                 mask_img=None):
    """Utility function to get only the rows which are occupied by sphere at
    given seed locations and the provided radius. Rows are in target_affine and
    target_shape space.
    Parameters
    ----------
    seeds : List of triplets of coordinates in native space
        Seed definitions. List of coordinates of the seeds in the same space
        as target_affine.
    niimg : 3D/4D Niimg-like object
        See :ref:`extracting_data`.
        Images to process.
        If a 3D niimg is provided, a singleton dimension will be added to
        the output to represent the single scan in the niimg.
    radius : float
        Indicates, in millimeters, the radius for the sphere around the seed.
    allow_overlap : boolean
        If False, a ValueError is raised if VOIs overlap
    mask_img : Niimg-like object, optional
        Mask to apply to regions before extracting signals. If niimg is None,
        mask_img is used as a reference space in which the spheres 'indices are
        placed.
    Returns
    -------
    X : 2D numpy.ndarray
        Signal for each brain voxel in the (masked) niimgs.
        shape: (number of scans, number of voxels)
    A : scipy.sparse.lil_matrix
        Contains the boolean indices for each sphere.
        shape: (number of seeds, number of voxels)
    """

    # Compute world coordinates of all in-mask voxels.
    if niimg is None:
        mask, affine = masking._load_mask_img(mask_img)
        # Get coordinate for all voxels inside of mask
        mask_coords = np.asarray(np.nonzero(mask)).T.tolist()
        X = None

    elif mask_img is not None:
        affine = niimg.affine
        mask_img = check_niimg_3d(mask_img)
        mask_img = image.resample_img(
            mask_img,
            target_affine=affine,
            target_shape=niimg.shape[:3],
            interpolation='nearest',
        )
        mask, _ = masking._load_mask_img(mask_img)
        mask_coords = list(zip(*np.where(mask != 0)))

        X = masking._apply_mask_fmri(niimg, mask_img)

    mask_coords = np.asarray(list(zip(*mask_coords)))
    mask_coords = image.resampling.coord_transform(
        mask_coords[0], mask_coords[1], mask_coords[2], affine
        )
    mask_coords = np.asarray(mask_coords).T

    clf = neighbors.NearestNeighbors(radius=radius, n_jobs= -1)
    clf.fit(mask_coords)

    A=clf.radius_neighbors_graph(seeds)

    print ('Sparse matrix computed')
    logging.info ('Sparse matrix computed')

    A = A.tolil()

    for i, seed in enumerate(seeds):
        if seed is None:
            continue

        A[i, seed] = True

    print ('Sparse to list computed')
    logging.info('Sparse to list computed')

    sphere_sizes = np.asarray(A.tocsr().sum(axis=1)).ravel()
    print(f'sphere size:{sphere_sizes}')
    logging.info(f'sphere size:{sphere_sizes}')
    return X, A

def _iter_signals_from_spheres(seeds, niimg, radius, allow_overlap,
                               mask_img=None):
    """Utility function to iterate over spheres.
    Parameters
    ----------
    seeds : :obj:`list` of triplets of coordinates in native space
        Seed definitions. List of coordinates of the seeds in the same space
        as the images (typically MNI or TAL).
    niimg : 3D/4D Niimg-like object
        See :ref:`extracting_data`.
        Images to process.
        If a 3D niimg is provided, a singleton dimension will be added to
        the output to represent the single scan in the niimg.
    radius: float
        Indicates, in millimeters, the radius for the sphere around the seed.
    allow_overlap: boolean
        If False, an error is raised if the maps overlaps (ie at least two
        maps have a non-zero value for the same voxel).
    mask_img : Niimg-like object, optional
        See :ref:`extracting_data`.
        Mask to apply to regions before extracting signals.
    """
    X, A = _apply_mask_and_get_affinity(seeds, niimg, radius,
                                        allow_overlap,
                                        mask_img=mask_img)

    for i, row in enumerate(A.rows):
        yield X[:, row]

class _ExtractionFunctor(object):

    func_name = 'nifti_spheres_masker_extractor'

    def __init__(self, seeds_, radius, mask_img, allow_overlap, dtype):
        self.seeds_ = seeds_
        self.radius = radius
        self.mask_img = mask_img
        self.allow_overlap = allow_overlap
        self.dtype = dtype

    def __call__(self, imgs):
        n_seeds = len(self.seeds_)
        imgs = check_niimg_4d(imgs, dtype=self.dtype)

        signals = np.empty(
            (imgs.shape[3], n_seeds),
            dtype=img_data_dtype(imgs)
        )
        print(f'signals_shape:{signals.shape}')
        logging.info(f'signals_shape:{signals.shape}')

        for i, sphere in enumerate(_iter_signals_from_spheres(
                self.seeds_, imgs, self.radius, self.allow_overlap,
                mask_img=self.mask_img)):
            signals[:, i] = np.mean(sphere, axis=1)
        return signals, None

class NiftiSpheresMasker(BaseMasker, CacheMixin):
    """Class for masking of Niimg-like objects using seeds.
    NiftiSpheresMasker is useful when data from given seeds should be
    extracted. Use case: Summarize brain signals from seeds that were
    obtained from prior knowledge.
    Parameters
    ----------
    seeds : :obj:`list` of triplet of coordinates in native space
        Seed definitions. List of coordinates of the seeds in the same space
        as the images (typically MNI or TAL).
    radius : :obj:`float`, optional
        Indicates, in millimeters, the radius for the sphere around the seed.
        Default is None (signal is extracted on a single voxel).
    mask_img : Niimg-like object, optional
        See :ref:`extracting_data`.
        Mask to apply to regions before extracting signals.
    allow_overlap : :obj:`bool`, optional
        If False, an error is raised if the maps overlaps (ie at least two
        maps have a non-zero value for the same voxel). Default=False.
    %(smoothing_fwhm)s
    standardize : {False, True, 'zscore', 'psc'}, optional
        Strategy to standardize the signal.
        'zscore': the signal is z-scored. Timeseries are shifted
        to zero mean and scaled to unit variance.
        'psc':  Timeseries are shifted to zero mean value and scaled
        to percent signal change (as compared to original mean signal).
        True : the signal is z-scored. Timeseries are shifted
        to zero mean and scaled to unit variance.
        False : Do not standardize the data.
        Default=False.
    standardize_confounds : :obj:`bool`, optional
        If standardize_confounds is True, the confounds are z-scored:
        their mean is put to 0 and their variance to 1 in the time dimension.
        Default=True.
    high_variance_confounds : :obj:`bool`, optional
        If True, high variance confounds are computed on provided image with
        :func:`nilearn.image.high_variance_confounds` and default parameters
        and regressed out. Default=False.
    detrend : :obj:`bool`, optional
        This parameter is passed to signal.clean. Please see the related
        documentation for details. Default=False.
    low_pass : None or :obj:`float`, optional
        This parameter is passed to signal.clean. Please see the related
        documentation for details.
    high_pass : None or :obj:`float`, optional
        This parameter is passed to signal.clean. Please see the related
        documentation for details.
    t_r : :obj:`float`, optional
        This parameter is passed to signal.clean. Please see the related
        documentation for details.
    dtype : {dtype, "auto"}, optional
        Data type toward which the data should be converted. If "auto", the
        data will be converted to int32 if dtype is discrete and float32 if it
        is continuous.
    memory : :obj:`joblib.Memory` or :obj:`str`, optional
        Used to cache the region extraction process.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.
    memory_level : :obj:`int`, optional
        Aggressiveness of memory caching. The higher the number, the higher
        the number of functions that will be cached. Zero means no caching.
        Default=1.
    verbose : :obj:`int`, optional
        Indicate the level of verbosity. By default, nothing is printed.
        Default=0.
    Attributes
    ----------
    n_elements_ : :obj:`int`
        The number of seeds in the masker.
        .. versionadded:: 0.9.2
    seeds_ : :obj:`list` of :obj:`list`
        The coordinates of the seeds in the masker.
    See also
    --------
    nilearn.maskers.NiftiMasker
    """
    # memory and memory_level are used by CacheMixin.

    def __init__(
        self,
        seeds,
        radius=None,
        mask_img=None,
        allow_overlap=False,
        smoothing_fwhm=None,
        standardize=False,
        standardize_confounds=True,
        high_variance_confounds=False,
        detrend=False,
        low_pass=None,
        high_pass=None,
        t_r=None,
        dtype=None,
        memory=Memory(location=None, verbose=0),
        memory_level=1,
        verbose=0,
    ):
        self.seeds = seeds
        self.mask_img = mask_img
        self.radius = radius
        self.allow_overlap = allow_overlap

        # Parameters for _smooth_array
        self.smoothing_fwhm = smoothing_fwhm

        # Parameters for clean()
        self.standardize = standardize
        self.standardize_confounds = standardize_confounds
        self.high_variance_confounds = high_variance_confounds
        self.detrend = detrend
        self.low_pass = low_pass
        self.high_pass = high_pass
        self.t_r = t_r
        self.dtype = dtype

        # Parameters for joblib
        self.memory = memory
        self.memory_level = memory_level

        self.verbose = verbose

    def fit(self, X=None, y=None):
        """Prepare signal extraction from regions.
        All parameters are unused; they are for scikit-learn compatibility.
        """
        if hasattr(self, 'seeds_'):
            return self

        error = (
            'Seeds must be a list of triplets of coordinates in '
            'native space.\n'
        )

        if not hasattr(self.seeds, '__iter__'):
            raise ValueError(
                error + 'Given seed list is of type: ' + type(self.seeds)
            )

        self.seeds_ = []
        # Check seeds and convert them to lists if needed
        for i, seed in enumerate(self.seeds):
            # Check the type first
            if not hasattr(seed, '__len__'):
                raise ValueError(
                    error + f'Seed #{i} is not a valid triplet '
                    f'of coordinates. It is of type {type(seed)}.'
                )
            # Convert to list because it is easier to process
            if isinstance(seed, np.ndarray):
                seed = seed.tolist()
            else:
                # in case of tuple
                seed = list(seed)

            # Check the length
            if len(seed) != 3:
                raise ValueError(
                    error + f'Seed #{i} is of length %{len(seed)} '
                    'instead of 3.'
                )

            self.seeds_.append(seed)

        self.n_elements_ = len(self.seeds_)

        return self

    def fit_transform(self, imgs, confounds=None, sample_mask=None):
        """Prepare and perform signal extraction.
        Parameters
        ----------
        imgs : 3D/4D Niimg-like object
            See :ref:`extracting_data`.
            Images to process.
            If a 3D niimg is provided, a singleton dimension will be added to
            the output to represent the single scan in the niimg.
        confounds : CSV file or array-like or :obj:`pandas.DataFrame`, optional
            This parameter is passed to signal.clean. Please see the related
            documentation for details.
            shape: (number of scans, number of confounds)
        sample_mask : Any type compatible with numpy-array indexing, optional
            Masks the niimgs along time/fourth dimension to perform scrubbing
            (remove volumes with high motion) and/or non-steady-state volumes.
            This parameter is passed to signal.clean.
            shape: (number of scans - number of volumes removed, )
                .. versionadded:: 0.8.0
        Returns
        -------
        region_signals : 2D :obj:`numpy.ndarray`
            Signal for each sphere.
            shape: (number of scans, number of spheres)
        """
        return self.fit().transform(imgs, confounds=confounds,
                                    sample_mask=sample_mask)

    def _check_fitted(self):
        if not hasattr(self, "seeds_"):
            raise ValueError('It seems that %s has not been fitted. '
                             'You must call fit() before calling transform().'
                             % self.__class__.__name__)

    def transform_single_imgs(self, imgs, confounds=None, sample_mask=None):
        """Extract signals from a single 4D niimg.
        Parameters
        ----------
        imgs : 3D/4D Niimg-like object
            See :ref:`extracting_data`.
            Images to process.
            If a 3D niimg is provided, a singleton dimension will be added to
            the output to represent the single scan in the niimg.
        confounds : CSV file or array-like or :obj:`pandas.DataFrame`, optional
            This parameter is passed to signal.clean. Please see the related
            documentation for details.
            shape: (number of scans, number of confounds)
        sample_mask : Any type compatible with numpy-array indexing, optional
            Masks the niimgs along time/fourth dimension to perform scrubbing
            (remove volumes with high motion) and/or non-steady-state volumes.
            This parameter is passed to signal.clean.
            shape: (number of scans - number of volumes removed, )
                .. versionadded:: 0.8.0
        Returns
        -------
        region_signals : 2D :obj:`numpy.ndarray`
            Signal for each sphere.
            shape: (number of scans, number of spheres)
        Warns
        -----
        DeprecationWarning
            If a 3D niimg input is provided, the current behavior
            (adding a singleton dimension to produce a 2D array) is deprecated.
            Starting in version 0.12, a 1D array will be returned for 3D
            inputs.
        """
        self._check_fitted()

        params = get_params(NiftiSpheresMasker, self)

        signals,_ = self._cache(
            _filter_and_extract,
            ignore=['verbose', 'memory', 'memory_level'])(
                imgs,
                _ExtractionFunctor(
                    self.seeds_, self.radius, self.mask_img,
                    self.allow_overlap, self.dtype
                ),
                # Pre-processing
                params,
                confounds=confounds,
                sample_mask=sample_mask,
                dtype=self.dtype,
                # Caching
                memory=self.memory,
                memory_level=self.memory_level,
                # kwargs
                verbose=self.verbose
        )
        return signals

#Function to process covariates across all variables
def proc_var(file, ROI, Radii):
    infile = os.path.join(indir,file)
    print(infile)
    logging.info(infile)
    img = nib.load(infile)
    img.set_data_dtype(np.float32)
    print(img.shape)
    logging.info(img.shape)
    Radii.sort(key=lambda x: x[0], reverse = True)
    print(Radii)
    logging.info(Radii)

    for roi in ROI:
        print(roi)
        logging.info(roi)
        for r, chunk_size in Radii:
            print(f'radius:{r} and chunk size:{chunk_size}')
            logging.info(f'radius:{r} and chunk size:{chunk_size}')
            check_empty = signal_extractor(r,chunk_size,roi,var,out_dir,img)

            if check_empty == True:
                print(f'Zero signal found for {var}_{roi}_{r}')
                logging.info(f'Zero signal found for {var}_{roi}_{r}')
                break

            gc.collect()
            print("cache cleared")
            logging.info("cache cleared")

    shutil.move(infile, dst_path)
    print(f'{infile} is processed')
    logging.info(f'{infile} is processed')

def signal_extractor(r,chunk_size,roi,var, out_dir,img):

    #Paramaeters (for more info, refer to https://nilearn.github.io/stable/modules/generated/nilearn.maskers.NiftiSpheresMasker.html#nilearn.maskers.NiftiSpheresMasker )
    rad=r
    std= False
    overlap= "True" 
    smooth= None
    mem_level=1
    mem='nilearn_cache'
    verb=2
    chunk = chunk_size


    #Seed list for coordinates to process
    roi_seed_fl = join(seed_dir,f'ROI_MNI_xyz_coord_{roi}_mask.csv')
    print(f'Variable:{roi},{roi_seed_fl}')
    logging.info(f'Variable:{roi},{roi_seed_fl}')
    data = genfromtxt(roi_seed_fl, delimiter=',', skip_header = 1, dtype = float)
    np.set_printoptions(suppress=True)
    print(data)
    p,q = data.shape
    df = np.empty([p,3])
    for i in range(p):
        df[i][0] = data[i][1]
        df[i][1] = data[i][2]
        df[i][2] = data[i][3]
    print(df.shape)
    logging.info(df.shape)

    coords = np.array_split(df,chunk)
    print(len(coords))
    logging.info(len(coords))

    #Initiate NiftiSphereMasker
    #from nilearn.maskers import NiftiSpheresMasker
    data_se = np.array([])
    counter = 1
    check_empty = False

    for coord in coords:
        print(len(coord))
        logging.info(len(coord))
        seeds = np.array(coord).tolist()

        masker = NiftiSpheresMasker(
            seeds, radius=rad, detrend=False, standardize=std, 
            allow_overlap= overlap, smoothing_fwhm = smooth, mask_img = mask,
            memory_level=mem_level, memory=mem, verbose=verb, dtype= "auto")

        #Extract signal
        signal_ext = masker.fit_transform(img)
        signal_ext
        data_se = np.array(np.concatenate((data_se, signal_ext[0])))
        # print(len(data_se))
        # print(data_se)
        print(f'chunk processed: {counter}')
        logging.info(f'chunk processed: {counter}')
        counter = counter + 1

    chk_mean = np.mean(data_se)
    print(f'Mean signal = {chk_mean}')
    logging.info(f'Mean signal = {chk_mean}')

    if chk_mean < 0.01:
        check_empty = True

    df_se =pd.DataFrame(data_se)
    print(df_se.shape)
    logging.info(df_se.shape)

    #Write data to file
    out_fol = join(out_dir,var)
    out_fl_name = f"bz_{var}_{roi}_{rad}.csv"
    print(out_fl_name)
    logging.info(out_fl_name)
    df_se.to_csv(join(out_fol,out_fl_name))

    del signal_ext, df_se, data_se, counter
    del coords, seeds

    return check_empty

warnings.filterwarnings('ignore')
t_start = datetime.today()

#Set-up log files
LOG_FILENAME = datetime.now().strftime('logfile_%H_%M_%S_%d_%m_%Y.log')
logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)

#Input
homedir = join((split(os.getcwd())[0]),"data")
print(homedir)
logging.info(homedir)
indir = join(homedir,"Input/Variables")
mask_file = join(homedir,"Input/Masks/MNI152b_brain_mask_0.1.nii")
print (mask_file)
logging.info (mask_file)
mask = nib.load(mask_file)
dst_path = join(homedir,"Input/Processed")
seed_dir = join(homedir,"Input/Seeds")

ROI = ("AM","HC","IPL","MFG","SMTG")
Radii = [[12.5,2],[15,4],[20,6]]

#Output
out_dir = join(homedir,"output")

files = os.listdir(indir)
files.sort()
print(files)
logging.info(files)
print(f'Number of files to be processed:{len(files)})')
logging.info(f'Number of files to be processed:{len(files)})')

#Loop through all variable files
index = 0
for file in files:
    name = file.split('.nii')[0]
    print(name)
    logging.info(name)
    var = name
    out_fldr = join(out_dir,name)
    if not os.path.isdir(out_fldr):
        os.makedirs(out_fldr)

    proc_var(file, ROI, Radii)
    index = index+1
    print(f'Number of files processed:{index}')
    logging.info(f'Number of files processed:{index}')

t_end = datetime.today()
t_run = t_end - t_start
print(f'total runtime is {t_run}')
logging.info(f'total runtime is {t_run}')