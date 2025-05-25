#!/usr/bin/env python
"""
Voltage Imaging Analysis Pipeline

This script performs a complete pipeline for processing voltage imaging data
including motion correction, segmentation, trace denoising, spike extraction,
and result saving. Dataset format must be compatible with CaImAn/VolPy.

Adapted and cleaned for Electra voltage imaging data analysis.
"""

import os
import glob
import h5py
import logging
import numpy as np
import matplotlib.pyplot as plt
import caiman as cm
from caiman.motion_correction import MotionCorrect
from caiman.source_extraction.volpy.volparams import volparams
from caiman.source_extraction.volpy.volpy import VOLPY
from caiman.summary_images import mean_image, local_correlations_movie_offline
from scipy.io import savemat

def clean_dict(d):
    """Recursively replace None with empty lists."""
    return {k: clean_dict(v) if isinstance(v, dict) else ([] if v is None else v) for k, v in d.items()}

def save_dict_to_matfile(d, path):
    savemat(path, clean_dict(d))

def run_volpy_pipeline(image_file, roi_file, flip_signal):
    file_dir = os.path.dirname(image_file)
    fr = 1000  # Frame rate (Hz)

    # Motion correction parameters
    opts_dict = {
        'fnames': image_file,
        'fr': fr,
        'pw_rigid': False,
        'max_shifts': (2, 2),
        'gSig_filt': (3, 3),
        'strides': (48, 48),
        'overlaps': (24, 24),
        'max_deviation_rigid': 3,
        'border_nan': 'copy'
    }
    opts = volparams(params_dict=opts_dict)

    # Set up cluster
    c, dview, n_processes = cm.cluster.setup_cluster(backend='multiprocessing', n_processes=None, single_thread=False)

    # Motion correction
    mc = MotionCorrect(image_file, dview=dview, **opts.get_group('motion'))
    mc.motion_correct(save_movie=True)

    # Memory mapping
    border_to_0 = 0 if mc.border_nan == 'copy' else mc.border_to_0
    base_name = 'memmap_' + os.path.splitext(os.path.basename(image_file))[0]
    fname_new = cm.save_memmap_join(mc.mmap_file, base_name=base_name, add_to_mov=border_to_0, dview=dview)

    # Segmentation - summary images
    img_mean = mean_image(mc.mmap_file[0], window=1000, dview=dview)
    img_corr = local_correlations_movie_offline(mc.mmap_file[0], fr=fr, window=fr*4, stride=fr*4, 
                                                 winSize_baseline=fr, remove_baseline=True, 
                                                 gaussian_blur=False, dview=dview).max(axis=0)

    img_norm = lambda x: (x - np.mean(x)) / np.std(x)
    summary_images = np.stack([img_norm(img_mean)] * 2 + [img_norm(img_corr)], axis=0).astype(np.float32)
    cm.movie(summary_images).save(image_file[:-5] + '_summary_images.tif')

    # Load manual ROIs
    with h5py.File(roi_file, 'r') as f:
        ROIs = f['mov'][()]
        if ROIs.ndim == 2:
            ROIs = np.expand_dims(ROIs, axis=0)

    # Restart cluster to free memory
    cm.stop_server(dview=dview)
    c, dview, n_processes = cm.cluster.setup_cluster(backend='multiprocessing', n_processes=None, single_thread=False, maxtasksperchild=1)

    # Trace extraction parameters
    opts_dict = {
        'fnames': fname_new,
        'ROIs': ROIs,
        'index': list(range(len(ROIs))),
        'weights': None,
        'template_size': 0.02,
        'context_size': 35,
        'visualize_ROI': False,
        'flip_signal': flip_signal,
        'hp_freq_pb': 1/30,
        'clip': 100,
        'threshold_method': 'adaptive_threshold',
        'min_spikes': 2,
        'pnorm': 0.5,
        'threshold': 3,
        'do_plot': False,
        'ridge_bg': 0.01,
        'sub_freq': 20,
        'weight_update': 'ridge',
        'n_iter': 2
    }
    opts.change_params(params_dict=opts_dict)

    # Run VolPy
    vpy = VOLPY(n_processes=n_processes, dview=dview, params=opts)
    try:
        vpy.fit(n_processes=n_processes, dview=dview)
    except:
        cm.stop_server(dview=dview)
        for log_file in glob.glob('*_LOG_*'):
            os.remove(log_file)
        return

    # Save results
    estimates = vpy.estimates
    estimates['ROIs'] = ROIs
    estimates['params'] = opts

    save_name = f"volpy_{os.path.basename(image_file)[:-3]}_adaptive_threshold"
    np.save(os.path.join(file_dir, save_name), estimates)

    keys_to_save = ['cell_n', 't', 'ts', 't_rec', 't_sub', 'spikes', 'templates', 'snr', 'thresh',
                    'context_coord', 'F0', 'dFF', 'ROIs']
    reduced_estimates = {k: estimates[k] for k in keys_to_save if k in estimates}
    savemat(os.path.join(file_dir, save_name + '.mat'), reduced_estimates)

    cm.stop_server(dview=dview)
    for log_file in glob.glob('*_LOG_*'):
        os.remove(log_file)

def get_flip_signal(folder):
    if "VS73" in folder or "VS89" in folder:
        return False
    elif "VS74" in folder or "VS91" in folder:
        return True
    else:
        raise ValueError(f"Unknown flip signal in folder: {folder}")

def main():
    base_dirs = ["<your_data_path_1>", "<your_data_path_2>"]  # <-- Replace with actual paths
    image_files, roi_files, flip_signals = [], [], []

    for base_dir in base_dirs:
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if file.endswith(".h5") and not file.endswith("_avg_projection.hdf5"):
                    image_path = os.path.join(root, file)
                    roi_path = image_path.replace('.h5', '_avg_projection.hdf5')
                    image_files.append(image_path)
                    roi_files.append(roi_path)
                    flip_signals.append(get_flip_signal(os.path.basename(root)))

    for image_file, roi_file, flip_signal in zip(image_files, roi_files, flip_signals):
        run_volpy_pipeline(image_file, roi_file, flip_signal)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    main()
