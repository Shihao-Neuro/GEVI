# Voltage Imaging Post-Analysis for Cultured Neurons

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat

# -----------------------------------------------------------------------------
# File Path Management
# -----------------------------------------------------------------------------
def get_file_paths():
    """Return grouped file paths for each voltage imaging dataset."""
    base_path = "<your_base_directory>"
    return {
        "VS73": [f"{base_path}/VS73/volpy_VS73-2001_adaptive_threshold.mat"],
        "VS74": [f"{base_path}/VS74/volpy_VS74-1001_adaptive_threshold.mat"],
        "VS89": [f"{base_path}/VS89/volpy_VS89-2001_adaptive_threshold.mat"],
        "VS91": [f"{base_path}/VS91/volpy_VS91-1001_adaptive_threshold.mat"]
    }

# -----------------------------------------------------------------------------
# Core Signal Processing Functions
# -----------------------------------------------------------------------------
def load_voltage_data(file_path):
    data = loadmat(file_path)
    return data['dFF'], data['spikes'].flatten()

def compute_snr(trace):
    trace = trace - np.median(trace)
    noise = np.sqrt(np.mean((-trace[trace < 0])**2))
    return trace / noise

def compute_individual_snr(snr, spikes):
    return snr[0, spikes - 1]

def export_mat(snr, spikes, individual_snr, output_path):
    savemat(output_path, {
        'snr': snr,
        'spike_positions': spikes,
        'individual_snr': individual_snr
    })

# -----------------------------------------------------------------------------
# Spike Waveform Analysis
# -----------------------------------------------------------------------------
def extract_waveforms(datasets, window=100, snr_threshold=None):
    """Extract and aggregate spike waveforms from time windows."""
    segments = [(1, 300000), (600001, 900000), (1200001, 1500000)]
    waveforms = [[] for _ in segments]

    for file_path in datasets:
        dFF, spikes = load_voltage_data(file_path)
        snr = compute_snr(dFF)
        spike_snr = compute_individual_snr(snr, spikes)
        spike_dict = dict(zip(spikes, spike_snr))

        for i, (start, end) in enumerate(segments):
            valid_spikes = [s for s in spikes if start <= s < end and s > window and (snr_threshold is None or spike_dict[s] >= snr_threshold)]
            for s in valid_spikes:
                wave = snr[0, s - window:s + window]
                waveforms[i].append(wave)

    return [np.array(w) for w in waveforms]

def compute_stats(waves):
    stats = {'mean': [], 'std': [], 'sem': []}
    for segment in waves:
        stats['mean'].append(np.mean(segment, axis=0))
        stats['std'].append(np.std(segment, axis=0))
        stats['sem'].append(np.std(segment, axis=0) / np.sqrt(len(segment)))
    return stats

def plot_waveforms(stats, window, save_path=None, title_prefix=""):
    time_axis = np.arange(-window, window)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    titles = ['1-300k', '600k-900k', '1.2M-1.5M']

    for i in range(3):
        axs[i].plot(time_axis, stats['mean'][i], label='Mean', color='blue')
        axs[i].fill_between(time_axis,
                            stats['mean'][i] - stats['sem'][i],
                            stats['mean'][i] + stats['sem'][i],
                            alpha=0.3, label='SEM', color='blue')
        axs[i].set_title(f"{title_prefix} {titles[i]}")
        axs[i].set_xlabel("Time (frames)")
        axs[i].set_ylabel("SNR")
        axs[i].legend()
        axs[i].grid(True)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

# -----------------------------------------------------------------------------
# Main Execution Block
# -----------------------------------------------------------------------------
def process_group(group_name, files, output_dir, window=100, snr_thresh=2.0):
    os.makedirs(output_dir, exist_ok=True)

    # SNR & spike export
    for i, file in enumerate(files):
        dFF, spikes = load_voltage_data(file)
        snr = compute_snr(dFF)
        individual_snr = compute_individual_snr(snr, spikes)
        export_mat(snr, spikes, individual_snr, os.path.join(output_dir, f"{group_name}_{i}.mat"))

    # Waveform analysis
    waveforms = extract_waveforms(files, window, snr_thresh)
    stats = compute_stats(waveforms)
    plot_waveforms(stats, window, save_path=os.path.join(output_dir, f"{group_name}_spike_waveform.png"), title_prefix=group_name)

    # Save waveforms to .mat
    savemat(os.path.join(output_dir, f"{group_name}_spike_waveform.mat"), {
        'spike_wave_form': np.array(stats['mean']),
        'spike_std': np.array(stats['std']),
        'spike_sem': np.array(stats['sem']),
        'individual_spikes': waveforms
    })

if __name__ == "__main__":
    output_dir = "<your_output_directory>"  # Replace with actual path
    file_groups = get_file_paths()

    for group_name, files in file_groups.items():
        process_group(group_name, files, output_dir, window=100, snr_thresh=2.0)
