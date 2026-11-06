# combining the definitions in one single fucntin do we can call it in another notebook
import h5py
from matplotlib import pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from pathlib import Path
import json


def regions_above_threshold(data, threshold):
    """Finds start and end indices of contiguous regions where data is above threshold."""
    mask = data > threshold
    segments = []
    in_segment = False
    segment_start = -1
    for k, val in enumerate(mask):
        if val and not in_segment:
            segment_start = k
            in_segment = True
        elif not val and in_segment:
            segments.append((segment_start, k - 1))
            in_segment = False
    if in_segment: # Handle case where data ends above the threshold
        segments.append((segment_start, len(data) - 1))
    return segments


    
# to dos: 
# - find dla crossing points above the threshold
# - plot the peaks above threshold with the exact pixel count
# - store the data in a list/array
# - plot the spctra with the detected peaks using the peak finder when finidng the crossing pintts
# - store all the peaks data like the id number, spectra id, x-coord (pixel count) and y-coord (log_tau) in json format
    
    


def arg_max_DLA_finding(file_path, n_spectra, threshold=6, save_json=False):
    # file_path_add = Path('/scratch/cavestru_root/cavestru0/mfho/priya/emu_full/')
    # file_path = file_path_add / file_path / path_append
    # fig, ax = plt.subplots(1,1, figsize=(15,6))

   
    all_peaks = []
    all_crossing_points = [] # To collect all crossing indices
    total_crossings = 0 # counter for total crossings

    with h5py.File(file_path,'r') as f:
        for i in range(n_spectra):
            tau_spectrum = f['tau/H/1/1215'][i,:]
            log_tau_spectrum = np.log10(tau_spectrum)
            
            signal_diff_from_threshold = log_tau_spectrum - 6

            sign_changes = np.diff(np.sign(signal_diff_from_threshold))

            upward_crossing_indices, _ = find_peaks(sign_changes, height=1.5)
            downward_crossing_indices, _ = find_peaks(-sign_changes, height=1.5) #found problem...implement arg_max

            crossing_indices = np.concatenate((upward_crossing_indices, downward_crossing_indices))
            crossing_indices = np.sort(crossing_indices)

            all_crossing_points.extend(crossing_indices) # Collect indices

            total_crossings += len(crossing_indices) # total number of crissing indices
            
            
            regions = regions_above_threshold(log_tau_spectrum, threshold) # to find peaks above the threshold of 6
            
            for start, end in regions:
                segment_data = log_tau_spectrum[start:end+1]
                peak_index_in_segment = np.argmax(segment_data)
                global_peak_index = start + peak_index_in_segment
                peak_value = segment_data[peak_index_in_segment]
                
                
                
                all_peaks.append({
                    'spectrum_id': i,
                    'peak_pixel_count': global_peak_index,
                    'peak_log_tau': peak_value,
                    'log_tau_spectrum_data': log_tau_spectrum # Store the entire spectrum for plotting
                })
                
             # ploting
            
    fig, ax = plt.subplots(1, 1, figsize=(15, 6))
    plotted_spectrum_ids = set() # To ensure each unique DLA-containing spectrum is plotted only once
    
    for dla_data in all_peaks:
        spectrum_id = dla_data['spectrum_id']
        if spectrum_id not in plotted_spectrum_ids:
            # Plot the entire spectrum if it hasn't been plotted yet
            ax.plot(np.arange(len(dla_data['log_tau_spectrum_data'])), dla_data['log_tau_spectrum_data'],
                    alpha=0.5, linewidth=1) # can add label if needed for the first plot 
            
            plotted_spectrum_ids.add(spectrum_id)

        # Mark the specific peak on the plotted spectrum (the red dts)
        ax.plot(dla_data['peak_pixel_count'], dla_data['peak_log_tau'], 'ro', markersize=5)

    ax.axhline(y=threshold, color='black', linestyle='--', label=f'Threshold (y={threshold})')
    ax.set_ylabel(r'$log_{10}(\tau)$', fontsize=20)
    ax.set_xlabel('line of sight, Pixel counts', fontsize=15)
    ax.set_title(f'Spectra with Detected DLAs (Threshold={threshold})')
    plt.savefig("Spectra with detected DLAs.pdf")
    
    # Print a summary of DLAs found, similar to output
    if len(all_peaks) > 0:
        pass
    else:
        print("No DLAs found.")

    print(f"\nNumber of spectra: {n_spectra}")
    print(f"Number of spectra with DLAs: {len(plotted_spectrum_ids)}")
    print(f"Number of DLAs detected: {len(all_peaks)}")
    
    #json dict
    dla_dict = []
    for i, dla in enumerate(all_peaks):
        dla_data = {
            'S.No': i + 1,
            'spect_id': dla['spectrum_id'],
            'x_coord': int(dla['peak_pixel_count']),
            'y_coord': float(dla['peak_log_tau'])
        }
        dla_dict.append(dla_data)


    print(json.dumps(dla_dict, indent=2))
    
    return {
        'all_peaks': all_peaks,
        'dla_dict': dla_dict,
        'total_crossings': total_crossings,
        'spectra_with_dlas': len(plotted_spectrum_ids)
     }
# example lpt
if __name__ == "__main__":
    file_path = 'ns0.859Ap1.29e-09herei3.92heref2.72alphaq1.87hub0.693omegamh20.141hireionz7.15bhfeedback0.0579/' 
    SNAP = 11
    path_append = 'output/' + f'SPECTRA_{SNAP:03d}' + '/lya_forest_spectra_grid_480.hdf5'
    file_path_add = Path('/scratch/cavestru_root/cavestru0/mfho/priya/emu_full/')
    file_path = file_path_add / file_path / path_append
    
    # Run analysis
    results = arg_max_DLA_finding(
        file_path=file_path,
        n_spectra=20,
        threshold=6,
        save_json=True
    )
    
    dla_data = results['dla_dict']
    peaks = results['all_peaks']
    print(f"Retrieved {len(peaks)} peaks from analysis")

