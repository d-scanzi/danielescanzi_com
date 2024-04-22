#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 17:38:58 2024

@author: daniele
"""

import mne
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

## SET UP
data_directory = "C:\\Users\\dsca347\\OneDrive - The University of Auckland\\PhD\\presentations\\labs\\2024\\lab_meeting_2024_04_22\\data"
file_name      = "example_data.set"

raw_data       = mne.io.read_raw_eeglab(os.path.join(data_directory, file_name), verbose=0)
raw_data.crop(tmin=60, tmax=350)


## 1 - DOWNSAMPLING

# Recorded signal
original_sampling_rate = 100
signal_duration_s      = 2
samples_vector         = np.arange(0, signal_duration_s, 1/original_sampling_rate)
sine_signal = 2*np.sin(2*np.pi*2*samples_vector)

# Downsample (MODIFY THE NEXT LINE)
new_sampling_rate = 50

number_of_points_to_retain = int((len(sine_signal)/original_sampling_rate) * new_sampling_rate)
step_downsampled_signal    = int(len(sine_signal) / number_of_points_to_retain)
downsampled_signal         = sine_signal[np.arange(0, len(sine_signal), step_downsampled_signal)]
downsampled_time           = samples_vector[np.arange(0, len(sine_signal), step_downsampled_signal)]

plt.plot(samples_vector, sine_signal, "-k", linewidth=8, alpha=0.3)
plt.scatter(downsampled_time, downsampled_signal, c="r", s=6, alpha=0.6)
plt.plot(downsampled_time, downsampled_signal, "-r", linewidth=1)
plt.rcParams["figure.figsize"] = [10, 5]
plt.show()

## 2 - FILTERING

# Create signal as sum of sines
sampling_rate     = 1000
signal_duration_s = 2
times             = np.arange(0, signal_duration_s, 1/sampling_rate)

signal_frequencies = np.arange(10, 25, 10)
signal_sines       = np.sin(2*np.pi*signal_frequencies[0]*times)
for frequency in signal_frequencies[1:]:
    signal_sines = np.vstack((signal_sines, np.sin(2*np.pi*frequency*times)))

composite_signal = np.sum(signal_sines, axis=0)

# Define Butterworth lowpass filter
polinomial_numerator, polinomial_denominator = signal.butter(10, 15, fs=sampling_rate)
filtered_half_signal = signal.filtfilt(polinomial_numerator, polinomial_denominator, composite_signal[int(len(composite_signal)/2):])

fig, ax = plt.subplots()
ax.plot(times, composite_signal, "k", linewidth=2, alpha=0.5)
ax.plot(times, np.append(composite_signal[0:int(len(composite_signal)/2)], filtered_half_signal), "r")
ax.vlines(times[int(len(times)/2)], min(composite_signal), max(composite_signal), colors="k", linestyles="dashed", linewidth=3)

# Replace second half of signal with filtered version
composite_signal[int(len(composite_signal)/2):] = filtered_half_signal


fft_window = signal.windows.blackman(500)
fft_hop    = int(len(fft_window)/20)
SFT = signal.ShortTimeFFT(fft_window, hop=fft_hop, fs=sampling_rate, scale_to="magnitude", mfft=2**11)
Sx = SFT.stft(composite_signal)  # perform the STFT

fig,ax = plt.subplots()
ax.imshow(abs(Sx), origin="lower", aspect="auto", cmap="viridis", extent=SFT.extent(len(composite_signal)))
ax.hlines((10, 20), 0, (2, 1), color="gray", linewidth=3, alpha=0.5)
ax.vlines(1, 0, 30, colors="r", linestyles="dashed", linewidth=3)
ax.set_ylim((0, 30))
ax.set_xlim((0, 2))

### Roll off and ringing

# Butterworth
cutoff_frequency = 50
filter_orders    = [2, 3, 4, 8, 10, 12]


fig, axs = plt.subplots(2, 1)
plt.tight_layout()
for filter_order in filter_orders:
    b,a = signal.butter(filter_order, cutoff_frequency, "low", analog=True)
    w,h = signal.freqs(b, a)
    t, y = signal.impulse((b,a))
    axs[0].plot(w, 20 * np.log10(abs(h)), label=f"Order: {filter_order}")
    axs[0].vlines(50, -50, 5, colors="k")
    axs[0].set_ylabel("Amplitude (dB)")
    axs[0].set_xlabel("Frequency (Hz)")
    axs[0].set_xlim((0, 100))
    axs[0].set_ylim((-50, 5))
    axs[0].legend()
    
    axs[1].plot(t, y)
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Amplitude (AU)")

# FIR
filter_types  = ["hamming", "kaiser", "blackman", "rectangular"]
filter_lenght = 80
fig, ax = plt.subplots()
for filter_type in filter_types:
    if filter_type != "kaiser":
        fir_coefficients = signal.firwin(filter_lenght, 50, fs=1000, window=(filter_type))
    else:
        fir_coefficients = signal.firwin(filter_lenght, 50, fs=1000, window=(filter_type, 8))

    
    fir_w, fir_h = signal.freqz(fir_coefficients, fs=1000)
    ax.plot(fir_w, 20 * np.log10(abs(fir_h)), label=filter_type)
    ax.set_xlim((0, 200))
    ax.legend()


# Filter components
sample_rate      = 1000
cutoff_frequency = 50
trans_bandwidth  = 10
filter_lenght_s  = 0.2 # Longer time = Steeper filter
filter_order     = int(round((sample_rate * filter_lenght_s))) + 1
kaiser_coefficients = signal.firwin(filter_order, 50, fs=1000, window=("blackman"))
kaiser_w, kaiser_h  = signal.freqz(kaiser_coefficients, fs=1000)

fig, ax = plt.subplots()
ax.plot(kaiser_w, 20 * np.log10(abs(kaiser_h)))
ax.fill_between((cutoff_frequency, cutoff_frequency+trans_bandwidth), 5, -110, alpha=0.35, color='red')
ax.set_xlim((0, 150))
ax.set_ylim((-110, 5))

# Filter EEG data
raw_filtered = raw_data.copy().filter(l_freq=0.5, h_freq=None)
raw_filtered = raw_filtered.filter(h_freq=30, l_freq=None)


## 3 - Find bad channels
raw_filtered.plot()