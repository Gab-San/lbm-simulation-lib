import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re


# LIST OF SCALES TO ANALYZE
scales = [1, 2, 3, 4, 5, 6, 7]
RE = 1000


def analyze_file(scale, plot_diagnostic=False):
    """
    Reads the CSV file for a given Reynolds number,
    computes the FFT of the signal and returns the Strouhal
    number calculated with the original formula.

    If plot_diagnostic=True, also shows time-domain and
    FFT spectrum plots for that single case.
    """

    # --------------------------------------------------------
    # 1. CSV FILE READING
    # --------------------------------------------------------
    filename = "../../benchmarks/strouhal/" + f"Eddy_Generation_Signal{scale}.csv"
    data = pd.read_csv(filename)

    # using only the x component for improved stability:
    time = data["Time"].to_numpy(dtype=float)
    signal = data["q1(velocity_y)"].to_numpy(dtype=float)

    N = len(signal)

    # --------------------------------------------------------
    # 2. SAMPLING FREQUENCY
    # --------------------------------------------------------
    dt = np.mean(np.diff(time))
    dt_variation = np.max(np.abs(np.diff(time) - dt))

    if dt_variation > 1e-10:
        print(f"[RE={RE}] WARNING: sampling is not perfectly uniform.")

    Fs = 1.0 / dt

    print(f"\n--- RE = {RE} ---")
    print(f"Number of samples: {N}")
    print(f"Sampling interval dt: {dt:.6f}")
    print(f"Sampling frequency Fs: {Fs:.6f}")

    # --------------------------------------------------------
    # 3. MEAN SUBTRACTION (after transient, 40%)
    # --------------------------------------------------------
    id_start = round(len(signal) * 0.1)
    signal = signal[id_start:]

    N = len(signal)

    signal_mean = np.mean(signal)
    signal_detrended = signal - signal_mean

    print(f"Signal mean value: {signal_mean:.6e}")

    # --------------------------------------------------------
    # 4. HANN WINDOW
    # --------------------------------------------------------
    window = np.hanning(N)
    signal_windowed = signal_detrended * window

    # --------------------------------------------------------
    # 5. TIME-DOMAIN SIGNAL PLOT (optional)
    # --------------------------------------------------------
    if plot_diagnostic:
        plt.figure(figsize=(10, 5))
        plt.plot(signal, label="Original signal", alpha=0.7)
        plt.plot(signal_windowed, label="Signal after mean subtraction + Hann window")
        plt.xlabel("Time [LBM]")
        plt.ylabel("Velocity magnitude")
        plt.title(f"Time-domain signal (RE={RE})")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    # --------------------------------------------------------
    # 6. FFT COMPUTATION
    # --------------------------------------------------------
    fft_signal = np.fft.rfft(signal_windowed)
    frequencies = np.fft.rfftfreq(N, d=dt)

    # --------------------------------------------------------
    # 7. AMPLITUDE SPECTRUM
    # --------------------------------------------------------
    coherent_gain = np.mean(window)
    amplitude = np.abs(fft_signal) / (N * coherent_gain)

    if N % 2 == 0:
        amplitude[1:-1] *= 2
    else:
        amplitude[1:] *= 2

    # --------------------------------------------------------
    # 8. PEAK DETECTION
    # --------------------------------------------------------
    min_bin = 2
    peak_index = np.argmax(amplitude[min_bin:]) + min_bin

    peak_frequency = frequencies[peak_index]
    peak_amplitude = amplitude[peak_index]

    print("FFT MAIN PEAK")
    print(f"Frequency: {peak_frequency:.8f} [1/time]")
    print(f"Amplitude:  {peak_amplitude:.6e}")

    # --------------------------------------------------------
    # 9. STROUHAL NUMBER CALCULATION (formula unchanged)
    # --------------------------------------------------------
    D = 16 * 2 * scale
    U = 0.1

    St = (peak_frequency / 500) * D / U


    print(f"D = {D}")
    print(f"U = {U}")
    print(f"St = {St:.6f}")

    # --------------------------------------------------------
    # 10. FFT SPECTRUM PLOT (optional)
    # --------------------------------------------------------
    if plot_diagnostic:
        plt.figure(figsize=(10, 5))
        plt.plot(frequencies, amplitude, label="FFT (mean subtracted + Hann)")
        plt.plot(peak_frequency, peak_amplitude, "ro",
                 label=f"Peak = {peak_frequency:.6f}")
        plt.xlabel("Frequency [1/time]")
        plt.ylabel("Amplitude")
        plt.title(f"FFT spectrum of the signal (RE={RE})")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return St



# LOOP OVER ALL REYNOLDS NUMBERS
Scale_values = []
St_values = []

for s in scales:
    try:
        St = analyze_file(s, plot_diagnostic=False)
        Scale_values.append(s)
        St_values.append(St)
    except FileNotFoundError:
        print(f"[WARNING] File for RE={s} not found, skipping this case.")

scales = np.array(scales, dtype=float)
Scale_values = np.array(Scale_values, dtype=float)
St_values = np.array(St_values, dtype=float)

print("\n===================================")
print("   SUMMARY St(Scale)")
print("===================================")
for re_, st_ in zip(Scale_values, St_values):
    print(f"Re = {re_:8.1f}  ->  St = {st_:.6f}")
print("===================================")

# Lower bound of the band
St_ref_low = 0.210
# Upper bound of the band
St_ref_up = 0.218
 

# PLOT St vs Re: computed points vs reference curve
plt.figure(figsize=(10, 6))
 
# Reference band from literature (lower and upper bounds)
St_ref_lower = np.full_like(scales, St_ref_low)
St_ref_upper = np.full_like(scales, St_ref_up)
plt.plot(scales, St_ref_lower, "-", color="tab:blue",
          label="Lower bound (Lienhard/Roshko)")
plt.plot(scales, St_ref_upper, "--", color="tab:blue",
          label="Upper bound (Lienhard/Roshko)")
plt.fill_between(scales, St_ref_lower, St_ref_upper,
                  color="tab:blue", alpha=0.15)
 
# Computed points from data (simple scatter, no interpolation)
plt.plot(Scale_values, St_values, "o", color="tab:red", markersize=8,
          label="Computed data (St from FFT)")

# Polynomial fit of St_values versus Re_values
grade = 1
title_fit = f"Polynomial fit (degree {grade})"
poly_coeffs = np.polyfit(Scale_values, St_values, grade)
poly_fit = np.polyval(poly_coeffs, Scale_values)
plt.plot(Scale_values, poly_fit, "-",
          color=f"tab:green", label=title_fit)
 
# plt.xscale("log")
# plt.xlim(20, 2000)
plt.ylim(0.15, 0.25)
 
plt.xlabel("Scale factor")
plt.ylabel("Strouhal number, St")
plt.title("Strouhal number as a function of Scale factor")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()
 
plt.tight_layout()
plt.show()