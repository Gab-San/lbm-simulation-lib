import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re


# LIST OF REYNOLDS NUMBERS TO ANALYZE
RE_list = [40, 100, 175, 300, 400, 500, 600, 700, 800, 1000, 1250, 1500, 1750, 1950]


def analyze_file(RE, plot_diagnostic=False):
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
    scales = [2, 3]
    found = False
    scale = 2  # default value

    for s in scales:
        filename = "../../benchmarks/strouhal/" + f"Eddy_Generation_RE{RE}_S{s}.csv"
        try:
            scale = s
            data = pd.read_csv(filename)
            print(f"[RE={RE}] File found: {filename}")
            found = True
            break
        except FileNotFoundError:
            continue

    if not found:
        print(f"[RE={RE}] No file found.")
        return None

    # using only the x component for improved stability:
    time = data["Time"].to_numpy(dtype=float)
    if RE < 800:
        signal = data["q1(velocity_y)"].to_numpy(dtype=float)
    else:
        signal = data["q1(velocity_magnitude)"].to_numpy(dtype=float)

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
    id_start = round(len(signal) * 0.4)

    init = round(len(signal) * (45.0/len(signal)))
    end =  round(len(signal) * (130.0/len(signal)))

    if RE <= 400:
        signal = signal[init:end]
    else:
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

    if RE < 800 or RE > 1800:
         St = ((peak_frequency) / 500) * D / U
    else:
         St = ((peak_frequency * 0.5) / 500) * D / U


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
Re_values = []
St_values = []

for RE in RE_list:
    try:
        St = analyze_file(RE, plot_diagnostic=False)
        Re_values.append(RE)
        St_values.append(St)
    except FileNotFoundError:
        print(f"[WARNING] File for RE={RE} not found, skipping this case.")

Re_values = np.array(Re_values, dtype=float)
St_values = np.array(St_values, dtype=float)

print("\n===================================")
print("   SUMMARY St(Re)")
print("===================================")
for re_, st_ in zip(Re_values, St_values):
    print(f"Re = {re_:8.1f}  ->  St = {st_:.6f}")
print("===================================")


# COMPARISON WITH REFERENCE CURVE (LITERATURE)
# Approximate reference values, digitized from the classic
# St-Re diagram for a circular cylinder (Lienhard 1966 / Roshko),
# valid in the range 40 <= Re <= 2e5 (laminar/turbulent wake
# region with laminar boundary layer on the cylinder).
Re_ref = np.array([
    40, 60, 80, 100, 150, 200, 300, 500,
    1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5
])
 
# Lower bound of the band
St_ref_lower = np.array([
    0.120, 0.140, 0.150, 0.160, 0.178, 0.190, 0.200, 0.205,
    0.210, 0.212, 0.210, 0.205, 0.198, 0.190, 0.185, 0.190
])
 
# Upper bound of the band
St_ref_upper = np.array([
    0.130, 0.152, 0.163, 0.172, 0.190, 0.200, 0.208, 0.213,
    0.218, 0.220, 0.218, 0.213, 0.207, 0.200, 0.198, 0.205
])
 

# PLOT St vs Re: computed points vs reference curve
plt.figure(figsize=(10, 6))
 
# Reference band from literature (lower and upper bounds)
plt.plot(Re_ref, St_ref_lower, "-", color="tab:blue",
          label="Lower bound (Lienhard/Roshko)")
plt.plot(Re_ref, St_ref_upper, "--", color="tab:blue",
          label="Upper bound (Lienhard/Roshko)")
plt.fill_between(Re_ref, St_ref_lower, St_ref_upper,
                  color="tab:blue", alpha=0.15)
 
# Computed points from data (simple scatter, no interpolation)
plt.plot(Re_values, St_values, "o", color="tab:red", markersize=8,
          label="Computed data (St from FFT)")

# Polynomial fit of St_values versus Re_values
grade = 5
title_fit = f"Polynomial fit (degree {grade})"
poly_coeffs = np.polyfit(Re_values, St_values, grade)
poly_fit = np.polyval(poly_coeffs, Re_values)
plt.plot(Re_values, poly_fit, "-",
          color=f"tab:green", label=title_fit)
 
# plt.xscale("log")
plt.xlim(20, 2000)
plt.ylim(0.1, 0.25)
 
plt.xlabel("Reynolds number, Re")
plt.ylabel("Strouhal number, St")
plt.title("Strouhal number as a function of Reynolds number")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()
 
plt.tight_layout()
plt.show()