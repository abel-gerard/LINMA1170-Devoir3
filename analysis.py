import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import os
import sys


if __name__ == "__main__":
    log = "log.txt"

    if not os.path.exists(log):
        print("Log file not found. Please run deformation to create it.")
        exit(1)

    if len(sys.argv) == 1:
        print("Missing argument.")
        exit(1)

    with open(log, "r") as f:
        lines = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

    analysis = sys.argv[1]
    if analysis == "energy":
        begin = 0
        end = -1

        if len(sys.argv) >= 3:
            try:
                begin = float(sys.argv[2])
            except:
                print("Expected float.")
                exit(1)
            if not 0 <= begin <= 1:
                print("begin parameter must be in the range [0,1].")
                exit(1)
            begin = int(begin*len(lines))

        if len(sys.argv) >= 4:
            try:
                end = float(sys.argv[3])
            except:
                print("Expected float.")
                exit(1)
            if not 0 <= end <= 1:
                print("end parameter must be in the range [0,1].")
                exit(1)
            end = int(end*len(lines))

        t = lines[begin:end, 0]
        Ep = lines[begin:end, 1]
        Ec = lines[begin:end, 2]
        E = lines[begin:end, 3]

        plt.plot(t, Ep, label="$E_p$", color="#43CCD0", linewidth=2, linestyle="--")
        plt.plot(t, Ec, label="$E_c$", color="#D043CC", linewidth=2, linestyle="--")
        plt.plot(t, E, label="$E$", color="#AAD043", linewidth=2)

        plt.xlabel("Time (s)")
        plt.ylabel("Energy (J)")
        plt.title("Energy vs Time")
        plt.legend()
        plt.grid()

        plt.show()
 
    elif analysis == "animation":
        pass 

    elif analysis == "frequency":
        t = lines[:, 0]
        Ep = lines[:, 1]
        Ec = lines[:, 2]
        E = lines[:, 3]

        dt = t[1] - t[0]
        N = len(t)

        mean = np.mean(Ec)
        print(f"Mean kinetic energy: {mean} [J]")
        spectrum = np.abs(fft(Ec-mean))[:N//2]
        frequencies = fftfreq(N, dt)[:N//2]

        plt.plot(frequencies, spectrum, color="#43CCD0", linewidth=2)
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude (J)")
        plt.title("Spectrum of (Detrended) Kinetic Energy")
        plt.grid()
        plt.show()

        with open(sys.argv[2], "r") as f:
            moves = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

        t = moves[:, 0]
        x = moves[:, 1]
        y = moves[:, 2]
        vx = moves[:, 3]
        vy = moves[:, 4]

        mean_x = np.mean(x)
        print(f"Mean x position of node 0: {mean_x} [m]")
        spectrum_x = np.abs(fft(x-mean_x))[:N//2]
        frequencies = fftfreq(N, dt)[:N//2]

        mean_y = np.mean(y)
        print(f"Mean y position of node 0: {mean_y} [m]")
        spectrum_y = np.abs(fft(y-mean_y))[:N//2]

        plt.plot(frequencies, spectrum_x, color="#D04394", linewidth=2, label="x")
        plt.plot(frequencies, spectrum_y, color="#43D07F", linewidth=2, label="y")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude (m)")
        plt.title("Spectrum of (Detrended) Position of Node 0")
        plt.legend()
        plt.grid()
        plt.show()

    else:
        print("Unrecognized analysis option, either 'energy', 'frequency'.")
        exit(1)
