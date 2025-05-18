import matplotlib.animation as animation
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq
from scipy.io import wavfile
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from tqdm import tqdm
import numpy as np
import sys
import os


if __name__ == "__main__":
    log = "log.txt"

    if not os.path.exists(log):
        print("Log file not found. Please run deformation to create it.")
        exit(1)

    if len(sys.argv) == 1:
        print("Missing argument.")
        exit(1)

    with open(log, "r") as f:
        magic = np.float64(f.readline().strip())
       
        row_ptr = np.array([int(x) for x in f.readline().strip().split()])
        col_idx = np.array([int(x) for x in f.readline().strip().split()])
        data = np.array([np.float64(x) for x in f.readline().strip().split()])
        Mcsr = csr_matrix((data, col_idx, row_ptr))
        Mcsr = (Mcsr + Mcsr.T) / 2.

        row_ptr = np.array([int(x) for x in f.readline().strip().split()])
        col_idx = np.array([int(x) for x in f.readline().strip().split()])
        data = np.array([np.float64(x) for x in f.readline().strip().split()])
        Kcsr = csr_matrix((data, col_idx, row_ptr))
        Kcsr = (Kcsr + Kcsr.T) / 2.

        lines = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

    print(f"Magic number: {magic:.15e}")

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

        t  = lines[begin:end, 0]
        Ep = lines[begin:end, 1]
        Ec = lines[begin:end, 2]
        E  = lines[begin:end, 3]

        t *= magic

        plt.plot(t, Ep, label="$E_p$", color="#43CCD0", linewidth=2, linestyle="--")
        plt.plot(t, Ec, label="$E_c$", color="#D043CC", linewidth=2, linestyle="--")
        plt.plot(t, E, label="$E$", color="#AAD043", linewidth=2)

        plt.xlabel("Time [$s$]", fontsize=12)
        plt.ylabel("Energy", fontsize=12)
        plt.title("Energy vs Time", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid()

        plt.show()
 
    elif analysis == "animation":
        if len(sys.argv) < 4:
            print("Missing dt and/or T parameters.")
            exit(1)

        try:
            dt = float(sys.argv[2])
            T = float(sys.argv[3])
        except:
            print("Expected floats.")
            exit(1)

        num_frames = int(T/dt)
        image_paths = [f"img/disp_{i}.png" for i in range(1, num_frames)]
        with imageio.get_writer('other_animation.gif', mode='I', duration=T) as writer:
            for i, filename in tqdm(enumerate(image_paths)):
                try:
                    image = imageio.imread(filename)
                    data = np.array(image)
                    if np.all(data == 0):
                        print(f"Warning: Image {filename} is empty. Skipping.")
                        continue
                    writer.append_data(image)
                except FileNotFoundError:
                    print(f"Warning: Image {filename} not found. Skipping.")

    elif analysis == "frequency":
        t = lines[:, 0]
        Ep = lines[:, 1]
        Ec = lines[:, 2]
        E = lines[:, 3]

        t *= magic

        dt = t[1] - t[0]
        N = len(t)

        # mean = np.mean(Ec)
        # print(f"Mean kinetic energy: {mean}")
        # spectrum = np.abs(fft(Ec-mean))[:N//2]
        # frequencies = fftfreq(N, dt)[:N//2]

        # peaks_e, _ = find_peaks(spectrum, height=1e-16, distance=10)
        # peaks_e = peaks_e[:2]
        # mid_freq = np.mean(frequencies[peaks_e])
        # diff_freq = frequencies[peaks_e][-1] - frequencies[peaks_e][0]

        # lines_xs = [frequencies[peaks_e][0]-300, frequencies[peaks_e][1]-300]
        # plt.plot(frequencies, spectrum, color="#D043CC", linewidth=2)
        # plt.hlines(spectrum[peaks_e], lines_xs, frequencies[peaks_e], color="#3E6BFF", linestyles="--")
        # for x, peak in zip(lines_xs, peaks_e):
        #     plt.text(x, spectrum[peak], f"{frequencies[peak]:.2f} [Hz]", fontsize=10, color="#3E6BFF", ha="left", va="bottom")

        # plt.plot(frequencies, spectrum, color="#43CCD0", linewidth=2)
        # plt.xlim(mid_freq-1.2*diff_freq, mid_freq+1.2*diff_freq)
        # plt.xlabel("Frequency ($Hz$)", fontsize=12)
        # plt.ylabel("Amplitude", fontsize=12)
        # plt.title("Spectrum of (Detrended) Kinetic Energy", fontsize=14)    
        # plt.grid()
        # plt.show()
        
        # if len(sys.argv) < 3:
        #     exit(0)

        if len(sys.argv) < 4:
            print("Missing the time file and/or node index.")
            exit(1)

        with open(sys.argv[2], "r") as f:
            moves = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

        t = moves[:, 0]
        x = moves[:, 1]
        y = moves[:, 2]
        vx = moves[:, 3]
        vy = moves[:, 4]

        t *= magic
        x *= magic
        y *= magic

        mean_x = np.mean(x)
        print(f"Mean x position of node {sys.argv[3]}: {mean_x} [m]")
        spectrum_x = np.abs(fft(x-mean_x))[:N//2]
        frequencies = fftfreq(N, dt)[:N//2]

        mean_y = np.mean(y)
        print(f"Mean y position of node {sys.argv[3]}: {mean_y} [m]")
        spectrum_y = np.abs(fft(y-mean_y))[:N//2]

        peaks_x, _ = find_peaks(spectrum_x, height=1e-5, distance=10)
        peaks_y, _ = find_peaks(spectrum_y, height=1e-5, distance=10)
        peaks = np.concatenate((peaks_x, peaks_y))
        peaks = np.unique(peaks)
        peaks = peaks[:3]
        mid_freq = np.mean(frequencies[peaks])
        diff_freq = frequencies[peaks][-1] - frequencies[peaks][0]

        eigs = eigsh(Kcsr, k=Mcsr.shape[0]-1, M=Mcsr, which="LM", return_eigenvectors=False)
        eigs_norm = np.abs(eigs)
        freqs = np.sqrt(eigs_norm)/(2.*np.pi)

        # Find closest eigenvalue to the first peak
        fundamental = peaks_x[0]
        closest_eigenvalue = np.argmin(np.abs(freqs - frequencies[fundamental]))
        print(f"First peak frequency: {frequencies[fundamental]} [Hz]")
        print(f"Closest eigenvalue to the first peak: {freqs[closest_eigenvalue]} [Hz]")

        lines_xs = [0, frequencies[peaks][0]*1.2, frequencies[peaks][1]*1.2]
        plt.plot(frequencies, spectrum_x, color="#D04394", linewidth=2, label="x")
        plt.plot(frequencies, spectrum_y, color="#43D07F", linewidth=2, label="y")
        plt.hlines(spectrum_y[peaks], lines_xs, frequencies[peaks], color="#3E6BFF", linestyles="--")
        for x, peak in zip(lines_xs, peaks):
            plt.text(x, spectrum_y[peak], f"{frequencies[peak]:.2f} [Hz]", fontsize=10, color="#3E6BFF", ha="left", va="bottom")
        
        plt.xlim(mid_freq-1.1*diff_freq, mid_freq+1.1*diff_freq)
        plt.xlabel("Frequency [$Hz$]", fontsize=12)
        plt.ylabel("Amplitude [$m$]", fontsize=12)
        plt.title(f"Spectrum of (Detrended) Position of Node {sys.argv[3]}", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid()
        plt.show()
    
    elif analysis == "state":
        if len(sys.argv) < 4:
            print("Missing the time file and/or node index.")
            exit(1)

        with open(sys.argv[2], "r") as f:
            moves = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

        try:
            node = int(sys.argv[3])
        except:
            print("Expected integer.")
            exit(1)

        t = moves[:, 0]
        x = moves[:, 1] 
        y = moves[:, 2]
        vx = moves[:, 3]
        vy = moves[:, 4]

        t *= magic
        x *= magic
        y *= magic

        # vx = np.convolve(vx, np.ones(20)/20, mode='same')
        # vy = np.convolve(vy, np.ones(20)/20, mode='same')

        factor = 10
        N = 20*1000#len(x)

        fig, (ax1) = plt.subplots(figsize=(5, 5))

        ax1.set_xlim(min(x) - np.std(x), max(x) + np.std(x))
        ax1.set_ylim(min(y) - np.std(y), max(y) + np.std(y))
        time_text = ax1.text(0.05, 0.95, '', transform=ax1.transAxes, fontsize=12, verticalalignment='top')
        ax1.set_title(f"Trajectory of Node {node}")
        line1, = ax1.plot([], [], 'bo', markersize=1, alpha=0.05)
        tracker1, = ax1.plot([], [], 'rx', markersize=7, alpha=1)

        # ax2.set_xlim(min(vx) - np.std(vx), max(vx) + np.std(vx))
        # ax2.set_ylim(min(vy) - np.std(vy), max(vy) + np.std(vy))
        # ax2.set_title(f"(Smoothed) Velocity of Node {node}")
        # line2, = ax2.plot([], [], 'ro', markersize=1, alpha=0.05)
        # tracker2, = ax2.plot([], [], 'rx', markersize=5, alpha=1)

        x_plot, y_plot = [], []
        # vx_plot, vy_plot = [], []

        def update(frame):
            x_plot.extend(x[factor*frame:factor*(frame+1)])
            y_plot.extend(y[factor*frame:factor*(frame+1)])

            # vx_plot.append(vx[factor*frame:factor*(frame+1)])
            # vy_plot.append(vy[factor*frame:factor*(frame+1)])

            time_text.set_text(f'Time = {1000*t[factor*(frame+1)-1]:.3f} ms')

            # if len(vx_plot) > 200:
            #     vx_plot.pop(0)
            #     vy_plot.pop(0)
            
            line1.set_data(x_plot, y_plot)
            tracker1.set_data([x_plot[-1]], [y_plot[-1]])

            # line2.set_data(vx_plot, vy_plot)
            # tracker2.set_data([vx_plot[-1]], [vy_plot[-1]])

            return line1, tracker1, time_text

        ani = animation.FuncAnimation(fig, update, frames=N//factor, interval=1, blit=True, repeat=False)

        plt.tight_layout()
        plt.show()

        query = input("Do you want to save the animation? This can take some time (y/n): ")
        if query.lower() == "y":
            ani.save(f"state_{node}.gif", writer='pillow', fps=30)
            print(f"Animation saved as state_{node}.gif")
        elif query.lower() == "n":
            print("Animation not saved.")
        
        plt.close(fig)

    elif analysis == "sound":
        if len(sys.argv) < 3:
            print("Missing the time file.")
            exit(1)

        with open(sys.argv[2], "r") as f:
            moves = np.array([list(map(np.float64, line.strip().split())) for line in f.readlines()])

        t = moves[:, 0]
        x = moves[:, 1] 
        y = moves[:, 2]

        t *= magic
        x *= magic
        y *= magic
        
        N = len(x)
        dt = t[1] - t[0]

        amplitudes = np.abs(fft(y))[:N//2]
        frequencies = fftfreq(N, dt)[:N//2]

        indices, _ = find_peaks(amplitudes, height=1e-6, distance=10)
        frequencies = frequencies[indices]
        amplitudes = amplitudes[indices]
        print(f"Frequencies: {frequencies}")

        # Normalize the amplitudes
        amplitudes /= np.max(amplitudes)
        sampling_rate = 44100
        duration = 3

        t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)
        y = 0.5 * np.sum([amp * np.sin(2 * np.pi * freq * t) for freq, amp in zip(frequencies, amplitudes)], axis=0)
        y_pcm = np.int16(y * 32767)

        plt.plot(t[:1000], y[:1000], color="#43CCD0", linewidth=2)
        plt.xlabel("Time [$s$]", fontsize=12)
        plt.ylabel("Amplitude", fontsize=12)
        plt.title("Generated Sound Wave", fontsize=14)
        plt.show()

        wavfile.write('tuning_fork.wav', sampling_rate, y_pcm)

    else:
        print("Unrecognized analysis option, either 'energy', 'animation', 'sound', 'state' or 'frequency'.")
        exit(1)
