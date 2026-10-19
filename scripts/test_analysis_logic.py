import numpy as np

def calculate_transient_ratio(signal, split_idx):
    rms_early = np.sqrt(np.mean(signal[:split_idx]**2))
    rms_late = np.sqrt(np.mean(signal[split_idx:]**2))
    return rms_early / max(rms_late, 1e-12)

def test_transient_ratio():
    # Signal with strong early transient
    t = np.linspace(0, 600, 601) # 1 min spacing
    signal = np.exp(-t/50) * np.sin(2*np.pi*t/10) + 0.1 * np.sin(2*np.pi*t/20)
    ratio = calculate_transient_ratio(signal, 100)
    print(f"Transient ratio (strong): {ratio:.2f}")
    assert ratio > 2.0

    # Steady signal
    signal_steady = 0.1 * np.sin(2*np.pi*t/20) + 0.01 * np.random.randn(601)
    ratio_steady = calculate_transient_ratio(signal_steady, 100)
    print(f"Transient ratio (steady): {ratio_steady:.2f}")
    assert ratio_steady < 2.0

def test_trimming():
    t = np.arange(0, 601, 1.0) # 0, 1, ..., 600
    trim_min = 120.0
    # Logic in analyze_mars_mountain_sweep.py will be like t >= (t[0] + trim_s)
    trim_s = trim_min * 60.0
    times_s = t * 60.0
    keep = times_s >= (times_s[0] + trim_s)
    t_trimmed = times_s[keep]
    assert t_trimmed[0] == 120.0 * 60.0
    print(f"Trimmed t[0]: {t_trimmed[0]} s ({t_trimmed[0]/60} min)")

if __name__ == "__main__":
    test_transient_ratio()
    test_trimming()
    print("All tests passed!")
