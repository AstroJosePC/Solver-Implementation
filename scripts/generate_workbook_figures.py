"""Generate the seven textbook figures for docs/workbook Part I in strict black and white styling.

Figures generated:
1. fig/clocks.pdf       - Ch 3 Fig 3.1: Wave arrival and round-trip return clocks
2. fig/ridge_field.pdf   - Ch 3 Fig 3.2: 2D linear stationary vertical velocity over Gaussian ridges
3. fig/ridge_spectrum.pdf - Ch 3 Fig 3.3: Continuous flux density and discrete mode shares
4. fig/critical.pdf      - Ch 3 Fig 3.4: WKB critical level approach and parcel displacement
5. fig/layer_refl.pdf    - Ch 5 Fig 5.3: Sponge damping profile and reflected flux shares
6. fig/cavity_phase.pdf  - Ch 5 Fig 5.4: Upward flux vs reflection phase
7. fig/toy.pdf           - Ch 6 Fig 6.2: Single-mode Richtmyer toy model decay and flux

Styling matches the monochrome, high-contrast, framed appearance of the precompiled workbook figures.
"""
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Typography and global rcParams matching book style
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['DejaVu Serif', 'Palatino', 'Times New Roman']
matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = 0.8
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['xtick.color'] = 'black'
matplotlib.rcParams['ytick.color'] = 'black'

LINE_WIDTH = 1.2


def apply_bw_box_style(ax):
    """Apply consistent black-and-white boxed frame styling with no background grid."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.8)
    ax.tick_params(colors='black', labelsize=9)
    ax.grid(False)


# -----------------------------------------------------------------------------
# 1. clocks.pdf
# -----------------------------------------------------------------------------
def generate_clocks_fig(out_dir):
    """Generate Chapter 3 Figure 3.1: Arrival and return clocks (B&W)."""
    lx = 2000.0e3
    u = 40.0
    zb = 60.0e3
    z = np.linspace(0.0, 60.0e3, 601)

    n_surf = 0.0176
    n_free = 0.0085
    z_inv = 3.0e3
    n_z = n_free + (n_surf - n_free) * np.exp(-z / z_inv)

    dz = z[1] - z[0]
    iz = np.cumsum(n_z) * dz
    iz_b = iz[-1]

    fig, ax = plt.subplots(figsize=(6.5, 4.1))

    # 12 to 17 h scoring window shaded in subtle light gray
    ax.axvspan(12.0, 17.0, facecolor='0.92', edgecolor='0.7', linestyle=':', zorder=0, label='12 to 17 h window')

    modes = [1, 2, 3, 4, 5, 6]
    apex_times = {}

    for n in modes:
        kn = 2.0 * np.pi * n / lx
        t_up = iz / (kn * u**2) / 3600.0
        t_rt = (2.0 * iz_b - iz) / (kn * u**2) / 3600.0
        apex_times[n] = t_up[-1]

        ax.plot(t_up, z / 1.0e3, color='black', lw=LINE_WIDTH, linestyle='-')
        ax.plot(t_rt, z / 1.0e3, color='black', lw=LINE_WIDTH, linestyle='--')

    # Mode arrival apex labels at domain top z = 60 km
    # Modes 6, 4, 3, 2, 1 sit neatly above their respective peaks
    # Mode 5 is offset slightly higher with a thin connector to provide generous breathing room
    for n in [6, 4, 3, 2, 1]:
        t_top = apex_times[n]
        ax.text(t_top, 61.0, f"$n={n}$", fontsize=8.5, ha='center', va='bottom', clip_on=False)

    t_top_5 = apex_times[5]
    ax.annotate("$n=5$", xy=(t_top_5, 60.0), xytext=(t_top_5, 63.6),
                fontsize=8.5, ha='center', va='bottom', clip_on=False,
                arrowprops=dict(arrowstyle='-', color='black', lw=0.6, shrinkA=3, shrinkB=0))

    # Legend items
    ax.plot([], [], color='black', linestyle='-', lw=LINE_WIDTH, label=r'Upgoing front $t_{\rm up}(z)$')
    ax.plot([], [], color='black', linestyle='--', lw=LINE_WIDTH, label=r'Returned from 60 km $t_{\rm rt}(z)$')

    ax.set_xlim(0.0, 32.0)
    ax.set_ylim(0.0, 60.0)
    ax.set_xlabel('Time [h]', fontsize=9.5, color='black')
    ax.set_ylabel('Altitude $z$ [km]', fontsize=9.5, color='black')
    apply_bw_box_style(ax)
    ax.tick_params(top=False)

    # Legend placed in wide open quadrant (t in [24, 32], z in [24, 40])
    ax.legend(loc='center right', bbox_to_anchor=(0.98, 0.48),
              frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.5)
    fig.tight_layout()

    out_path = os.path.join(out_dir, 'clocks.pdf')
    fig.savefig(out_path, format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 2. ridge_field.pdf
# -----------------------------------------------------------------------------
def generate_ridge_field_fig(out_dir):
    """Generate Chapter 3 Figure 3.2: 2D stationary vertical velocity over ridges (B&W)."""
    u = 40.0
    n_val = 0.0085
    c_val = 230.0
    h_scale = 7000.0
    h_peak = 100.0

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)

    cases = [
        ('Hydrostatic (sigma = 50 km)', 50.0e3, 200.0e3, 60.0e3, axes[0]),
        ('Dispersive (sigma = 5 km)', 5.0e3, 30.0e3, 40.0e3, axes[1])
    ]

    for title, sigma, x_span, z_max, ax in cases:
        nx = 256
        nz = 160
        x = np.linspace(-x_span, x_span, nx)
        z = np.linspace(0.0, z_max, nz)
        dx = x[1] - x[0]

        k = 2.0 * np.pi * np.fft.fftfreq(nx, dx)
        h_tilde = h_peak * sigma * np.sqrt(2.0 * np.pi) * np.exp(-0.5 * k**2 * sigma**2)

        m2 = (n_val / u)**2 - k**2 + (k * u / c_val)**2 - 1.0 / (4.0 * h_scale**2)
        m = np.zeros(nx, dtype=complex)
        prop = m2 > 0
        evan = m2 <= 0
        m[prop] = np.sign(k[prop]) * np.sqrt(m2[prop])
        m[evan] = 1j * np.sqrt(-m2[evan])

        w_field = np.zeros((nz, nx))
        for iz, z_val in enumerate(z):
            w_hat = -1j * k * u * h_tilde * np.exp(1j * m * z_val)
            w_field[iz, :] = np.real(np.fft.ifft(w_hat))

        norm_w = w_field / (np.max(np.abs(w_field)) + 1.0e-9)

        x_km = x / 1.0e3
        z_km = z / 1.0e3
        levels_pos = np.linspace(0.15, 0.95, 6)
        levels_neg = np.sort(-np.linspace(0.15, 0.95, 6))

        ax.contour(x_km, z_km, norm_w, levels=levels_pos, colors='black', linestyles='-', linewidths=0.9)
        ax.contour(x_km, z_km, norm_w, levels=levels_neg, colors='black', linestyles='--', linewidths=0.9)

        # Draw terrain shape at bottom in solid grayscale
        h_shape = h_peak * np.exp(-0.5 * (x / sigma)**2)
        h_scaled = (h_shape / h_peak) * (0.07 * (z_max / 1.0e3))
        ax.fill_between(x_km, 0, h_scaled, facecolor='0.75', edgecolor='black', lw=0.9, zorder=5)

        ax.set_xlim(-x_span / 1.0e3, x_span / 1.0e3)
        ax.set_ylim(0.0, z_max / 1.0e3)
        ax.set_xlabel('Distance $x$ [km]', fontsize=9.5, color='black')
        ax.set_title(title, fontsize=10, color='black')
        apply_bw_box_style(ax)

    axes[0].set_ylabel('Altitude $z$ [km]', fontsize=9.5, color='black')
    fig.tight_layout()

    out_path = os.path.join(out_dir, 'ridge_field.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 3. ridge_spectrum.pdf
# -----------------------------------------------------------------------------
def generate_ridge_spectrum_fig(out_dir):
    """Generate Chapter 3 Figure 3.3: Ridge spectrum and discrete mode shares (B&W)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Panel (a): Continuous dimensionless flux density
    ksigma = np.linspace(0.0, 2.5, 300)
    flux_density = ksigma * np.exp(-ksigma**2)
    flux_density_norm = flux_density / np.max(flux_density)

    ax1.plot(ksigma, flux_density_norm, color='black', lw=LINE_WIDTH)
    peak_idx = 1.0 / np.sqrt(2.0)
    ax1.axvline(peak_idx, color='0.5', linestyle=':', lw=0.9)
    ax1.plot([peak_idx], [1.0], marker='o', markerfacecolor='black', markeredgecolor='black', markersize=5.0)
    ax1.text(peak_idx + 0.08, 0.95, r'$k\sigma = 1/\sqrt{2}$', fontsize=8.5, color='black')

    ax1.set_xlim(0.0, 2.5)
    ax1.set_ylim(0.0, 1.08)
    ax1.set_xlabel(r'Dimensionless wavenumber $k\sigma$', fontsize=9.5, color='black')
    ax1.set_ylabel(r'Flux density $k|\tilde{h}|^2$ [normalized]', fontsize=9.5, color='black')
    ax1.set_title('(a) Continuous Gaussian spectrum', fontsize=10, color='black')
    apply_bw_box_style(ax1)

    # Panel (b): Discrete mode shares for Lx = 2000 km
    lx = 2000.0e3
    u = 40.0
    n_val = 0.0085
    z_target = 55.0e3
    modes = np.arange(1, 13)

    width = 0.36
    # Contrast: light gray fill vs striped hatching
    k_modes = 2.0 * np.pi * modes / lx

    raw_50 = k_modes * np.exp(-k_modes**2 * (50.0e3)**2)
    share_50 = raw_50 / np.sum(raw_50)

    raw_100 = k_modes * np.exp(-k_modes**2 * (100.0e3)**2)
    share_100 = raw_100 / np.sum(raw_100)

    ax2.bar(modes - width/2, share_50, width=width, label=r'$\sigma = 50\,$km',
            facecolor='0.85', edgecolor='black', lw=0.7)
    ax2.bar(modes + width/2, share_100, width=width, label=r'$\sigma = 100\,$km',
            facecolor='none', edgecolor='black', hatch='///', lw=0.7)

    # Arrival times above mode bars
    for n in range(1, 7):
        kn = 2.0 * np.pi * n / lx
        t_arr = (n_val * z_target) / (kn * u**2) / 3600.0
        ax2.text(n - width/2, 0.28, f"{t_arr:.1f}h", ha='center', va='bottom', fontsize=7.0, rotation=90, color='black')

    ax2.set_xlim(0.3, 12.7)
    ax2.set_ylim(0.0, 0.38)
    ax2.set_xticks(modes)
    ax2.set_xlabel('Horizontal mode $n$', fontsize=9.5, color='black')
    ax2.set_ylabel('Fraction of hydrostatic flux', fontsize=9.5, color='black')
    ax2.set_title('(b) Discrete mode shares and arrival at 55 km', fontsize=10, color='black')
    ax2.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.5, loc='upper right')
    apply_bw_box_style(ax2)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'ridge_spectrum.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 4. critical.pdf
# -----------------------------------------------------------------------------
def generate_critical_fig(out_dir):
    """Generate Chapter 3 Figure 3.4: Critical level WKB scaling and displacement (B&W)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Panel (a): WKB scaling near critical level
    z_zc = np.linspace(0.0, 0.95, 200)
    u_rel = 1.0 - z_zc
    w_rel = np.sqrt(u_rel)
    u_pert_rel = 1.0 / np.sqrt(u_rel)
    lam_rel = u_rel

    ax1.plot(w_rel, z_zc, color='black', linestyle='-', lw=LINE_WIDTH, label=r'$|w| \propto U^{1/2}$')
    ax1.plot(np.clip(u_pert_rel, 0, 4.0), z_zc, color='black', linestyle='--', lw=LINE_WIDTH, label=r"$|u'| \propto U^{-1/2}$")
    ax1.plot(lam_rel, z_zc, color='black', linestyle=':', lw=LINE_WIDTH, label=r'$\lambda_z \propto U$')
    ax1.axhline(1.0, color='0.5', linestyle=':', lw=0.9)
    ax1.text(0.1, 0.98, 'Critical level $z_c$', color='black', fontsize=8.5)

    ax1.set_xlim(0.0, 3.5)
    ax1.set_ylim(0.0, 1.05)
    ax1.set_xlabel('Normalized magnitude', fontsize=9.5, color='black')
    ax1.set_ylabel('Normalized height $z/z_c$', fontsize=9.5, color='black')
    ax1.set_title('(a) WKB scaling near critical level', fontsize=10, color='black')
    ax1.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.0, loc='center right')
    apply_bw_box_style(ax1)

    # Panel (b): Parcel displacement
    t_hours = np.linspace(0.0, 3.0, 300)
    t_sec = t_hours * 3600.0
    w_imposed = 0.024
    k_mode = 1.0e-5

    u_fast = 40.0
    zeta_fast = (w_imposed / (k_mode * u_fast)) * np.sin(k_mode * u_fast * t_sec)
    zeta_secular = w_imposed * t_sec

    ax2.plot(t_hours, zeta_secular, color='black', linestyle='-', lw=1.4, label=r'Calm air ($U \leq 1\,$m/s): $\zeta = w t$')
    ax2.plot(t_hours, zeta_fast, color='black', linestyle='--', lw=1.2, label=r'Wind $U = 40\,$m/s: bounded oscillation')

    ax2.set_xlim(0.0, 3.0)
    ax2.set_ylim(-150.0, 300.0)
    ax2.set_xlabel('Time [h]', fontsize=9.5, color='black')
    ax2.set_ylabel('Vertical displacement $\zeta$ [m]', fontsize=9.5, color='black')
    ax2.set_title('(b) Displacement under imposed $w$', fontsize=10, color='black')
    ax2.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.0, loc='upper left')
    apply_bw_box_style(ax2)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'critical.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 5. layer_refl.pdf
# -----------------------------------------------------------------------------
def generate_layer_refl_fig(out_dir):
    """Generate Chapter 5 Figure 5.3: Sponge damping profile and reflected share (B&W)."""
    sys.path.insert(0, os.path.abspath('scripts'))
    from analyze_linear_flux_constancy import rayleigh_sponge_reflection

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Panel (a): Damping profile
    z_sponge = np.linspace(60.0, 95.0, 200)
    z_base = 60.0
    z_top = 95.0
    k0_val = 0.05
    k_profile = k0_val * np.sin(0.5 * np.pi * (z_sponge - z_base) / (z_top - z_base))**2

    ax1.plot(k_profile, z_sponge, color='black', lw=1.4, label=r'$K(z)$ ($K_0 = 0.05\,$s$^{-1}$)')

    lx = 2000.0e3
    u = 40.0
    modes = [1, 3, 6]
    linestyles = [':', '--', '-.']

    for n, ls in zip(modes, linestyles):
        kn = 2.0 * np.pi * n / lx
        omega_hat = kn * u
        ax1.axvline(omega_hat, color='black', linestyle=ls, lw=1.1, label=rf'$|\hat{{\omega}}|$ mode {n}')

    ax1.set_xlim(0.0, 0.052)
    ax1.set_ylim(60.0, 95.0)
    ax1.set_xlabel(r'Rate [s$^{-1}$]', fontsize=9.5, color='black')
    ax1.set_ylabel('Altitude $z$ [km]', fontsize=9.5, color='black')
    ax1.set_title('(a) Damping rate vs intrinsic frequency', fontsize=10, color='black')
    ax1.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.0, loc='lower right')
    apply_bw_box_style(ax1)

    # Panel (b): Reflected share R at 40 km
    z_grid = np.linspace(0.0, 95.0e3, 951)
    n_const = np.full_like(z_grid, 0.0085)
    h_scale = 7000.0
    rho_const = np.exp(-z_grid / h_scale)

    n_eval = [1, 2, 3, 4, 6, 8, 10, 12]
    k_eval = [2.0 * np.pi * n / lx for n in n_eval]
    idx_40 = int(np.searchsorted(z_grid, 40.0e3))

    r_shares_k05 = []
    r_shares_k005 = []
    for k_val in k_eval:
        r_prof_05 = rayleigh_sponge_reflection(k_val, z_grid, n_const, rho_const, u, 0.05,
                                               z_base * 1.0e3, z_top * 1.0e3, top="radiation")
        r_prof_005 = rayleigh_sponge_reflection(k_val, z_grid, n_const, rho_const, u, 0.005,
                                                z_base * 1.0e3, z_top * 1.0e3, top="radiation")
        r_shares_k05.append(r_prof_05[idx_40])
        r_shares_k005.append(r_prof_005[idx_40])

    ax2.plot(n_eval, r_shares_k05, marker='o', markerfacecolor='black', markeredgecolor='black',
             markersize=4.5, color='black', lw=1.2, linestyle='-', label=r'$K_0 = 0.05\,$s$^{-1}$')
    ax2.plot(n_eval, r_shares_k005, marker='s', markerfacecolor='white', markeredgecolor='black',
             markersize=4.5, color='black', lw=1.2, linestyle='--', label=r'$K_0 = 0.005\,$s$^{-1}$')

    ax2.set_yscale('log')
    ax2.set_ylim(1.0e-4, 1.0)
    ax2.set_xlabel('Horizontal mode index $n$', fontsize=9.5, color='black')
    ax2.set_ylabel(r'Reflected flux share $|r|^2$ at 40 km', fontsize=9.5, color='black')
    ax2.set_title('(b) Downward-integrated reflection', fontsize=10, color='black')
    ax2.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.5, loc='upper right')
    apply_bw_box_style(ax2)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'layer_refl.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 6. cavity_phase.pdf
# -----------------------------------------------------------------------------
def generate_cavity_phase_fig(out_dir):
    """Generate Chapter 5 Figure 5.4: Cavity upward flux vs reflection phase (B&W)."""
    phi = np.linspace(-np.pi, np.pi, 300)
    phi_deg = phi * 180.0 / np.pi

    fig, ax = plt.subplots(figsize=(5.5, 3.2))

    rg_cases = [
        (0.05, ':', 'black', 1.0, r'$|r_g| = 0.05$ (0.25% flux share)'),
        (0.10, '-', 'black', 1.5, r'$|r_g| = 0.10$ (1.0% flux share)'),
        (0.15, '--', 'black', 1.1, r'$|r_g| = 0.15$ (2.25% flux share)'),
    ]

    for rg, ls, col, lw, label in rg_cases:
        ratio = 1.0 / (1.0 + rg**2 + 2.0 * rg * np.cos(phi))
        ax.plot(phi_deg, ratio, linestyle=ls, color=col, lw=lw, label=label)

    ax.axhline(1.0, color='0.6', linestyle='-', lw=0.7)

    # Pure black annotations
    ax.annotate(r'$\approx 1.23$ (constructive)', xy=(180, 1.23), xytext=(85, 1.31),
                arrowprops=dict(arrowstyle='->', color='black', lw=0.8), fontsize=8.0, color='black')
    ax.annotate(r'$\approx 0.83$ (destructive)', xy=(0, 0.83), xytext=(-85, 0.72),
                arrowprops=dict(arrowstyle='->', color='black', lw=0.8), fontsize=8.0, color='black')

    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(0.65, 1.45)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
    ax.set_xlabel(r'Reflection phase $\phi = \arg(r_g)$ [rad]', fontsize=9.5, color='black')
    ax.set_ylabel(r'Upward flux ratio $F_\uparrow / F_0$', fontsize=9.5, color='black')
    apply_bw_box_style(ax)
    ax.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.0, loc='upper left')
    fig.tight_layout()

    out_path = os.path.join(out_dir, 'cavity_phase.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# 7. toy.pdf
# -----------------------------------------------------------------------------
def solve_single_mode_toy(dt, path_type='face', t_end=36000.0):
    """Solve the 1D single horizontal Fourier mode acoustic-gravity system.

    Two-step Richtmyer Lax-Wendroff vertical scheme with buoyancy source.
    """
    u_wind = 40.0
    n_val = 0.01
    c_val = 230.0
    lx = 2000.0e3
    k = 2.0 * np.pi * 6.0 / lx

    dz = 500.0
    z_max = 60.0e3
    nz = int(z_max / dz)
    z = (np.arange(nz) + 0.5) * dz

    zb = 35.0e3
    k_rate = np.where(z >= zb, 0.05 * np.sin(0.5 * np.pi * np.clip((z - zb) / (z_max - zb), 0, 1))**2, 0.0)

    u = np.zeros(nz, dtype=complex)
    w = np.zeros(nz, dtype=complex)
    b = np.zeros(nz, dtype=complex)
    p = np.zeros(nz, dtype=complex)

    nt = int(t_end / dt)
    t0 = 7200.0
    w_amp = 0.01

    rot_full = np.exp(-1j * k * u_wind * dt)

    def flux_G(w_arr, p_arr):
        return p_arr, c_val**2 * w_arr

    for step in range(nt):
        t_n = step * dt
        wf_n = w_amp * np.tanh(t_n / t0)
        wf_half = w_amp * np.tanh((t_n + 0.5 * dt) / t0)

        # Sources at cell centres
        Su = -1j * k * p - k_rate * u
        Sw = b - k_rate * w
        Sb = -n_val**2 * w - k_rate * b
        Sp = -c_val**2 * (1j * k * u) - k_rate * p

        # Fluxes at cell centres
        Gw, Gp = flux_G(w, p)

        # Predictor at interior faces
        w_half = 0.5 * (w[:-1] + w[1:]) - (dt / (2.0 * dz)) * (Gw[1:] - Gw[:-1]) + 0.25 * dt * (Sw[:-1] + Sw[1:])
        p_half = 0.5 * (p[:-1] + p[1:]) - (dt / (2.0 * dz)) * (Gp[1:] - Gp[:-1]) + 0.25 * dt * (Sp[:-1] + Sp[1:])

        if path_type == 'ghost':
            w_ghost = 2.0 * wf_n - w[0]
            p_ghost = p[0]
            Gw_g, Gp_g = flux_G(w_ghost, p_ghost)
            Sw_g = b[0] - k_rate[0] * w_ghost
            Sp_g = -c_val**2 * (1j * k * u[0]) - k_rate[0] * p_ghost

            w_face_0 = 0.5 * (w_ghost + w[0]) - (dt / (2.0 * dz)) * (Gw[0] - Gw_g) + 0.25 * dt * (Sw_g + Sw[0])
            p_face_0 = 0.5 * (p_ghost + p[0]) - (dt / (2.0 * dz)) * (Gp[0] - Gp_g) + 0.25 * dt * (Sp_g + Sp[0])
        else:
            w_face_0 = wf_half
            p_face_0 = p[0]

        w_face_top = w[-1]
        p_face_top = p[-1]

        w_faces = np.concatenate(([w_face_0], w_half, [w_face_top]))
        p_faces = np.concatenate(([p_face_0], p_half, [p_face_top]))

        Gw_faces, Gp_faces = flux_G(w_faces, p_faces)

        # Corrector update
        dGw = (Gw_faces[1:] - Gw_faces[:-1]) / dz
        dGp = (Gp_faces[1:] - Gp_faces[:-1]) / dz

        w_new = w - dt * dGw + dt * Sw
        p_new = p - dt * dGp + dt * Sp
        u_new = u + dt * Su
        b_new = b + dt * Sb

        u = u_new * rot_full
        w = w_new * rot_full
        b = b_new * rot_full
        p = p_new * rot_full

    return z, w, u, p


def generate_toy_fig(out_dir):
    """Generate Chapter 6 Figure 6.2: Single-mode toy model results (B&W)."""
    t_start = time.time()

    dz = 500.0
    dt1 = 0.74
    dt2 = 0.37

    print("Running toy model simulations...")
    z, w_g1, u_g1, _ = solve_single_mode_toy(dt1, 'ghost', t_end=36000.0)
    z, w_f1, u_f1, _ = solve_single_mode_toy(dt1, 'face', t_end=36000.0)
    z, w_g2, u_g2, _ = solve_single_mode_toy(dt2, 'ghost', t_end=36000.0)
    z, w_f2, u_f2, _ = solve_single_mode_toy(dt2, 'face', t_end=36000.0)

    wall_clock = time.time() - t_start
    print(f"Toy model simulations completed in {wall_clock:.2f} s wall-clock time.")

    def filter_2dz(f):
        amp = np.zeros(len(f))
        for j in range(1, len(f) - 1):
            amp[j] = np.abs(0.25 * f[j-1] - 0.5 * f[j] + 0.25 * f[j+1])
        amp[0] = amp[1]
        amp[-1] = amp[-2]
        return amp

    amp_g1 = filter_2dz(np.abs(w_g1))
    amp_f1 = filter_2dz(np.abs(w_f1))
    amp_g2 = filter_2dz(np.abs(w_g2))
    amp_f2 = filter_2dz(np.abs(w_f2))

    norm_g1 = amp_g1 / (amp_g1[0] + 1.0e-12)
    norm_g2 = amp_g2 / (amp_g2[0] + 1.0e-12)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    # Panel (a): Decay of 2dz mode near ground (B&W)
    z_km = z / 1.0e3
    idx_plot = z_km <= 7.0

    ax1.semilogy(z_km[idx_plot], norm_g1[idx_plot], 'o-', color='black',
                 markerfacecolor='black', markeredgecolor='black', markersize=4, lw=1.0,
                 label=r'Ghost ($\Delta t = 0.74\,$s, $L \approx 0.71\,$km)')
    ax1.semilogy(z_km[idx_plot], norm_g2[idx_plot], 's--', color='black',
                 markerfacecolor='white', markeredgecolor='black', markersize=4.5, lw=1.0,
                 label=r'Ghost ($\Delta t = 0.37\,$s, $L \approx 1.45\,$km)')

    l1 = 0.705
    l2 = 1.455
    z_ref = np.linspace(0.25, 6.0, 100)
    ax1.plot(z_ref, np.exp(-z_ref / l1), ':', color='0.4', lw=1.2, label='Analytic slope $-1/L_1$')
    ax1.plot(z_ref, np.exp(-z_ref / l2), '-.', color='0.4', lw=1.2, label='Analytic slope $-1/L_2$')

    ax1.set_xlim(0.0, 6.5)
    ax1.set_ylim(1.0e-4, 1.5)
    ax1.set_xlabel('Altitude $z$ [km]', fontsize=9.5, color='black')
    ax1.set_ylabel(r'Relative $2\Delta z$ amplitude of $w$', fontsize=9.5, color='black')
    ax1.set_title('(a) Grid-scale mode vertical decay', fontsize=10, color='black')
    ax1.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=7.5, loc='lower left')
    apply_bw_box_style(ax1)

    # Panel (b): Band-mean momentum flux (B&W)
    band = (z >= 5.0e3) & (z <= 25.0e3)
    def momentum_flux(u_arr, w_arr):
        return 0.5 * np.real(u_arr * np.conj(w_arr))

    mf_f2_ref = np.mean(momentum_flux(u_f2[band], w_f2[band]))
    rel_f1 = np.mean(momentum_flux(u_f1[band], w_f1[band])) / mf_f2_ref
    rel_f2 = 1.0
    rel_g1 = np.mean(momentum_flux(u_g1[band], w_g1[band])) / mf_f2_ref
    rel_g2 = np.mean(momentum_flux(u_g2[band], w_g2[band])) / mf_f2_ref

    dt_vals = [dt2, dt1]
    face_vals = [rel_f2, rel_f1]
    ghost_vals = [rel_g2, rel_g1]

    ax2.plot(dt_vals, face_vals, 's-', color='black', markerfacecolor='black', markeredgecolor='black',
             markersize=5, lw=LINE_WIDTH, label='Face path (invariant)')
    ax2.plot(dt_vals, ghost_vals, 'o--', color='black', markerfacecolor='white', markeredgecolor='black',
             markersize=5.5, lw=LINE_WIDTH, label=r'Ghost path ($\Delta t$ dependent)')
    ax2.axhline(1.0, color='0.6', linestyle=':', lw=0.8)

    ax2.set_xlim(0.2, 0.9)
    ax2.set_ylim(0.97, 1.02)
    ax2.set_xlabel(r'Time step $\Delta t$ [s]', fontsize=9.5, color='black')
    ax2.set_ylabel('Band flux relative to face path', fontsize=9.5, color='black')
    ax2.set_title('(b) 5 to 25 km band-mean flux', fontsize=10, color='black')
    ax2.legend(frameon=True, facecolor='white', edgecolor='0.8', fontsize=8.0, loc='lower left')
    apply_bw_box_style(ax2)

    fig.tight_layout()
    out_path = os.path.join(out_dir, 'toy.pdf')
    fig.savefig(out_path, format='pdf')
    plt.close(fig)
    print(f"Generated {out_path}")


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
def main():
    import sys
    out_dir = os.path.join('docs', 'workbook', 'fig')
    os.makedirs(out_dir, exist_ok=True)

    target = sys.argv[1] if len(sys.argv) > 1 else 'all'

    generators = {
        'clocks': generate_clocks_fig,
        'ridge_field': generate_ridge_field_fig,
        'ridge_spectrum': generate_ridge_spectrum_fig,
        'critical': generate_critical_fig,
        'layer_refl': generate_layer_refl_fig,
        'cavity_phase': generate_cavity_phase_fig,
        'toy': generate_toy_fig,
    }

    if target in generators:
        print(f"Generating {target} figure...")
        generators[target](out_dir)
        print(f"{target} figure generated successfully.")
    else:
        print("Generating all black-and-white workbook figures...")
        for name, gen in generators.items():
            gen(out_dir)
        print("All black-and-white workbook figures generated successfully.")


if __name__ == '__main__':
    main()
