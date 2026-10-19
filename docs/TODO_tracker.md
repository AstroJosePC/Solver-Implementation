# MAGNUS-P Solver: TODO Tracker & Decisions Log

## Outstanding TODOs

### HIGH PRIORITY
- [ ] **Obtain Eckermann et al. (2011) CO2 radiative damping coefficients**
  - Paper: "Scale-dependent infrared radiative damping rates on Mars..." (Icarus)
  - Need: cubic polynomial fit coefficients a(z), b(z), c(z), d(z) for tau_r^{-1}(z,m)
  - These are altitude- and wavenumber-dependent; essential for realistic Mars AGW damping
  - Currently using constant Newtonian cooling as placeholder (Eq. 25-26)
  - Location in code: `src/damping.jl`

- [ ] **Resolve residual sponge reflection discrepancy (REDUCED — 60 km sponge gives r=3.74%)**
  - With 20 km sponge: r=19.7% (K0=0.005). With 40 km: r=7.4%. With 60 km: r=3.74%.
  - Paper Table 4 target: K0=0.005 → r=0.55%. Still ~7x higher at 60 km.
  - Klemp et al. (2008) rule (cited by paper): sponge must be ≥1 vertical wavelength thick.
    Our lambda_z=46.1 km (T=239K, N=0.020), so 20 km sponge (0.43 lambda_z) violates rule.
    60 km sponge = 1.30 lambda_z is above the minimum; further improvement possible with domain
    extension (e.g., 200 km domain with 80 km sponge = 1.74 lambda_z).
  - Paper's N/T inconsistency is a likely contributor: if paper used N=0.0182 (implying T≈289K),
    their effective lambda_z for sponge test would differ from ours.
  - BC approaches explored and rejected: Sommerfeld gradient BC (failed); Thompson/Poinsot-Lele
    characteristic BC (failed). Both crash at step ~14-400 due to near-vacuum rho~1e-8 at
    z=120 km amplifying floating-point momentum noise into enormous velocities. Literature
    survey confirms no AGW simulation paper uses characteristic/radiation BCs — all use sponge.
  - Eq. 19 outflow BC confirmed correct and stable. It is intentionally imperfect; sponge does
    the work (paper explicitly states this in Section 2.2.6).
  - Next option: extend domain to 200 km with 80 km sponge and re-run.

- [ ] **Resolve paper N inconsistency for dispersion test**
  - Paper Section 3.1 states T=239K (twice, verbatim) but quotes N=0.0182 rad/s
  - T=239K + paper's own formula N=sqrt(gamma-1)*g/cs → N=0.02002 (not 0.0182)
  - N=0.0182 implies cs=341 m/s → T≈289K (not 239K)
  - Consequence: our lambda_z=14.93 km (from N=0.020) vs paper's 16.7 km (from N=0.0182)
  - Our simulation validated internally (measured 15.02 km vs predicted 14.93 km, 0.6% error)
  - Paper comparison blocked by paper's own parameter inconsistency

- [ ] **Decide how to treat the breaking regime found in the 6 h U=40/h=3 extension**
  - The extension does not reach a steady state; it enters overturning with a decelerated
    layer descending from ~55 km through the AMEC band (section "Arsia Mons Stage 1
    U=40/h=3 time extension"). Long Stage 1 integrations of the h=3 km modes therefore
    have no quasi-steady endpoint to report against.
  - Three separable questions, roughly in order of cost:
    (a) Is the late-time deceleration maximum influenced by the sponge? Partly bounded as
        of 2026-08-31 (`docs/ext6h_instability_mechanism.md`): the sponge-adjacent 56-60 km
        band accounts for 0.2-4.8% of the deposition while 58-67% occurs over 25-40 km, but
        the 57.5-60 km band carries 3-5x more static instability per unit resolved wave slope
        than 45-57 km, and the first overturning altitude sits at the top of the sponge-free
        column. Narrowed 2026-09-01 by the 12 h record: that excess falls from a factor of
        4.6 at 6 h to 1.06 at 12 h as the interior instability grows to match it, so it
        appears to be a startup-phase feature rather than a persistent one. A repeat at the
        documented `K0 = 0.005` and with the sponge base raised (taller domain) would close
        it; ~2400 s each at 6 h based on the measured 6253 s for 12 h.
    (b) Does radiative damping bound the growth? The extension runs with
        `enable_cooling = false`; the Eckermann et al. (2011) coefficients above are the
        blocking item.
    (c) Is the breaking resolved? At dx = 4 km, dz = 500 m it is not. Quantified 2026-08-31:
        the vertical wavelength stays resolved by 19-54 points, but the horizontal scale
        needed to carry the measured flux at saturation amplitude falls from ~980 km at 3 h
        to 45 km at 6 h at z = 55 km, which is 11 points at dx = 4 km, so the post-5 h state
        is marginally resolved. Revised 2026-09-01: 6 h is the worst-resolved point of the
        12 h record, not the endpoint. As the wind at 55 km recovers from 14.5 to 26.0 m/s
        between 6 h and 12 h the same measure improves from 12.8 to 98.6 points, so the
        marginal state is transient. A dx = 2 km repeat over 6 h would still show whether the
        indicators at that worst point are grid-sensitive. The vertical side is supported as
        of 2026-09-01 by the linear control (`docs/linear_flux_constancy.md`), which recovers
        lambda_z from the measured phase tilt to a band mean of 0.996 of 2 pi U / N(z) over
        10-55 km at 53-64 points per wavelength; the horizontal side is untested at h = 0.1 km
        because no saturation occurs there.
    (d) Does the layer base stall near 26-34 km as the saturation-level calculation predicts?
        Answered 2026-09-01 by the 12 h extension, partially. The base occupies the predicted
        band from 6.5 h (33.25 km) to 10.5 h (26.25 km) with a descent rate that reverses
        sign in four of nine samples, then falls to 11.25 km by 12 h. The premise behind the
        25.8 km asymptote is not borne out: the flux at 25 km plateaus at 24-29% of the
        steady linear value rather than approaching it. Superseded by (e).
    (e) Is anything after ~7 h interpretable in Stage 1? Two independent diagnostics say the
        domain-mean becomes unreliable there. The x-mean budget closes to within 20% from
        2-6 h then degrades to a ratio of 0.25-0.70 over 7-12 h, and the advective transit
        from the inflow edge to the ridge is LX/(2U) = 6.94 h. Separately the domain spans
        only 2.87 flux-weighted wavelengths, and the ridge-window minus far-field mean wind
        reaches 26.2 m/s at 49.75 km at 12 h, alternating sign with height. Stage 2 is
        periodic in x and closes the budget without a lateral term, so extending Stage 2 is
        now the highest-value next run rather than a further Stage 1 extension.
    (f) What is the solver's own error on MF(z)? Opened 2026-09-01 by the linear control
        (uniform wind, periodic-x, h = 0.1 km, 12 h, `docs/linear_flux_constancy.md`). Where
        the wave is conservative Eliassen-Palm requires dMF/dz = 0, and molecular viscosity
        accounts for only 0.03% of the flux by 55 km, so any drift is numerical. Measured
        over 10-55 km in the final hour: 50.6% peak-to-peak with a fitted drift of -1.24% per
        km, so the flux is not flat at the 5% and 0.1%/km reporting thresholds. That is an
        upper bound rather than the error bar, because every measure is still decreasing in
        time (peak-to-peak 6.82 at 1 h to 0.476 at 12 h) and the spin-up is height-dependent:
        d ln|MF|/dt over the final 3 h is -1.5%/h at 4.75 km and +1.2%/h at 9.75 km, against
        +11.1%/h at 39.75 km and +20.0%/h at 54.75 km. Arrival of the wave front, which the
        group velocity places at 55 km by 9.33 h, is not the same as settling of the flux.
        Narrowed 2026-09-01 without a further run, by benchmarking against the arrived part
        of the ridge spectrum rather than a fully established wave (section 8 of
        `scripts/analyze_linear_flux_constancy.py`). Arrival time scales as 1/k, and in the
        2000 km periodic domain five discrete modes carry 95% of the flux; at 55 km after
        12 h the two longest (lambda_x = 2000 and 1000 km, together 45.2% of the flux) have
        not arrived, reaching that height only at 26.8 h and 13.4 h. Dividing by the arrived
        fraction, snapshot by snapshot, drops the fitted drift from -1.2426 to +0.0372% per
        km, and the band-fit spurious drag from -0.008 / -0.029 / -0.114 / -0.466 m/s/sol at
        9.75 / 24.75 / 39.75 / 54.75 km to +0.0002 / +0.001 / +0.004 / +0.026. The front
        factor that would null the drift is f = 0.989 against the parameter-free f = 1, so
        the WKB group velocity places the fronts where the record puts them without tuning.
        The height dependence of MF in this control is therefore consistent with an arrival
        transient. What remains open is the size of the numerical component, which the
        correction does not resolve below about 1%/km, the spread over the plausible range
        of the front-sharpness assumption (-0.88%/km at f = 0.5 to +3.30 at f = 2).
        Narrowed again 2026-09-02 with a sigma = 50 km twin (2587.66 s wall-clock, 58609
        steps, converged, timestep flat at 0.737-0.740 s). Correction to the earlier
        framing: narrowing the ridge does not move any mode's arrival time, since
        t_arr = I(z)/(k U^2) contains k and not sigma; it reweights the flux from the slow
        modes to the fast ones, raising the arrived fraction at 55 km at 12 h from 0.548 to
        0.862. Because m = N/U and the steady amplitude N h (rho_0 N_0/rho N)^0.5 both carry
        no k, lambda_z, the sponge thickness in wavelengths, dz per wavelength and the
        linearity margin are identical between the two runs, which makes the pair a
        controlled difference. Confirmed in the twin: lambda_z band mean measured/predicted
        1.0200 against 0.996, min_x u = 0.543 U with no overturning, max |w| below the
        sponge 8.53 m/s. Model-free outcome: reweighting the spectrum alone moved the
        measured drift over 10-55 km from -1.2426 to -0.1138%/km (a factor of 10.9),
        peak-to-peak from 50.6 to 8.66%, accumulated departure from -55.9 to -5.12%, and the
        band-fit spurious drag at 54.75 km from -0.466 to -0.080 m/s/sol. A fixed numerical
        flux divergence would have been reproduced by that change and was not. Model-
        dependent outcome, weaker: the fill-in model predicted a separation of -0.9508%/km
        against a measured -1.1288, a mismatch of -0.178 or 2.1 times the combined
        snapshot scatter of 0.086, so at the parameter-free front factor the residual is not
        common to the pair (D = +0.0211 at sigma = 100 km, +0.1992 at 50 km). The mismatch is
        plausibly in the sharp-front assumption: the modes still unarrived at 55 km at 12 h
        are the same two in both runs (n = 1 at 26.8 h, n = 2 at 13.4 h) and differ only in
        flux share (0.4524 against 0.1381), so in the narrower case the whole correction
        reduces to one 14% step at the top of the band, and the twin's own correction
        overshoots (drift -0.1138 to +0.2003, peak-to-peak 8.66 to 13.13%, nulling front
        factor 0.534 against 0.989 for the wider ridge). Solving for the common front factor
        that reconciles the pair gives f = 1.112 with D = +0.257; with two measurements and
        two unknowns that is a solve rather than a test. Net: the numerical drift over
        10-55 km is bounded between 0.02 and 0.26%/km, improved from about 1%/km.
        Narrowed again 2026-09-02 by extending that twin from 12 to 15 h (3694.87 s
        wall-clock, 73266 steps, converged, dt flat, bit-identical to the 12 h twin over
        all 113 shared snapshots, so an extension and not a repetition). Over 10-55 km in
        the final hour the drift fell from -0.1138 to -0.0332%/km, peak-to-peak from 8.66
        to 4.20% and pointwise drag band RMS from 0.1427 to 0.0815 m/s/sol, so the record
        is flat at the 5% and 0.1%/km thresholds for the first time. Spectrum, grid,
        column, wind and sponge did not move, so a numerical flux divergence fixed in
        height would not have fallen with the horizon and -0.0332%/km is a model-free
        upper bound on it. Both accounts of the drift's time dependence were fitted to
        t <= 12 h and scored out of sample; both are excluded. The relaxation
        a exp(-t/tau) + D_inf is refused by the hour-on-hour increments, needing +0.020 to
        +0.046%/km over 14-15 h at its fitted tau of 2.3-2.5 h against a measured +0.0052
        with a per-snapshot spread of 0.0030. The fill-in account is refused by the level,
        with an out-of-sample error of -0.070%/km at 14 and 15 h, about 23 times that
        spread, and its correction now moves the 15 h drift from -0.0332 to +0.1293. The
        nulling front factor is 0.567 here against 0.534 for the same twin at 12 h and
        0.989 at sigma = 100 km, so it is behaving as a fitted parameter: the f = 0.989
        null is no longer support for the model, and the statistic that 97.0% of the drift
        is arrival transient is not quotable alone, its companions being -76% and -290%.
        A separate check bears on what the residual is. Over 10-55 km, entirely below the
        60 km sponge base, the flux divergence is matched by the measured acceleration of
        the x-mean flow at every hour from 9 to 15 h (band RMS ratios 1.013 to 1.123,
        correlations 0.9964-0.9992, residuals 6.6-13.5%; over the final hour the residual
        is 0.0058 m/s/sol rms, the acceleration slightly exceeding the divergence). That
        bounds loss and not misplacement, since a scheme moving momentum to the wrong
        height conserves it and deposits into du/dt exactly as physical deposition does,
        so the drift in MF(z) is still the only handle and still an upper bound; but the
        band-integrated divergence should not be described as spurious drag while the
        column is filling. Still open: at 15 h the flux median is falling about 1.5%/h and
        the mean flow is accelerating at 0.086 m/s/sol band RMS, so whether -0.035%/km is
        a numerical floor or a slower transient is unsettled, and neither account survives
        to extrapolate it.
        Narrowed differently 2026-09-02 by decomposing the same 15 h record into
        upward- and downward-propagating flux, at no simulation cost (analysis 4.7 s
        against 1-2 s, section 11 of `analyze_linear_flux_constancy.py`, rereading 60 of
        the 733 snapshots). The drift is one projection of the flux profile and the level
        is the other, and a height-uniform loss moves the level while leaving the drift
        untouched, so the "flat at both thresholds" reading above rests on the projection
        that cannot see it. Over 12.0-14.7 h the band flux falls at 1.42%/h while the
        drift stalls at -0.035%/km, which is that combination. Method: with rho_0 the
        x-mean density, Psi = FFT_x[sqrt(rho_0) w'] and Xi = FFT_x[sqrt(rho_0) u'] satisfy
        dPsi/dz = -ikXi + Psi/(2H), and writing Psi = A exp(imz) + B exp(-imz) gives
        Re(Xi Psi*) = -(m/k)(|A|^2 - |B|^2) with the AB* cross terms cancelling
        identically, so the split is pointwise in z and needs no vertical transform. That
        matters because the band holds only about 1.6 vertical wavelengths at the
        band-median lambda_z of 30.7 km and a windowed FFT split would leak the upward
        peak into the downward half. Floor measured before use on a synthetic WKB wave:
        8.4e-05 in flux for a purely upward wave, with imposed fractions of 0.05 and 0.10
        returned as 0.05004 and 0.10002; on the run, R = 0.0486 to 0.0527 across m x 0.90
        to x 1.10 and 0.0491 with background N and U = 40. Closure: the two parts sum to
        <rho_0 u'w'> identically, and against <rho u'w'> the largest level difference is
        1.14% of the flux, the triple correlation <rho' u'w'>; evanescent wavenumbers carry
        1.9e-04. Result: a downward component is present and grows monotonically,
        |MF_dn|/|MF_up| = 0.0248 at 8 h, 0.0410 at 12 h and 0.0486 at 14.74 h, with the
        half-hourly trace running 0.0187 to 0.0494 over 6-15 h at 0.0037 per hour and no
        sign of settling. Its height structure points above the band, rising from 0.037 at
        9.75 km to 0.079 at 54.75 km and 0.18 at 59.75 km just below the sponge base. Two
        independent checks agree: the variance cross term peaks at 15.0 km against the
        predicted lambda_z/2 = 15.3 km at 8 and 12 h (45 km at 14.74 h, so two of three),
        and the sigma = 100 km twin carries R = 0.0138 at 11.74 h with its flux at 0.644 of
        the steady linear value, against 0.0410 and 0.837 for the narrower ridge, which is
        the direction a top-boundary origin predicts. Consequences for the bound above.
        The drift of the upward component alone is +0.0445%/km against -0.0348 for the
        net, so the sign of the residual depends on which wave is being asked about and
        the two are not separable on this record. The sub-band split is weaker than it
        looked: over 40-55 km the drift is -0.195%/km weighting by rho and -0.337
        weighting by its x-mean, a definitional choice with no content at this order, so a
        0.142%/km spread sits on the number, comparable to the whole residual; over 10-55
        km that spread is 0.0005, so the full-band -0.035 is unaffected. Reflection is not
        the main channel: over 12.0-14.7 h the incident flux falls at 1.123%/h and the net
        at 1.417, so the returning wave accounts for 16.9% of the loss, and it accounts for
        21.4% of the 19.47% amplitude shortfall. The other 83% is not a weakening source,
        the x-mean wind at 0.25-9.75 km changing by at most 0.0010 m/s across 8, 12 and
        15 h. On the sponge, K(z) = K0 sin^2(pi/2 (z-60)/30) with K0 = 0.05 1/s set at
        `scripts/sweep_arsia_mons_stage2.jl:129`, against the 0.005 default recorded at
        `src/parameters.jl:87` as the Srivastava optimal; K0/omega_hat is 33-199 over the
        flux-carrying modes and the damping rate overtakes the intrinsic frequency 1.4-3.3
        km above the base, that is within 0.04-0.11 lambda_z, so on the scale the wave
        responds to the layer is closer to a step than to its nominal 0.98 lambda_z
        thickness. Those are scale arguments, not a reflection calculation.
        Qualified 2026-09-03 by resolving that decomposition per horizontal
        wavenumber, again at no simulation cost (`--updown-modes`,
        `--updown-mode-heights`, `--updown-mode-share-floor`; whole analysis about 13 s,
        the per-mode reporting itself free at 13.1 s against 13.9 s with the table off,
        the increase being a second partition formed for every snapshot). fup and fdn
        were already resolved by k inside updown_split() and only their sum was being
        reported; its rfft grid is the grid ridge_spectrum() builds, so ridge mode n is
        column n of both. Each mode then has its own clock: it reaches z at
        t_up = I(z)/(k U^2) and anything reflected at the sponge base reaches the same
        height only at t_rt = (2 I(z_spg) - I(z))/(k U^2), and since cg_z scales with k
        the record ends inside that gap for the long modes. At 9.75 km the round trips of
        n = 1, 2, 3 are 53.15, 26.57 and 17.72 h against a 14.99 h record. Predictions
        fixed first: a sponge origin needs R_k at the floor in those three and orders R_k
        with K0/omega_hat, which runs 33.2 to 397.9 and increases with lambda_x; an
        interior origin needs R_k nonzero wherever t_up has elapsed. Two floors were
        established before reading the result. The per-mode method floor, remeasured on
        the synthetic purely upward WKB wave because the 8.4e-05 quoted above is a
        share-weighted sum that bounds no single column, is 2e-05. That synthetic is
        steady and so does not exercise the one approximation the partition contains,
        dPsi/dz from steady continuity dropping d(rho')/dt in a compressible solver, so
        the derivative was also taken by differencing Psi in z, accurate to about 0.2% at
        dz = 500 m against lambda_z = 31 km; writing the relative difference between the
        two estimates as eps, an error of that size manufactures about eps^2/4 in every
        mode. Measured eps is 0.100-0.169 at 9.75 km, giving a bound of 0.0029, and it
        falls rather than grows over the record. Result, share-weighted over modes and
        with class membership recomputed at each hour: at 9.75 km the modes that cannot
        have received a sponge reflection carry R_k = 0.0073 at 8 h, 0.0162 at 12 h and
        0.0202 +- 0.004 at 14.74 h, a factor of 7 above the closure bound and three orders
        above the floor, and the partition rebuilt on the differenced derivative gives
        0.0160 for the same class. Individually those three read 0.0263, 0.0117 and
        0.0242, so it is not one mode. The scale ordering also runs the wrong way for a
        damping layer, the share-weighted d r_n / d log10 lambda_x over the returned modes
        being -0.107 at 9.75 km, +0.008 at 29.75 and +0.065 at 54.75, where at 54.75 km
        the round trip is only 1.2 times the upward transit and the classes are not
        separated. Replicated on the sigma = 100 km twin, which shares the k grid but
        weights n = 1 at 0.182 against 0.048: its unreturned class carries 0.0079 at 8 h
        and 0.0193 +- 0.008 at 11.74 h against 0.0073 and 0.0162 for the narrower ridge,
        with the two derivative estimates agreeing to 1%, and its four slopes are negative
        at all three heights. So the attribution of the downward component to the damping
        layer is not supported by the timing on either run. Three further corrections.
        The step at the round trip does not replicate: at 9.75 km the returned class
        exceeds the unreturned one on sigma = 50 km (0.0288 against 0.0162) and the
        ordering reverses on sigma = 100 km (0.0137 against 0.0193), so the 41% that a
        face-value reading would put at the round trip is an upper bound, not a
        measurement. The twin comparison quoted above as "the direction a top-boundary
        origin predicts" measures composition rather than origin: the unreturned class is
        the same on both runs and the band ratios differ because 86.6% of the twin's flux
        sits in modes that could not have returned by 11.74 h against 39.1% for the
        narrower ridge at 11.99 h. And R is filling rather than growing, at a lower level
        than the aggregate correction gives. Because R_k is already a per-mode ratio, a
        returned mode carries no filling of its own, and the returned-class mean is close
        to constant in time at every height (0.0321/0.0288/0.0345 at 9.75 km,
        0.0310/0.0330/0.0385 at 29.75, 0.0520/0.0523/0.0535 at 54.75 across 8, 12 and
        14.74 h). Dividing the band R by the share that has had time to return, which
        assumes every mode reflects alike, implies r = 0.0454/0.0538/0.0585; measuring r_n
        directly on the modes that did return gives 0.0384/0.0381/0.0422, about 40% lower.
        The R(z) rise of 0.037 to 0.079 quoted above is also partly composition: modes not
        yet arrived contribute to the numerator of the band ratio and almost nothing to
        its denominator, n = 1 at 54.75 km not having arrived at all (t_up = 26.77 h) yet
        reading R_k = 2.30. The returned-class asymptote runs 0.0345/0.0385/0.0535, a
        factor of 1.55 across the band against 2.12 for the band ratio, which still leaves
        a height dependence a single reflector at 60 km would not give, about 0.98%/km in
        r or 0.49%/km of symmetric attenuation per leg. What survives unchanged from the
        previous entry: a downward component is present, well above the floor of the
        method, and the band ratio grows because the column is filling. What does not: the
        localisation of its source. New figure `linear_flux_permode.png`. Section 11 stays
        additive, the 15 h analysis reproducing all 1039 recorded leaves with 0 changed and
        only `updown/per_mode` added.
        Qualified further 2026-09-03 by testing the assumption the per-mode
        classification rests on, again at no simulation cost (`--transient-stride`,
        `--transient-window`, `--transient-bins`, section 12 of the same script; whole
        analysis 27 s against 15 s with the section off, the increase being 245 further
        snapshot reads at stride 3). t_rt is the round trip of a stationary mountain
        component, ground-based frequency zero and omega_hat = -kU, and the split's
        m = sqrt(N^2/U^2 - k^2) is the same assumption; a component at another frequency
        travels at omega_hat^2/(Nk) instead of kU^2/N and is misassigned by the split,
        which reports a downward share ((1-r)/(1+r))^2 with r = kU/|omega_hat| from a
        purely upward wave. So the unreturned class is unreturned only for waves at the
        mountain frequency. The record tests this: the x-transform of exp(i(kx - omega t))
        carries exp(-i omega t), so a stationary field holds A_k and B_k fixed and
        anything else rotates them at omega = -d arg/dt. Three results. First, the
        measured rates do not open the loophole. At 9.75 km on sigma = 50 km the
        unreturned modes rotate at omega/kU of -0.287, 1.470 and 0.692, giving
        |omega_hat|/kU of 1.279, 0.470 and 0.301 and corrected round trips of 32.1, 120.5
        and 187.3 h against 14.99 h; 3 of the 4 unreturned entries are slowed rather than
        hastened. On the twin the shortest corrected round trip is 21.4 h against 11.99 h,
        with 4 of 6 slowed. Second, the flux splits exactly as <|B|^2> = |<B>|^2 + var(B),
        and only the first part is governed by the stationary clock. Over the last 3 h the
        unreturned class at 9.75 km gives R = 0.0146 on sigma = 50 km and 0.0139 on the
        twin, of which 0.0097 and 0.0095 are steady, against a per-mode floor of 1.7e-05
        and closure bounds of 0.0029 and 0.0045, so the steady part stands at 3.3 and 2.1
        times the closure bound and about 570 times the floor. Across the ten class
        entries on the two runs the fluctuating part of the downward branch is 0.21 to
        0.97 of the fluctuating part of the upward branch, median 0.42, while the steady
        parts stand in a ratio of 0.0027 to 0.0504: entry by entry the fluctuating
        component is divided between the branches 8 to 125 times more evenly, median 39, which is what a partition does with something whose vertical wavenumber is
        not the one it assumes and is not what a reflection does. Third, the growth
        recorded above is not resolved as a trend. Binning the whole record at one hour
        with membership recomputed in each bin gives, for the share-weighted ratio at
        9.75 km, mean 0.01419, sd 0.00545, range 0.00783-0.02281 and a trend of
        +0.00103 +- 0.00052 per hour (2.0 se) on sigma = 50 km, and mean 0.01686,
        sd 0.00614, range 0.00878-0.02569 and -0.00061 +- 0.00068 per hour (0.9 se) on the
        twin; the sigma = 50 km sequence is non-monotone and the two hours quoted above
        are its second-lowest and its highest. The one monotone quantity is the
        fluctuating part, falling 0.00304 to 0.00014 and 0.01110 to 0.00106 across the ten
        bins while tracking the upward branch's fluctuating part, which is the start-up
        settling. Two negatives alongside. No descent: the fluctuating class-summed
        downward flux is simultaneous within one 221 s sample interval over 39.25-29.25 km
        at coherences 0.96-0.98, bounding any propagation from below at 45 m/s against the
        1.16-1.27 m/s stationary group speed of the one mode in that class there, and
        below 29 km the coherence falls to 0.33-0.86 with lags that change sign between
        adjacent levels on the twin. No association with the overturning layer: with
        everything slower than 55 min removed the correlation of the downward flux at
        9.75 km with the unstable area fraction is 0.052 and 0.112 for the class sum and
        -0.385 to 0.186 for individual modes, while the raw coefficients (0.334, 0.293)
        are exceeded by the upward branch (0.705, 0.889), which is filling and cannot be
        driven by the layer, so what raw correlation exists is a shared trend; the
        unstable fraction also peaks at 13.27 h and falls 44% by the end of the sigma = 50
        km record while the ratio does not follow. Net effect on the previous paragraph:
        its central finding stands at about half the level, 0.0097 steady rather than
        0.0202 share-weighted at one hour, the open item on the time dependence is closed
        as scatter rather than explained, and about a third of the quoted ratio is
        reclassified as a component the partition cannot assign a direction to. What is
        left is a downward component steady in the ground frame, at the same level on two
        runs whose weighting of the affected modes differs by 3.8, not growing, not
        descending and not tracking the instability. Also fixed in passing: a local named
        `share` inside the level-decay loop of section 11 shadowed the ridge flux-share
        array of section 8, harmless until section 12 read it; renamed `refl`. New figure
        `linear_flux_transient.png`. Additive, the 15 h analysis reproducing all 1760
        recorded leaves with none removed and none changed beyond 1.7e-11 relative, only
        `transient` added.
        Settled 2026-09-03 by the amplitude control. h = 30 m at sigma = 50 km and 15 h,
        everything else identical, ran in 3684.6 s over 73055 steps against 3695.1 s and
        73266 steps at h = 100 m, confirming dt is set by the vertical acoustic limit and
        not by the wave. R is dimensionless and every candidate left after the previous
        paragraph except one is a linear operator on the wave, so sponge reflection, the
        lower boundary, the start-up residue and truncation on either axis all scale
        numerator and denominator together and leave R fixed; only a source nonlinear in
        the amplitude does not. Prediction fixed before the run finished: the steady
        unreturned-class ratio at 9.75 km either holds at 0.0097 within +-0.0017, making
        the origin linear, or falls toward the 0.0029 closure bound. It fell. Writing each
        quantity as h^p with the incident flux at h^2, so p = 2 means the ratio is
        unchanged: the steady unreturned class gives 0.00969 -> 0.00199 at 9.75 km
        (p = 3.31) and 0.00707 -> 0.00100 at 29.75 km (p = 3.62), both h = 30 m values
        below the closure bounds of 0.0029 and 0.0038 and so unmeasurable there. The
        controls confirm the sensitivity: the fluctuating parts of both branches give
        p = 1.96 to 2.26, five of six within 0.07 of 2, so the component the previous
        paragraph identified as one the partition divides evenly is amplitude-independent
        as that reading requires. The returned class behaves differently, p = 2.14 at
        9.75 km, that is amplitude-independent to 14%, and 2.53 and 2.39 at 29.75 and
        54.75 km. Two amplitudes fix a two-parameter reading, an amplitude-independent
        part plus one quadratic in h, which is what a source quadratic in the wave gives
        since its radiated flux would go as h^4 against the incident h^2: at 9.75 km
        a = 0.00123 and the quadratic part is 0.00846 at h = 100 m, at 29.75 km a = 0.00040
        and 0.00667, so 87% and 94% of what was measured at h = 100 m is amplitude-
        dependent and both offsets sit below the closure bound. Consequence: a downward
        component whose steady part scales as h^3.3 to h^3.6 cannot come from the sponge,
        the lower boundary, the start-up adjustment or truncation on either axis, which is
        the whole candidate list, and what is left is a nonlinearity not yet located.
        Supporting numbers: the steady linear reference scales as h^2 exactly (11.111)
        while the measured incident flux scales as 10.28, so the shortfall against it
        falls from 15.3% to 8.5% and about half of that shortfall is itself amplitude-
        dependent rather than a numerical loss; band R 0.0486 -> 0.0269; R(z) at 10/30/55
        km 0.0371/0.0480/0.0786 -> 0.0257/0.0219/0.0739, so the height dependence inside
        the band largely goes while the value below the sponge barely moves; peak mean-flow
        deviation -2.150 -> -0.206 m/s against -0.194 predicted as h^2; min u/U 0.363 ->
        0.762 against 0.809 predicted by linear scaling, the 6% shortfall being the sense a
        saturating h = 100 m run gives; 40-55 km rho versus rho_0 weighting spread 0.142 ->
        0.057 %/km against 0.043 predicted as h. One prediction fails and withdraws an
        interpretation: the N^2 < 0 area fraction was expected to go to zero and instead
        rose, 0.00207 -> 0.00750, and resolving it in height shows why. In both runs those
        points sit in a narrow 58.75-61.75 km layer straddling the sponge base at 60 km,
        seven levels at h = 100 m and five at h = 30 m, with none below, so that diagnostic
        measures a static instability generated at the top of the damping layer rather
        than wave overturning in the free column and it does not scale with the wave. The
        null result of the previous paragraph is unaffected but the layer should not be
        called the overturning of the mountain wave, and it is not the nonlinearity this
        run detects. Separately, the floor of the method was corrected. synthetic_wkb()
        builds its field from the same m = sqrt(N^2/U^2 - k^2) the split assumes, so it
        cannot see an error in m, and there is one: the anelastic vertical structure
        equation for Psi = sqrt(rho_0) w' at constant N and constant scale height carries
        m^2 = N^2/U^2 - k^2 - 1/(4H^2), and with the band-mean N = 0.00849 1/s and
        H = 9.61 km that term is 6.00% of (N/U)^2, so the assumed m is low by 3.0%. Passing
        the exact solution through the split returns 2.438e-04 over the band and
        2.394e-04 to 2.556e-04 across modes 1-12, matching ((1-r)/(1+r))^2 with
        r = m_exact/m to five digits, while that field's own flux is height-independent to
        0.0000%. This is a systematic bias of one sign in every mode, so it does not
        average down, and it rather than the 1.7e-05 is the floor a per-mode ratio should
        be read against; added as `updown/exact_check` and reported next to the WKB check.
        Nothing reverses, but the margins should be read as 84 times the floor for the
        0.0202 and 40 times for the 0.0097, with the 0.0029 closure bound binding in either
        case. New figure `linear_flux_transient.png` for the h = 30 m run too. Additive,
        the 15 h analysis reproducing all 1760 previously recorded leaves with none removed
        and none differing by more than 1.7e-11 relative, only `transient` and
        `updown/exact_check` added.
        Remaining: locate the nonlinearity, which is now the open question and is the
        cheapest item left. Three parts, none needing a new configuration. A third
        amplitude at h = 10 m, about 3690 s, separates the two readings of the pair: the
        offset-plus-quadratic form predicts 0.00131 for the steady unreturned class at
        9.75 km and a single exponent, R going as h^1.32, predicts 0.00047, a factor of
        2.8, the offset itself being 0.00123 +- 0.00112 and so unresolved on two points;
        both sit above the 2.4e-04 exact floor and below the 0.0029 closure bound, so
        the reading turns on how tight that bound is, and the h = 30 m point indicates it
        is not tight since its measured 0.00199 already sits below it. Going up is not
        available: extrapolating the measured perturbation, which grows slightly more
        slowly than h between 30 and 100 m, min u/U reaches about 0.04 near h = 150 m and
        passes through zero shortly above, so anything much above 100 m overturns; and at
        h = 60 m the two accounts differ by 14%, so an intermediate point is not worth its
        hour. Second, a wave-overturning diagnostic that excludes the damping layer, since
        the one in hand does not: restricting the N^2 < 0 fraction to below 55 km and
        reporting min N^2 rather than the fraction below zero would make it a measurement
        of the wave, at no simulation cost. Third, resolve the triple correlation
        <rho' u'w'>_x by horizontal wavenumber as the split was resolved, to see whether it
        sits in the same modes as the unexplained downward component; one further pass over
        the snapshots. The sponge variation follows, now aimed at the returned class rather
        than the unreturned one and run on the h = 30 m configuration: the returned class
        at 9.75 km has p = 2.14, amplitude-independent to 14%, which is the signature a
        linear return channel has, and at 0.0237 on the h = 30 m run it is the quantity to
        read K0_sponge 0.05 -> 0.005 on. Comparison at equal filled fraction rather than
        equal clock is required, the per-mode class means supply it, and the bin-to-bin
        scatter fixes what counts as a change at +-0.0017 on a ten-bin mean. dz = 250 m
        follows at about four times the 3695 s, ahead of dx because the vertical axis is
        the less resolved for the modes carrying the signal: n = 1, 2, 3 sit at 500, 250
        and 167 points per horizontal wavelength at dx = 4 km, against 61 points for
        lambda_z = 30.7 km at dz = 500 m. dx = 2 km is last and a null result is the
        expected one: horizontal truncation enters as (k dx)^2 and would spread the
        downward spectrum 12^2 = 144 times more than the upward one across n = 1-12, while
        the measured downward spreads are 7.81, 17.75 and 26.75 at 9.75, 29.75 and 54.75 km
        against upward spreads of 10.16, 10.51 and 18.38, a largest excess of 1.7; vertical
        truncation enters as (m dz)^2 and is spectrum-preserving here since k^2/m^2 is at
        most 0.029, and is not excluded by the same test. Its cost is about twice rather
        than four times the 3695 s, the timestep being set by the vertical acoustic limit
        (dz/c = 2.17 s and dx/(c+U) = 14.8 s for c of 173-230 m/s, against dt = 0.737 s).
        Neither run settles the level: the incident flux is still decaying at 1.12%/h under
        a launch wind constant to 0.0010 m/s, height-uniform and so invisible to
        d ln|MF|/dz, and MF_up rather than the drift is the projection that sees it. Note
        also that the h = 30 m run, not the h = 100 m one, is the linear control this
        record has been treating itself as, and flux error bars should be quoted on it.
        Consequence in the meantime: the 12 h nonlinear flux profiles above ~30 km describe
        a filling column rather than a steady state, and flux ratios against the steady
        linear value should be quoted with the arrived fraction (1.00 below ~20 km, 0.82
        over 25-45 km, 0.55 over 50-55 km in the final hour).
        Update 2026-09-22: the sponge variation is under way, with strength (K0 0.005 and
        0.0158), thickness (60-120 km) and base height (70-100 km) run at h = 30 m, 15 h,
        against a linear Rayleigh-layer prediction fixed beforehand and a second estimate
        of the reflected share from the variance ripple that does not use the split. See
        `docs/reflection_convergence_tests.md`.
        Update 2026-09-23: the measured share (0.017 to 0.142) does not follow the layer
        theory; a single-mode test (n = 6) gives 0.03 to 0.05 for two layers whose
        theoretical returns differ by 160 times, and linear throughout. The source of the
        return remains open; the start-up transient is the next candidate (ramped
        forcing, doc Section 8). L1 (36 h) shows the 15 h flat band profile was a
        cancellation between modes: peak-to-peak 9 to 19 % after 17 h. The 1.1 %/h
        incident decay belongs to h = 100 m, not h = 30 m.
        Update 2026-09-23 (later): background audit, doc Section 9. The well-balanced flux
        subtracted P_bar = (gamma-1)*Q_bar[4], which includes the mean-wind kinetic
        energy. With any mean wind this drifts the column (-0.07 %/km in d ln rho/dz),
        loses 0.74 % of the mass through the bottom face, and leaves a near-neutral or
        statically unstable N^2 interface at the sponge base. New switch
        `SimConfig.pbar_mode` (:total is the default and bit-identical; :pressure is the fix;
        `pbar=1` in the Stage 2 sweep script). In the h = 0 pair the fix holds the column
        to round-off. Fixed reruns: the weak-layer excess disappears (C1 0.047 -> 0.005),
        and the strong layer returns 0.019 against a theory value of 0.30 (open). The
        single-mode net-flux peak-to-peak stays at 35 %. All earlier records with wind,
        including Stage 1 and 2, carry the term; Stage 2 eps(z) reaches 2 % in the AMEC band.
        Update 2026-09-23 (closing the linear corner, doc Section 11): :pressure is now the
        SimConfig default (bit-identical to :total at zero wind); the four legacy mean-wind
        scripts pin :total, and the Stage 2 script keeps pbar=0 as its default (new work
        passes pbar=1). Snapshots now carry pbar_mode, K0_sponge and sponge_frac attributes.
        Decisions: ripple is the estimator of record (split floor about 0.006); quote levels
        against band-mean N with the 2.08 N(0) factor alongside; production top 60-120 km,
        K0 = 0.005; K0 = 0.05 retired; strong-layer discrepancy parked with a pre-specified
        isothermal test. Final batch launched with criteria fixed in the doc: ramped pair
        C0R/C1R (tanh, t0 = 2 h, 18 h) and a pbar=1 A/B of Stage 2 U40 h1.0/h3.0 (3 h).
        Stage 1 conclusions carry the same :total flag until they are revisited.
        Results 2026-09-23 (doc Section 11): A/B (U40 h1.0/h3.0, 3 h) keeps the 4.5 km
        absorption (MF(1.25) within 0.07 %) and the AMEC-band flux at 50.25 km (+0.5 %,
        +3.7 %); under :total the 53.3 km critical level fell to 52.0-52.3 km (spurious
        +6 to +10 m/s at 50-55 km) and the flux at 55.25 km was 28-34 % low. Ramped pair
        (C0R/C1R): net peak-to-peak 12-17 h mean 21-22 % against 38 % impulsive, 13-16 %
        at 16-17 h and still falling; not below 10 %, so criterion 1 not met; C0R return
        0.0150 with the same phase (met); C1R level 0.868 against 0.951 (not met).
        Next: tanh ramp on MountainParams forcing (not implemented), then Stage 2
        records >= 18 h with pbar=1, 60-120 km K0 = 0.005 top and a 2 h ramp.
        Update 2026-09-24 (doc Section 12): MountainParams.t_ramp added (0 = bit-identical);
        nine 24 h Stage 2 runs (pbar=1, 60-120 km K0 0.005, 2 h ramp, ~8.8 h wall each; dt
        halved by the -327 m/s MCD wind at 120 km). h = 0.1 km runs complete, steady to
        6-7.5 % in 18-24 h (drop at ~20 h coincides with the jet's 19.5 h wrap time),
        critical level held, ~30 % of the band flux crosses 53.3 km; the opposite-sign
        flux above 4.5 km persists and grows (x1.8) rather than decaying. h = 1 km: U20,
        U40 complete but unsteady (24-33 %), critical level descends to 52.4/51.2 km,
        near-ground mean wind reverses in U40. Aborted (timestep collapse): all h = 3 km
        (at 100-120 km, where the MCD wind is supersonic, |u|/c up to 1.73) and U57 h1.0
        (near 42 km). OPEN: top treatment above ~95 km (lower top / limited wind in the
        layer / stronger top damping) and a smaller-CFL rerun of U57 h1.0.
        Update 2026-09-24 (doc Section 13): top moved to 95 km (layer 60-95 km, K0 0.005).
        C2R single-mode control (11532 s wall): ripple 0.0019 (<= 0.02, met; theory 0.011,
        factor 5.8 below), background at 18 h 0.026 % / min N2 ratio 1.0005 (met), level
        0.910 against theory 0.911 (ratio to C1R 1.048, inside the 10 % band but opposite
        sign). Nine 24 h Stage 2 cases on the 95 km top launched 13:03:44; CFL 0.2 rerun of
        U57 h1.0 (16 h) running since 09:48. OPEN: per-step cost of C2R 1.5x C1R at normal
        priority, cause unknown.
        Update 2026-09-24 evening (doc Section 13 results): 95 km batch ended 18:09 (wall
        16636-18276 s, 0.13 s/step); 7 of 9 complete, U40 h3 and U57 h3 abort at 92-95 km.
        Linear runs steady (7.4-8.1 %), 31-34 % of flux crosses 53.3 km. Top comparison
        FAILS (MF(50.25) 16-43 % lower than 120 km top), but the CFL test shows it tracks
        the time step: CFL 0.2 U57 h1 (dt 0.29 s) collapses at 15.65 h while CFL 0.4 (dt
        0.66 s) completes; at dt 0.29 s the 95 km run follows the 120 km run (dt 0.38 s)
        within 8 % to 8 h. OPEN (blocking Stage 2 amplitudes): time-step convergence,
        proposed U40 h0.1 at CFL 0.4/0.2/0.1 for 9 h (+ C2R at CFL 0.2).
        Update 2026-09-30 (doc Section 14): pre-registration recorded (launch-layer hypothesis
        H_L, readings P1-P6). Driver options dx=, dz=, outint=, ublend= (defaults reproduce
        S04 and C2R bit for bit); launcher scripts/run_dt_convergence_S2T95.sh. Zero-cost
        checks Z1-Z5 all consistent with H_L. Suite S14A-F + S14Z launched 00:34:14, ~8 h.
        NEXT: archive, score (analyze_dt_convergence.py suite; compare_sponge_reflection.py
        --cavity C2R S14C), record readings; then write Section 15 (h = 1 km convergence).
        Update 2026-09-30 (Section 14 "Suite results"): suite complete 05:37 (5.05 h wall,
        19.4 CPU h). S14Z regression met against a matched 500 m run S14Z0. H_L predictions
        failed: P1 MF30 +14.9 % (neither band), P2 fails on R10/30 (zero crossing near 9 km),
        P3 L = 0.52-0.66, P4 uniform-wind fluxes -6.5 % (magnitude up) with reflection
        estimators unchanged, P5 +12-14 %, P6 MF30 -6.4 %. Branch: intrinsic solver dt
        dependence. NEXT (needs approval): switch-level on/off pairs on C2R (split sponge,
        reversed splitting) plus a candidate lower-boundary pair (forced w imposed on the
        face flux, not via the ghost cell before the LW predictor); write predictions first.
        Section 15 (h = 1 km) deferred. Stage 2 absolute fluxes at CFL 0.4 are 29-33 % below
        CFL 0.1 at 30-50 km.
        Update 2026-09-30 (Section 14 "Addendum 2026-09-30, ground-forcing test"): deviation
        recorded (forcing test ahead of the sponge and splitting pairs). SimConfig.bottom_forcing
        (:ghost default, bit for bit; :face = forced w at t^{n+1/2}, ghost rebuilt after the
        x-sweep, ground-face w imposed in the corrector), driver bcface=1. Criteria: G1/G2 change
        in level and Psi0/imp <2 % forcing explains P4, 2-5 % partial, >5 % not the source.
        Expected switch-on level 1.03 (0.99-1.07). Order bands for C2R/S14C/G3. Check passed
        (hash 49e29559450b495a). G1-G3 launched 12:49:47. NEXT: score against the addendum.
        Update 2026-09-30 22:08 ("Ground-forcing test results"): verdict PARTIAL by the
        registered rule (G1->G2 level +0.03 %, Psi0/imp +2.69 %). The face path makes the
        transmitted flux dt-independent (MF10-50 within 0.06 %). The ghost path converges
        toward it (level p_add 1.22, at the first-order edge; extrapolated 1.015-1.021 vs 1.019).
        Psi0/imp is non-convergent (p 0.13) because of a near-ground 2dz structure that follows
        dt, not the forcing path. R_rip -1.19 % in G1 (outside 1 %). The C2R-theory agreement
        is a compensation: the converged level is +11.9 % above theory. NEXT (needs the user's
        decision): the sponge-step and reversed-splitting pairs (registered); ablation of the
        three switch components; a Stage 2 MCD pair with bcface=1; the near-ground structure.
        Update 2026-10-01 (Section 14 "Addendum 2"): bcface=1 adopted for all new runs (default
        stays :ghost); no absolute Stage 2 flux or F/F_lin quoted until steps 1-4 are done.
        Step 1 (zero-cost): face-path level rings over 12-17 h (+0.33 %/h fit, halves 1.011 and
        1.028; G2/G1 constant at +0.03 %); 2dz component has a dt-independent ground amplitude
        (~0.30 of smooth w) and a decay length proportional to 1/dt (0.73, 1.56, 3.68 km);
        theory with the inversion N(z) moves F by -0.4 to +0.3 %, so the +11.9 % excess stays
        unattributed. Two rows added to the Section 11 budget. Code: sponge_step (:split,
        spsplit=1), sweep_order (:zx, zx=1), face_probe (probe=N); check passed (hash
        1fd7c912d485b418). Launched 2026-10-01 14:04:10 (`launch a2`): S14M1/M2 MCD face pair,
        S14K1 split sponge, S14X1 z-then-x, S14R1 dz=250, S14P1-P3 probes; cost history in
        ~/MAGNUSP_runs/logs/S14A2.cost.tsv. Step 4 (H = 0, F/F_lin, Queney on the face path)
        follows before any production or Section 15.
        Update 2026-10-01 20:14 (Addendum 2 results): MCD face pair MF(30.25) +0.35 % per
        halving (ghost path +29 %) -> registered reading: Stage 2 production at CFL 0.4 on
        the face path, after step 4. Split sponge (K1) and z-then-x (X1) nulls met (<=0.05 %).
        dz=250 (R1): level -1.9 % vs G2, at the edge of "dz not the main source"; the +11.9 %
        excess stays unattributed. 2dz mode is set by the vertical Courant number (decay
        0.36 km at dz 250). Probes: component 3 (w imposed in the corrector) carries the
        first-order ghost-path deficit; size bands and stationarity failed. Failed expectation:
        MCD +21 % (measured +0.35 %); the ghost path lost 39 % of the MCD flux at CFL 0.4.
        Residual: MF(1.25) below the 4.5 km critical level still +120 % per halving.
        Slowdowns traced to idle system sleep (pmset sleep 1); `monitor` now holds
        caffeinate -s -i. NEXT: step 4 (H = 0 face equilibrium, F/F_lin, Queney on face path;
        Queney driver needs a bcface/pbar option), then Section 15.
        Update 2026-10-02 (Section 14 "Addendum 3"): step 1 zero-cost checks written first
        (face pair +0.20/+0.34/+0.45 % per hour at MF30; ghost path also reshapes the MCD
        profile, R10/30 -15 %, R50/30 +7 %; 2dz decay length follows L = dz^2/(2 c dt) within
        -9..+17 % in uniform wind; MF(1.25) on the MCD column declared a non-convergent
        diagnostic; two Section 11 budget rows; supersession notes on Sections 11-13 and the
        Stage 1 docs). Code: run_mountain_benchmark_strict.jl --bcface=1 (default ghost, pbar
        pinned :total; --outdir SubString bug fixed), analyze_dt_convergence.py zalt/zflin/
        zequil/zqueney, launcher set a3 (check adds S14M1 0.1 h and the Queney default path;
        hash 98f7fe52db0bfcf7). Batch a3 launched 15:02:37 from the launcher copy, ended
        20:08:00 (5.09 h wall, 17.4 CPU h, no sleep). Outcomes: V1 ridge F/F_lin 1.121 (CFL
        0.4) and 1.122 (CFL 0.2): excess like C2R, attributed to the boundary closure by the
        registered rule (expected 0.86 failed; the attribution rests on the theory's layer
        factor, against the column theory the ridge is within 4 %); V2 H = 0 failed literally
        (uniform at round-off; MCD energy drift +10 %, face-path mass drift 13x, open item);
        V3 Queney face/ghost 1.0499, at the edge, acoustic Courant reading (ghost deficit
        9.3 %; Table 11 metric 2.49 % on both paths vs Srivastava 0.74 %); S1 30 h level 1.020,
        range 0.40 %, met; S2 MCD dz=250 FAILED (MF30 -7.7 %, R10/30 -5.2 %, growing):
        face-path Stage 2 flux not resolved in dz at 500 m; S3 MF(1.25) -6.7 % with MF10/30
        < 0.3 %: not the calm layer alone; S4 z-then-x null met. NEXT (needs approval): Section
        15 (h = 1 km on the face path) and the 24 h nine-case Stage 2 rerun with bcface=1;
        proposals in the session recap, both to carry a dz row or a dz = 250 m twin.
        Update 2026-10-19 (Section 15 pre-registration; addendum 3 V2 note): zero-cost checks
        first. (a) Section 11 row added: dz 500 vs 250 m on the MCD face path, +8.3 % at
        MF(30.25) relative to 250 m. (b) V2 audit (`zequil --budget`): the face-path H = 0 drift
        closes on top outflow plus sponge to 1 %; the 13x mass and +10 % energy difference is a
        ground-face inflow on the ghost reference S14Z0 (velocity, not momentum, mirrored by
        impose_ghost_w! at h = 0 over a stratified ghost background); proposal (not applied):
        S14V2m becomes the face-path equilibrium reference. (c) `zoverturn` on the ghost 24 h
        h = 1 km records: the literal "below 55 km" onset fires at 3.0-5.3 h (calm layer at
        1 km, critical-level region), so the Section 15 window is keyed to the 20-45 km free
        column (ghost onset 19.4 h U40, 12.2 h U57, none U20, at 44-44.75 km, 3-4 km below the
        WKB Lindzen height). (d) Queney Table 11 metric: 1.29 % floor from the x-mean of the
        infinite-domain reference (0.74 % unreachable under this metric; paper definition to
        re-read); the wave part sits at 25-40 km where the steady state is not established
        at 10 h; :total drift <= 0.006 m/s. Registered: F-a to F-c (S15F1/S15F2 h = 1 km dz
        pair, Lindzen, linearity), T1 (24 h CFL 0.2 twin, 1 % at 18-24 h), H3 timings, P
        (PR1 production fluxes, descriptive), V3p (pbar attribution, expected null). Code:
        Queney --pbar=1; analyze_dt_convergence.py zoverturn, zequil --budget, suite --window
        --pair, zqueney --metric-tops; launcher set s15 (run-ID logs, provenance header,
        docs/run_matrix.csv rows, git diff saved, deduplicated archive). dry and check passed
        (hash c503a6c106887e07, 29 snapshots at 0.0). Batch s15: 12 processes, ~13.6-18 h
        wall (S15F2), ~54 CPU h, ~121 GB. OPEN: git state (src/ uncommitted on main; docs/
        and scripts/ git-ignored, so the prereg commit and tag rules cannot be met as is);
        launch waits on that decision.
    (g) Which N does each downstream quantity use? Opened 2026-09-01. N_REF = 0.018 s^-1 is
        the near-surface inversion value (N = 0.01763 at z = 0.25 km), not a lower-atmosphere
        mean; over 10-55 km in this column N is 0.0078-0.0097 with a mean of 0.00849, a
        factor of 2.08 smaller. Audited: `analyze_ext6h_mechanism.py` uses the height-resolved
        N[j] for F_sat = rho k U^3/(2N), for cg_z and for lambda_z, so the Lindzen saturation
        level and the 26-34 km stall prediction are not affected; both sweep drivers and
        `analyze_arsia_stage2.py` use 0.018 only for the printed Nh/U label. Two quantities
        are affected. The first, closed 2026-09-01, is the U=57/h=3 ring-down identification
        in `docs/arsia_stage1_U57_h3_deep_dive.md` Section 5: BUOYANCY_N_REF also fed
        `buoyancy_period_2pi_over_Nref_s` and the wavelet reference lines, and the measured
        388 s at z = 49.75 km was read as the buoyancy period 2 pi/N_ref = 349 s. From the
        run's own background at that height, N = 0.00962 (2 pi/N = 653 s, N/2 pi = 1.532 mHz)
        and the acoustic cutoff gamma g/(2c) is 1.963 mHz, against a measured 2.574 mHz with
        0.257 mHz bins. omega < N holds only below about 1.75 km, so at the station the
        disturbance is on the acoustic branch, not the gravity branch, and the 349 s
        agreement is a coincidence at the 11 per cent level. The timing conclusion is
        unchanged; the mechanism attached to the frequency is now recorded as two open
        readings, an acoustic startup transient or Doppler-shifted advection of a 22.1 km
        structure past the station. `background_frequencies()` was added to
        `analyze_arsia_stage1_deepdive.py`, the wavelet lines now use the local N and local
        cutoff at the station height with the surface value kept as a labelled third line,
        and the metrics carry `local_N_per_s`, `local_acoustic_cutoff_mHz` and
        `band_at_station`. The second is open: the steady linear reference flux scales as
        rho_0 N_0 and is evaluated at N_0 = 0.01763, so MF/MF_linear of 0.637 over the band
        would be 1.32 against a band-N reference. WKB does not hold across the inversion
        (N falls by 2.08 over ~4 km while lambda_z is 14-31 km), so the transmitted reference
        lies between the two and requires solving the vertical structure equation per k with
        a radiation condition aloft. Until then the amplitude ratio must be quoted with its
        N_0. The constancy and fill-in results do not inherit this: for hydrostatic waves
        m = N/U is independent of k (k^2/m^2 ~ 1.8e-3 here), so any partial reflection scales
        every wavenumber alike and cancels from a normalized shape.
  - The Stage 2 U=40/h=3 run was also trending at its 3 h t_end and has not been extended.
    Its imposed 53.3 km critical level sits where the Stage 1 self-induced deceleration
    develops, so a Stage 2 6 h run would test whether the two mechanisms are distinguishable.
  - `scripts/sweep_arsia_mons_stage1.jl ext u= h= hours=` runs any single Stage 1 mode for
    an arbitrary duration; Stage 2 has no equivalent mode yet.

### MEDIUM PRIORITY
- [ ] **Tabulated MCD atmosphere profiles** — replace isothermal with Mars Climate Database rho_0(z), T_0(z), P_0(z)
- [ ] **DSMC coupling file writer** — NetCDF format per planA Phase 2 (n, T, u, w at z_couple)
- [ ] **Mars mountain wave test case** — orographic forcing (Srivastava Eq. 49, Fig. 10)
- [ ] **Mars acoustic pulse test case** — pressure perturbation (Srivastava Eq. 43, Fig. 43)
- [ ] **Critical layer test** — Gaussian x-envelope forcing (Srivastava Eq. 33)
- [ ] **Re-run dispersion test with t_end > 8000s** — with T=239K, c_gz=16.6 m/s; wave center
  reaches sponge (z=140km) at t≈9618s; current 8000s run ends before sponge absorption begins.
  Extend to t=12000s to see full absorption and validate Fig. 6a wave travel diagram.

### LOW PRIORITY / FUTURE
- [ ] **GPU acceleration** — evaluate AMDGPU.jl compatibility; only if CPU simulations take hours
- [ ] **Multi-core parallelization** — threading for Euler sweeps if needed
- [ ] **Tangential velocity in DSMC coupling** — shift mean of tangential Maxwellian samples by u_horizontal
- [ ] **Dispersion test: switch to wavelet analysis** — zero-crossing method gives altitude-dependent
  apparent lambda_z (measurement artifact); wavelet power spectrum matches paper Fig. 6b method

---

## Decisions Log

| Date       | Decision | Rationale |
|------------|----------|-----------|
| 2026-03-03 | Julia 1.12.5 as implementation language | Modern, fast, good scientific computing ecosystem |
| 2026-03-03 | CPU-only initially | GPU compatibility uncertain; parallelize later if needed |
| 2026-03-03 | Earth isothermal validation first | Match Srivastava Fig. 6-7 to verify solver correctness before Mars |
| 2026-03-03 | Newtonian cooling placeholder for CO2 damping | Eckermann coefficients not yet available |
| 2026-03-03 | Julia stdlib sparse + LU for implicit solves | Sufficient for moderate grids (~500x240); avoid extra dependencies |
| 2026-03-03 | NetCDF output via NCDatasets.jl | CF-convention, good for atmospheric science, portable for DSMC coupling |
| 2026-03-03 | 0-120 km vertical domain | Standard thermopause; matches DSMC coupling altitude in planA/planB |
| 2026-03-04 | Sponge damping inside LW z-sweep source term | Avoids operator-split instability; K*dt<<1 so splitting would be fine too but integrated form cleaner |
| 2026-03-04 | Well-balanced Euler: P-P_bar z-flux, (rho-rho_bar)*g source | Makes hydrostatic background exact fixed point; eliminates hydrostatic drift |
| 2026-03-04 | T=239K for all Earth isothermal validation cases | Matches paper Section 3 verbatim; paper's quoted N=0.0182 is inconsistent with this T (see TODO) |
| 2026-03-04 | Dispersion test: Gaussian forcing (Eq. 29/31), not linear | Paper uses Gaussian envelope for dispersion test too (Table 6 lists t0, sigma_t); linear forcing gives secular KE growth |
| 2026-03-04 | Sommerfeld gradient BC rejected | Acoustic instability at T=239K: near-vacuum top cells (rho~1e-10) amplify noise in <400 steps; Eq. 19 retained |
| 2026-03-05 | Viscosity diffusion bug fixed | Domain-average nu caused exponential energy injection: uniform kappa smooths large high-altitude velocities (u∝exp(z/2H)) down to dense cells, inflating density-weighted KE. Fixed by computing nu_z=mu(T)/rho_bar(z) per altitude level and building CN matrix with altitude-dependent coefficients. Confirmed: K0=0 KE now stays bounded after forcing ends (no secular growth). |
| 2026-03-05 | Characteristic BC rejected | Thompson/Poinsot-Lele NSCBC attempted; crashes at step ~14 (t≈17s). Near-vacuum rho~1e-8 at z=120 km means momentum/rho noise amplifies by ~430x per step into multi-km/s velocities → negative P. Perturbation-only decomposition also unstable. Literature confirms no AGW paper uses characteristic BCs; all use sponge layers (physically appropriate). |
| 2026-03-05 | Thicker sponge reduces reflection significantly | Tested 20/40/60 km sponge thickness with K0=0.005: r=19.7%/7.4%/3.74%. Confirms the 20 km sponge violates Klemp et al. (2008) minimum thickness rule (≥1 lambda_z = 46.1 km). For next validation run, extend domain or sponge to satisfy this rule. |
| 2026-03-05 | Literature survey: sponge is the primary BC mechanism | Srivastava (2022) explicitly states Eq. 19 is insufficient alone; sponge does the work. No paper in AGW corpus uses radiation/characteristic BCs at top. Fedorenko (2021) confirms viscosity dominates above ~120 km. Fritts & Alexander (2003) full dispersion relation confirmed consistent with our implementation. |
| 2026-06-30 | Arsia Mons MY31 background extracted from MCD | New scenario for AMEC mountain-wave study (lat -9.5, lon 237, Ls 226, LT 6h, dust=31). `scripts/prep_mars_arsia_mons_atmosphere.py` parses the cached MCD column (from `fetch_web_form_data` preset) into `external/mcd_profiles/mars_arsia_mons_my31.csv` (cols z_m,rho,P,T,nu,u_zonal; 0-120 km @ 500 m). Wire-in: `atm_file=` for z,rho,P,T; pass u_zonal as `u_init` (Vector over grid.z_full) and its surface value as the mountain `U_bar`. |
| 2026-06-30 | MCD z=0 temperature is ground skin temp; corrected to air-consistent value | Per-level R_eff=P/(rho*T) is ~191.1 J/kg/K at every level except z=0 (211.6, an outlier). Reset T(z=0) from 167.2 K (cold pre-dawn surface skin) to 187.3 K = P/(rho*R) so the wave-launching bottom layer is air-consistent; restores the surface density match with raw MCD. |
| 2026-06-30 | Keep solver R=188.9 (pure CO2) for Arsia Mons background | MCD effective R~191.3 (realistic CO2+N2+Ar mix) not adopted, to avoid a mars_config change. Well-balanced hydrostatic reconstruction with R=188.9 preserves the MCD thermal structure; rebuilt density diverges from raw MCD only near the top (up to ~57% at 120 km, inside the sponge). Zonal wind reverses sign at ~4.5/53/78 km (critical levels for a stationary mountain wave). |
| 2026-06-30 | Isothermal/constant-wind cap profile generated to 156 km (reserved for a later thermospheric stage) | A 155 km domain was first attempted so the ~120 km DSMC coupling level would sit below the sponge. The profile `external/mcd_profiles/mars_arsia_mons_my31_cap156km.csv` (313 levels, ztop=156 km to cover top ghost cells) holds T and U flat above the MCD top (120 km) and integrates P hydrostatically (clamping P/rho would break hydrostatic balance and inject a spurious pressure gradient). This 155 km run blew up (see smoke-test row below); the cap profile is retained for a future thermospheric stage that controls wave-amplitude growth. |
| 2026-06-30 | Smoke test revealed thermospheric wave breaking; Stage 1 top lowered to 90 km (AMEC focus), 120 km coupling deferred | A 1 h smoke run (U=40, h=1 km) on the 155 km domain blew up at t~991 s. Cause is physical, not a bug: the stationary mountain wave amplifies as rho^-1/2 into the thermosphere (peak \|w\| ~5 m/s at 50 km -> 375 m/s at 110 km -> supersonic ~800 m/s at 120 km) and breaks, because the viscosity cap (nu<=1e4) throttles the molecular viscosity that physically damps these waves above ~100 km (cf. Fedorenko 2021). The mid-atmosphere (<=80 km) is clean, including a healthy ~5 m/s updraft at the 45-53 km AMEC level. Decision: cut the domain to 90 km (sponge 60-90 km) for the AMEC study; defer the 120 km coupling to a later stage that raises the nu cap to restore realistic thermospheric viscosity. Added a dt-collapse abort guard in `src/run.jl` so a breaking run bails out instead of churning ~1e-12 s steps. |
| 2026-06-30 | Arsia Stage 1 baseline: uniform eastward steering wind + outflow x-BC, real T/rho/P; sweep U_bar x h_max | `scripts/sweep_arsia_mons_stage1.jl`. Domain LX=2000 km (dx=4 km, nx=500), LZ=90 km (dz=500 m, nz=180); fixed Gaussian hill sigma=100 km (Arsia-scale); sponge_frac=0.333 (sponge 60-90 km, K0=0.05); CFL=0.4; viscosity on, cooling off; atm_file=`mars_arsia_mons_my31.csv` (0-120 km real MCD covers the 90 km domain). Sweep U_bar in {20,40,57} m/s and h_max in {0.1,1,3} km (linear -> nonlinear onset, Nh/U~1 near h_max=2.8 km for U=50, N=0.018). The flat-bottom forcing w=U_bar*dh/dx cannot represent the literal ~20 km obstacle (Nh/U~7, blocked); h_max is an effective forcing amplitude. U(0)=0 in the profile, so U_bar is a deliberate choice (eastward magnitude; the real westward jet is mirrored for outflow-BC consistency). |
| 2026-06-30 | Real reversing u_zonal shear deferred to Stage 2 | The outflow x-BC assumes a unidirectional eastward wind at all heights; the real U(z) reverses sign (4.5/53/78 km), which it mishandles. Stage 2 will pass u_init as a Function closure (interpolating u_zonal onto grid.z_full; a Vector must match nz_total incl. ghosts) with periodic-x and a wider domain so wakes do not wrap before sponge absorption, activating the critical levels to test whether jet-launched waves deposit momentum near the ~45-53 km AMEC cloud level. |
| 2026-06-30 | Arsia Stage 1 smoke test and full 3x3 sweep validated on the 90 km domain | The `output/arsia_mons_stage1_smoke/` directory inherited from the prior session was stale 155 km data (NetCDF `z=310`, `ng=2` per `src/parameters.jl`, so interior nz=306, not the 90 km config's nz=180; confirmed via `ncdump`), not an actual 90 km result as the prior session's notes assumed. Cleared it and re-ran. Smoke test (U_bar=40, h_max=1 km, t_end=1 h): converged, dt stayed near 0.73 s throughout (no collapse), 4920 steps, wall-clock 219 s. Full sweep (U_bar in {20,40,57} m/s, h_max in {0.1,1,3} km, t_end=3 h): all 9 runs converged, including the most nonlinear corner (U_bar=20, h_max=3 km, Nh/U=2.70) that was expected to be the likeliest to break. Per-run wall-clock ranged 1214-1622 s (13519-17510 steps); the 9 runs executed concurrently across 16 threads, so total elapsed time was close to the slowest single run (~27 min), not the summed 12367 s. Output: `output/arsia_mons_stage1/run_U<U_bar>_h<h_max>/` + `sweep_results.csv`. |
| 2026-06-30 | Arsia Stage 1 sweep analyzed; hill-region AMEC amplitude tabulated and an inflow-edge artifact isolated | `scripts/analyze_mars_mountain_sweep.py` extended with an opt-in `--composite` / `--composite-only` mode (`compute_band_response`, `plot_sweep_composite`, fast final-snapshot loader). Outputs in `output/arsia_mons_stage1_analysis/` (per-run figures + `metrics.json`, `sweep_composite_w.png`, `amec_response.csv/.md`). Hill-region (\|x_rel\|<=500 km) AMEC-band mean \|w\| increases monotonically with both U_bar and h_max; the composite shows upstream phase tilt with height (vertically propagating orographic gravity wave). A persistent inflow-edge feature (~10-50 m/s, scales with U_bar, ~independent of h_max) sits ~900 km upwind of the summit; it inflated the unrestricted peak only, is outside the hill window, and is an additional reason for the Stage 2 periodic-x BC. The U=57/h=3 km corner (Nh/U=0.95) does not settle in 3 h (amplitude still rising at t_end, \|w\|~50 m/s ~ U), consistent with near-overturning. Full results in the section below. |
| 2026-06-30 | Arsia Stage 2 code authored (real reversing shear + periodic-x); to run on the M2 | `scripts/sweep_arsia_mons_stage2.jl` (runner) and `scripts/analyze_arsia_stage2.py` (analysis). Keeps the Stage 1 eastward forcing (U_bar in {20,40,57} m/s, h_max in {0.1,1,3} km) so any change is attributable to the background, but replaces the uniform wind with the real reversing `u_zonal(z)` installed as a `u_init` Function closure (interpolates CSV col 6 onto `grid.z_full`, clamped at endpoints) and switches `bc_x_type` to `:periodic` (removes the Stage 1 inflow-edge artifact). No solver change required: the sheared background is a discrete steady state (x-uniform gives zero x-sweep tendency; w=0 and P=P_bar give zero z-sweep flux/source) and `diffusion_step!` diffuses only u' = u - u_bar(z), so the jet persists. Domain width LX is a separate sensitivity axis {2000,3000,4000} km (the `width` run mode at a representative U=40, h=1 km), not crossed into the U x h grid. Vertical column unchanged (90 km, sponge 60-90 km); the 53.3 km critical level sits below the sponge base. The analyzer adds the momentum-flux profile MF(z)=<rho u'w'>_x (EP stress; ~constant where freely propagating, drops where deposited) with the 4.5/53.3/77.7 km critical levels overlaid. Both validated (Julia parse OK; analyzer run end-to-end on synthetic data). Details in the Stage 2 plan section below. |
| 2026-07-13 | Stage 2 analysis executed; critical-level absorption at 53.3 km observed in all 9 runs, ~85% of jet-level momentum flux deposited across the 45-56 km layer | `analyze_arsia_stage2.py` run over `output/arsia_mons_stage2/` (miniforge base python; analysis wall-clock 4.3 s). Composite shows wave activity terminating at the 53.3 km critical level in every forced panel (quiet above, unlike Stage 1 where waves filled the column to the sponge). Most launched flux is absorbed below ~5 km at the 4.5 km critical level, so the analyzer's `mf_deposited_frac` (referenced at 10 km, where MF is near zero or sign-flipped) is unreliable for the h>=1 km runs; a multi-level MF table referenced at 30 km (mid-jet) gives a uniform 0.83-0.86 deposition fraction across 45-56 km and >0.98 by 60 km. Hill-region AMEC mean \|w\| is reduced relative to Stage 1 at matched forcing for U>=40 (e.g. U=57/h=3: 0.86 vs 3.99 m/s; hill peak 16.4 vs 49.6 m/s) and the run now appears settled rather than still growing. Periodic-x removed the inflow-edge artifact (global peak = hill peak in all runs). Width sweep (`stage2_width`) not yet run. Full results in the Stage 2 section below. |
| 2026-07-13 | Stage 2 width-sensitivity sweep run and analyzed; hill-region results width-independent to ~5-12% across LX=2000-4000 km | `sweep_arsia_mons_stage2.jl width` (3 runs at U=40/h=1, LX in {2000,3000,4000} km, t_end=10800 s; all converged; wall-clock 770/1080/1498 s, ~linear in LX; total elapsed 1506 s with 3 concurrent runs). Analyzer run with `--width-mode` (12.5 s) -> `output/arsia_mons_stage2_width_analysis/`. Hill-window RMS w agrees to ~5% (0.211-0.221 m/s, non-monotonic); AMEC hill mean \|w\| 0.206/0.198/0.184 m/s (~11% spread, weakly decreasing with LX); deposition fraction (30->56 km, multi-level MF) 0.90/0.82/0.85. MF as a domain mean dilutes as 1/LX; the width-integrated flux MF*LX agrees to ~10%. No wrap-around signal: far-field (\|x_rel\|>700 km) RMS w is 3-5x below the hill window at all widths. Conclusion: the LX=4000 km baseline results are width-robust, and LX=2000 km reproduces the same physics within ~10% at half the cost (candidate grid for longer runs). Results appended to the Stage 2 section below. |
| 2026-07-13 | Deep dive on the Stage 1 U=57/h=3 km corner; station oscillations identified as the startup transient, late wavelet burst as Doppler-shifted small scales | Report: `docs/arsia_stage1_U57_h3_deep_dive.md`; new script `scripts/analyze_arsia_U57h3_deepdive.py` (renamed `scripts/analyze_arsia_stage1_deepdive.py` on 2026-08-31) (deepdive_* figures + metrics in `output/arsia_mons_stage1_analysis/run_U57.0_h3.0/`, 16.5 s wall-clock). Key findings: the `:mountain` forcing is a steady BC active from t=0 for the whole run (no ramp, no shutoff); the station ring-down (dominant period 388 s ~ 2pi/N_ref = 349 s) is the impulsive-start transient and settles by ~3200-4000 s, after which the stationary wave appears as a near-constant offset. The `wavelet_w.png` feature near minute 165 is a 2-8 mHz burst made prominent by the 120-min trim + max-normalization; apparent frequencies above the acoustic cutoff correspond to lambda_x ~ 7-28 km structures advected at U=57 m/s (Doppler shift), matching the final-snapshot striations. Mean flow decelerates by up to 22 m/s at 57.75 km with MF increasingly negative and still trending at t_end, consistent with early-stage steepening/deposition below the sponge base. Snapshot cadence 62.65 s caps resolvable frequencies at ~8 mHz. [Revised 2026-09-01: the frequency identification in that finding is superseded. At z = 49.75 km the local N/2 pi is 1.532 mHz and the local acoustic cutoff 1.963 mHz, both below the measured 2.574 mHz, so the ring-down is not a buoyancy oscillation at the station. The timing conclusion is unaffected. See item (g) and Section 5 of the report.] |
| 2026-08-31 | U=40/h=3 confirmed non-converged at the 3 h sweep t_end in both stages; 6 h Stage 1 extension run, which does not settle but continues into overturning with a descending decelerated layer | Target changed from U=57/h=3 to U=40/h=3 at the user's request (cleaner wave field at U=40, h=1 and 3 km). Baselines already existed in `output/arsia_mons_stage1/run_U40.0_h3.0` and `output/arsia_mons_stage2/run_U40.0_h3.0`, so no new 3 h runs were needed. A momentum-budget check on both showed MF and mean-flow drift still trending at t_end, so the same motivation as U=57/h=3 applies. New `ext` mode added to `scripts/sweep_arsia_mons_stage1.jl` (args `u=`, `h=` in km, `hours=`) writing to `output/arsia_mons_stage1_ext6h/`; run wall-clock 2434.1 s for 32253 steps to t=21600 s, single-threaded combination, 1.17 GB of snapshots at 74.0 s cadence. `scripts/analyze_arsia_U57h3_deepdive.py` generalized and renamed to `scripts/analyze_arsia_stage1_deepdive.py` (new `--run-dir/--compare-run/--u-bar/--label/--profile-hours` args and a new `deepdive_convergence.png` fixed-height time-series panel with per-window linear trends); analysis wall-clock 193.3 s over 478 snapshots. Outcome: the extension is not a converged steady state. The extension reproduces the 3 h baseline to within ~1% over the shared interval, then between 3 h and 6 h the AMEC-band hill mean \|w\| grows from 0.79 to 4.06 m/s, the AMEC-band mean zonal wind falls from 40.8 to 22.2 m/s, and the static-instability fraction in the hill window below 60 km rises from 0.004 to 0.086. MF at 30/45/52 km approximately plateaus after ~4 h (last-hour normalized trend +0.16 to +0.29 per hour with CoV 0.06-0.16, against -0.67 to -0.79 per hour over 2-3 h) while the mean-flow deceleration continues at ~3 m/s per hour, a combination consistent with wave saturation rather than with a steady non-breaking solution. Full results in the section "Arsia Mons Stage 1 U=40/h=3 time extension" below. |
| 2026-08-31 | Time-series evidence figures built for the 6 h extension; two documented values corrected in the process | Three figures added under `output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/`, each from `deepdive_series.npz` rather than from the snapshots: `deepdive_trendfit.png` (least-squares fits over the 2-3 h and 5-6 h windows drawn and annotated, the evidence for the quoted trends), `deepdive_hovmoller.png` (time-height mean-flow deviation, MF and hill RMS w, with labelled deceleration contours and the 5 m/s layer base traced), and `deepdive_breaking.png` (full 323-snapshot indicator record plus a time-height instability map). Scripts `scripts/viz_ext6h_{trendfit,hovmoller,breaking}.py`, about 5 s combined. Two corrections resulted. First, the static-stability diagnostic had used kappa = 0.2248 from assumed constants; recomputed with the solver's own Mars values (gamma = 1.31, R_gas = 188.9, so kappa = 0.2366) the unstable fraction is 0.004 at 3 h and 0.086 at 6 h rather than 0.005 and 0.094, and the U=57/h=3 comparison value is 0.063 rather than 0.067. Second, the descent of the decelerated layer had been read off a profile figure by eye as reaching ~20 km at 6 h; measured against a 5 m/s threshold with z < 2 km excluded, the layer base is 57.3 km at 3 h and 35.3 km at 6 h, so the AMEC band is enclosed from about 5 h rather than 4.5 h. The exclusion is necessary because the 0.75 km level holds a persistent deviation reaching -5.7 m/s that is a separate near-surface feature. The instability map independently supports the descent, showing instability from ~58 km at 3.3 h spreading to ~30 km by 6 h. It also shows the unstable fraction at 6 h is 0.512 at 58.75 km against 0.132 at 49.75 km, about 3.9x higher next to the sponge base, which quantifies the sponge caveat recorded above. First crossings from the full record: local flow reversal at 1.68 h, unstable fraction above 0.001 at 2.55 h, hill-window max \|w\| above U_bar at 4.94 h. |
| 2026-08-31 | Descending instability layer in the 6 h extension attributed to wave saturation; two configuration influences bounded | `scripts/analyze_ext6h_mechanism.py` (1.5 s, reads only `deepdive_series.npz` and the MCD column; metrics in `mechanism_metrics.json`) tests the 6 h U=40/h=3 record against analytical results using the solver's own constants. Findings: (i) hill-window RMS w grows at the conservative (rho N)^-1/2 rate to within 1% over 45-59 km at 1 h, and at 50-66% of it below 45 km, consistent with a wave field still filling the column; (ii) the stationary-wave group velocity cg_z = kU^2/N is 1.5-1.8 m/s, so the steady mountain wave reaches only ~37 km by 6 h, while measured flux fill-in times run inverted with height (1.16 h at 57 km against 1.55 h at 10 km), placing the early upper-column signal in the impulsive startup transient; (iii) the overturning criterion u_hat >= U (equivalently min_x u <= 0) is first met at 59.75 km at 1.68 h and descends monotonically to 29.75 km by 4.00 h, with N^2 < 0 following each level by 0.6-1.9 h; (iv) d<u>/dt equals -(1/rho) d<rho u'w'>/dz to within 2-18% at the level of peak deceleration from 2 h to 5 h, with a positive mid-level residual whose sign matches the lateral influx through the non-periodic outflow x-BC; (v) the layer base is reproduced to within 0.5-2.0 km from 3-5 h by the Lindzen saturation level \|MF\| = rho k U^3/(2N) with k_eff = 9.01e-6 1/m from the ridge flux spectrum, and its descent rate to within 10-30% by the kinematic contour speed -(du/dt)/(du/dz); (vi) F_sat evaluated against the supplied low-level flux places the asymptotic breaking height near 33 km (25.8 km for the full steady linear flux of 22.2 mPa, of which the run reaches 49% at 6 h), giving a testable prediction that the base stalls at 26-34 km rather than descending further. Bounded configuration influences: 58-67% of the flux entering at 25 km is deposited between 25 and 40 km at all times while the sponge-adjacent 56-60 km band takes 0.2-4.8%, so the sponge cannot account for the bulk response; but `src/run.jl:57` sets `Q_bar = copy(Q)` so the sponge pins <u> at 40.0 m/s above ~65 km at 6 h against 14.4 m/s at 56.75 km, the 57.5-60 km band carries 3-5x more static instability per unit resolved wave slope than 45-57 km (so the 0.604 peak at 58.25 km is not a breaking measurement), the first overturning altitude is bounded by the top of the sponge-free column, and the horizontal scale carrying the saturated flux falls to 11 points at dx = 4 km by 6 h. Details in `docs/ext6h_instability_mechanism.md`. |
| 2026-09-01 | 12 h Stage 1 extension run to test the 26-34 km stall prediction; prediction partially supported, and a ~7 h validity limit on Stage 1 domain-mean diagnostics identified | The `ext` mode of `scripts/sweep_arsia_mons_stage1.jl` now names its output directory after the horizon (`hours=6` still resolves to `..._ext6h`, so existing paths are unaffected), and `scripts/analyze_ext6h_mechanism.py` derives its sampled hours from the series instead of nine hardcoded lists capped at 6 h, with one snapshot interval of tolerance so a run ending at 5.991 h still reports a 6 h row. Regression on the 6 h series: no keys removed, no values changed, three added. Run: `julia --project=. scripts/sweep_arsia_mons_stage1.jl ext u=40 h=3 hours=12`, wall-clock 6253.5 s for 72922 steps to t=43200 s, converged, 730 snapshots, 2.5 GB. No dt collapse (dt 0.43-0.55 s at the end), so no zero-wind critical level formed. Deep-dive extraction 382.9 s over 730 snapshots (against 193.3 s over 323, i.e. linear in snapshot count); figures and mechanism report about 6 s combined. Throughput fell from 21.3 to 11.7 steps/s over the run while per-step cost is dt-independent, which points at snapshot output rather than the solver. The solver has no `@threads`/`@spawn` in `src/`, so single-mode runs are serial and deterministic; the 12 h series is bit-identical to the 6 h series over the full 5.991 h overlap on all five cached fields and on the time axis (max abs diff 0.000e+00). Outcome on the prediction: the layer base occupies the predicted 26-34 km band from 6.5 h (33.25 km) to 10.5 h (26.25 km), with dz/dt reversing sign in four of nine samples (+0.40 to +1.51 km/h at 8-10 h), then falls to 11.25 km by 12 h while the saturation-level diagnostic moves the other way (36.25 km at 9 h to 44.75 km at 11 h). The premise behind the 25.8 km asymptote fails: the flux at 25 km plateaus at 5.3-6.4 mPa, i.e. 24-29% of the 22.23 mPa steady linear value, rather than approaching it, and the equilibrium crossing height stays at 32.2-36.2 km from 4-11 h. Reframing: the altitude of peak deceleration stays at 50-55 km from 5-12 h while the base descends, so after ~5 h the structure is better described as a layer deepening downward from a roughly stationary top than as a layer descending as a whole. Levels that decelerated most then recover (54.75 km from 14.54 to 26.00 m/s and 39.75 km from 30.70 to 36.44 m/s between 6 h and 12 h, both well below the 60 km sponge base) while 9.75 and 24.75 km lose 5.5 and 5.1 m/s. The overturning criterion continues descending (29.75 km at 4.00 h, 19.75 km at 6.23 h, 9.75 km at 9.91 h) and the unstable fraction becomes column-wide (0.204 at 5 km, 0.132 at 10 km, 0.124 at 25 km at 12 h). Validity limit: the x-mean momentum budget closes to within 20% from 2-6 h then degrades to a ratio of 0.25-0.70 over 7-12 h, with the residual in the sign expected from lateral influx of undecelerated air. The advective transit from the inflow edge to the ridge is LX/(2U) = 6.94 h, coinciding with the failure between the 6 h and 7 h samples. Independently, the domain spans only 2.87 flux-weighted wavelengths (2000 km against lambda_x = 697 km), and the mean over a 700 km ridge window minus the far-field mean reaches 26.2 m/s at 49.75 km and 18.0 m/s at 39.75 km at 12 h, alternating sign with height on the ~15 km vertical wavelength, against 6.5 and 9.7 m/s at 3 h. Conclusions from the 6 h record sit inside this window; the post-10.5 h base collapse does not and is treated as a diagnostic artifact. Local extrema also become implausible late: the minimum hill-window u reaches -113.5 m/s at 11 h (Mach 0.58 against c_s = 195.7 m/s) and isolated samples exceed the local sound speed from 9.54 h, consistent with grid-scale noise where breaking is unresolved and with the dt reduction from 0.55 to 0.43 s. The minimum zonal-mean wind below 60 km, by contrast, bottoms at 13.4 m/s at 5.7 h and recovers to a steady ~22 m/s from 8-12 h without approaching zero, so no self-induced critical level forms in the mean flow. Two earlier qualifications weaken: the 57.5-60 km instability excess per unit wave slope falls from 4.6x at 6 h to 1.06x at 12 h, and the horizontal resolution of the flux-carrying scale at 55 km improves from 12.8 to 98.6 points as the wind recovers, so 6 h is the worst-resolved point of the record rather than the endpoint. Deposition stays bottom-heavy throughout (49-69% over 25-40 km, 0.2-5.9% in the sponge-adjacent 56-60 km band). Outputs in `output/arsia_mons_stage1_ext12h/` (snapshots, `sweep_results.csv`, `run_log.txt`) and `output/arsia_mons_stage1_ext12h_analysis/run_U40.0_h3.0/` (`deepdive_series.npz`, `deepdive_{trendfit,hovmoller,breaking}.png`, `deepdive_metrics.json`, `mechanism_metrics.json`, `mechanism_report.txt`). Analysis in `docs/ext6h_instability_mechanism.md` sections 11-16. Next run indicated: extend Stage 2, which is periodic in x and closes the budget without a lateral term. |
| 2026-09-01 | Linear control corner run to measure the solver's own momentum-flux error bar; flux is not height-independent, but the record is not converged and the departure is bounded above rather than measured | Rationale: every flux figure quoted so far (the 24-29% plateau at 25 km, the 58-67% deposition partition over 25-40 km, the 0.83-0.86 Stage 2 fractions) was reported without an error bar, and the production configuration violates the assumptions of the analytic results it is compared against in several ways at once. At h = 0.1 km the steady linear solution gives u_hat/U = 0.45 at 55 km, so the overturning criterion is never met below the sponge, the wave is conservative, and Eliassen-Palm requires dMF/dz = 0; molecular viscosity from the MCD nu column removes only 0.03% of the flux by 55 km, so any drift is numerical. Neither existing h = 0.1 km run serves: the Stage 1 one (uniform wind, outflow-x) carries max\|w\| = 20.8 m/s against a linear 0.65 and MF of -80.4 uPa at 16.8 km against a steady linear -24.70 uPa, and the Stage 2 one (periodic, real shear) has its flux absorbed at the 4.5 km critical level, settling near 1% of the linear value. New `linear` mode added to `scripts/sweep_arsia_mons_stage2.jl` (uniform u_init = U_bar with periodic-x; args `u=`, `h=` in km, `hours=`, `lx=`; horizon- and width-named output directory because MF is a domain x-mean and scales as 1/LX, giving -24.703 uPa at 2000 km against -12.509 uPa at 4000 km). `ridge_spectrum()` and `linear_steady_flux()` in `scripts/analyze_ext6h_mechanism.py` parameterized for reuse; regression on the 12 h series: 235 keys, none added, removed or changed. New product `scripts/analyze_linear_flux_constancy.py`. Run: wall-clock 2500.7 s for 58564 steps to t=43200 s, converged, 586 snapshots, 2.0 GB, LX = 2000 km at the Stage 1 grid; whole analysis chain under 72 s. An earlier launch was lost when the SMB share unmounted at about 33% of the integration; all four edited files were re-verified intact afterwards. dt held 0.7370-0.7400 s across 117 diagnostic samples, a spread of 0.41%, which supports attributing the 0.55 to 0.43 s reduction in the nonlinear 12 h run to the flow rather than to the configuration, since everything except ridge height and x-BC was held fixed. Validity of the control: min_x u below 60 km never falls under +29.735 m/s (0.743 U), max \|<u>_x - U\| at 55 km is 0.227 m/s, and measured lambda_z from the phase tilt of w matches 2 pi U / N(z) to a band mean of 0.996 (individual levels 0.89-1.17), with the one-sided cospectrum summing back to mean_x(u'w') to 4.4e-16. The 8.36 m/s record maximum of hill-window \|w\| occurs at 0.02-0.10 h and is the impulsive switch-on transient; at t_end the amplitude is 0.69-0.85 of the linear prediction. Outcome: over 10-55 km in the final hour MF varies by 50.58% peak-to-peak about a median of -16.377 uPa with a fitted drift of -1.2426% per km, so the flux is not flat at the stated 5% and 0.1%/km thresholds. The record does not support attributing that to the solver: peak-to-peak falls monotonically from 6.822 at 1 h to 0.476 at 12 h and the drift from -7.09 to -1.12% per km, and d ln\|MF\|/dt over the final 3 h is -1.5%/h at 4.75 km and +1.2%/h at 9.75 km against +5.4%/h at 29.75 km, +11.1%/h at 39.75 km and +20.0%/h at 54.75 km. Arrival of the wave front (55 km by 9.25 h from cg_z = k_eff U^2/N) is therefore not the same as settling of the flux, and 12 h is sufficient below about 10-20 km and insufficient above about 30 km. Spurious drag -(1/rho) dMF/dz over the final hour is -0.0095 m/s/sol at 10 km, -0.0119 at 25 km, -0.128 at 40 km and -0.868 at 55 km (band RMS 0.321), the same ordering as the spin-up. At the settled levels MF/MF_linear is 0.784 at 4.75 km and 0.770 at 9.75 km, consistent with the 0.853 \|w\| ratio squared; this report does not separate a solver amplitude loss from an error in the analytic normalization, so that ratio is provisional in a way the constancy measurement is not. Most useful result, because it needs neither run converged: differencing the control against the nonlinear 12 h run at the same time on the same grid, with each profile normalized by its own steady linear value and then by its own 10 km value, gives an excess that scatters around unity below 30 km (0.71-1.37, no trend) and falls monotonically above 35 km to 0.61 at 39.75 km, 0.52 at 44.75 km, 0.21 at 49.75 km and 0.099 at 54.75 km. The results are consistent with saturation acting above about 35 km and do not support attributing the flux decrease below about 30 km to saturation from this record, which bears on the 58-67% deposition partition over 25-40 km since that window spans the boundary. Broader implication: the 12 h nonlinear flux profiles above ~30 km describe a filling column rather than a steady state. Outputs in `output/arsia_mons_linear12h_LX2000{,_analysis}/`. Analysis in `docs/linear_flux_constancy.md`; open item (f) above. |
| 2026-09-01 | Height dependence of the control's momentum flux benchmarked against the arrived part of the ridge spectrum; the drift is consistent with an arrival transient rather than a numerical flux divergence | Motivation: extending the 12 h control would not have settled the question, because the fill-in is not an exponential approach to a steady state. For a stationary hydrostatic component cg_z = k U^2/N is linear in k, so the transit time t_arr = I(z)/(k U^2) with I(z) = int_0^z N dz' scales as 1/k, and the Gaussian ridge spectrum k\|h_k\|^2 has a low-k tail with no natural cutoff; the drift therefore decays like a spectral tail, and a run three times longer would leave the same ambiguity at three times the cost. No new simulation was run. Method, added as section 8 of `scripts/analyze_linear_flux_constancy.py`: substituting the WKB solution already used by `linear_steady_flux()` into the one-sided cospectrum gives per-wavenumber contributions F_k = -rho_0 N_0 U w_k k \|h_k\|^2, each height-independent on its own and summing to the steady value to 0.0e+00 in this run, so each mode can be given its own arrival clock and the benchmark becomes MF_fill(z,t) = sum over k >= I(z)/(U^2 t) of F_k. `arrival_profile()` was corrected to integrate N rather than average cg_z between levels, which is the exact form of int dz/cg_z and matters because N falls by 2.08 across the near-surface inversion; arrival times shift by under 2%, 55 km moving from 9.25 to 9.33 h. Spectrum: in the 2000 km periodic domain five bins carry 95% of the flux (lambda_x = 2000, 1000, 667, 500, 400 km with shares 0.182, 0.271, 0.248, 0.166, 0.085), and at 55 km after 12 h the two longest, 45.2% of the flux, have not arrived, reaching that height at 26.8 h and 13.4 h. Arrived fraction over the band in the final hour: 1.00 below ~20 km, 0.818 over 25-45 km, 0.548 over 50-55 km. Outcome: correcting each snapshot by its own arrived fraction before averaging moves the fitted drift over 10-55 km from -1.2426 to +0.0372% per km (accumulated -55.9% to +1.7%) and the band-fit spurious drag from -0.008 / -0.029 / -0.114 / -0.466 m/s/sol at 9.75 / 24.75 / 39.75 / 54.75 km to +0.0002 / +0.001 / +0.004 / +0.026. Multiplying every arrival time by f gives -0.877% per km at f = 0.5, +0.037 at f = 1 and +3.303 at f = 2 (the last with part of the band below the arrived-fraction floor), and the f that nulls the drift is 0.989; since f is not a free parameter, a null within 1% of unity indicates the WKB group velocity places the fronts where the record puts them. Reservations: the residual is still an upper bound, containing whatever the sponge, the dispersive spread about each front and the WKB group velocity contribute; and the residual peak-to-peak of 0.279 is not a solver figure, since the measured profile is smooth in height while the benchmark steps at four or five discrete arrival heights, so within this section the fitted slope carries information and the scatter does not. The excess table differencing the control against the nonlinear 12 h run is unchanged, because the arrived fraction depends on the launched spectrum, column, wind and elapsed time, all shared between the two runs, so it cancels from a ratio of ratios; saturation moving the nonlinear run's flux-carrying scale to higher k would make it the more filled of the two, so the shortfall attributed above 35 km is if anything understated. Separately, the amplitude ratio was found to depend on the reference more than on the measurement: MF_linear scales as rho_0 N_0 and uses N_0 = 0.01763 s^-1 at the first cell centre, so the band mean of 0.637 would be 1.32 against the band-mean N of 0.00849; WKB does not hold across the inversion, and net Eliassen-Palm flux is conserved through a non-dissipating partially reflecting layer, so a shortfall against a single-N launch formula is not a departure from dMF/dz = 0. The constancy and fill-in results do not inherit this, because m = N/U is k-independent for hydrostatic waves (k^2/m^2 ~ 1.8e-3) and any reflection factor cancels from a normalized shape. N audit recorded as open item (g). Analysis chain rerun in 4.9 s; no simulation. Sections 8 to 13 of `docs/linear_flux_constancy.md` revised, with sections 10 (fill-in) and the sigma = 50 km proposal new; cross-references in `docs/ext6h_instability_mechanism.md` sections 14 and 16 updated. |
| 2026-09-02 | sigma = 50 km twin of the linear control run as a differencing partner rather than a longer horizon; most of the flux drift is confirmed to be spectrum-dependent, and the numerical component is bounded between 0.02 and 0.26% per km | Correction to the earlier framing: narrowing the ridge does not move any mode's arrival time, since t_arr = I(z)/(k U^2) contains k and not sigma. It reweights the flux from the slow modes to the fast ones, raising the arrived fraction at 55 km at 12 h from 0.548 to 0.862, while m = N/U and the steady amplitude N h (rho_0 N_0/rho N)^0.5 both carry no k, so lambda_z, the sponge thickness in wavelengths, dz per wavelength and the linearity margin are identical between the two runs. The pair is therefore a controlled difference of the same kind that already worked for the nonlinear excess table. Predicted before launching from `per_k_steady_flux()` and `fillin_flux()` unchanged (2 s): transient drift -1.2637%/km at sigma = 100 km against -0.3130 at 50 km, a separation of -0.9508 against a per-snapshot scatter of 0.085, so measurable. `sigma=NN` (km) added to the linear mode of `scripts/sweep_arsia_mons_stage2.jl`; the default resolves to 100e3 exactly and the suffix is appended to the output directory only when it differs, so existing run directories do not move. Run: 2587.66 s wall-clock, 58609 steps, final_t = 43200 s, converged, timestep flat at 0.737-0.740 s (no reduction, as in the sigma = 100 km run, which supports reading the nonlinear run's timestep reduction as a property of the flow), 587 snapshots, 2.0 GB into `output/arsia_mons_linear12h_LX2000_sig50/`. A first launch was lost at step 37700 when the network share disconnected; the partial output was removed and the run repeated. Extraction 28.6 s, constancy analysis 2.2 s. Controls confirmed unchanged in the twin: lambda_z band mean measured/predicted 1.0200 against 0.996, min_x u = 0.543 U with no overturning, max \|w\| below the sponge 8.53 m/s, mean-flow departure at 55 km 0.310 m/s. Model-free result: the measured drift over 10-55 km moved from -1.2426 to -0.1138%/km (factor 10.9), peak-to-peak 50.6 to 8.66%, accumulated -55.9 to -5.12%, band-fit spurious drag at 54.75 km -0.466 to -0.080 m/s/sol, MF/MF_linear band mean 0.637 to 0.842. A fixed numerical flux divergence would have been reproduced by a change that moved only the flux weighting across k, and was not. Model-dependent result, weaker: predicted separation -0.9508 against measured -1.1288, a mismatch of -0.178 or 2.1 times the combined snapshot scatter of 0.086, so at f = 1 the residual is not common (D = +0.0211 and +0.1992). The sharp-front assumption is the likely source: the modes still unarrived at 55 km at 12 h are the same two in both runs (n = 1 at 26.8 h, n = 2 at 13.4 h) and differ only in flux share (0.4524 against 0.1381), so the whole correction in the narrower case reduces to one 14% step at the top of the band, and the twin's own correction overshoots (-0.1138 to +0.2003%/km, peak-to-peak 8.66 to 13.13%, nulling front factor 0.534 against 0.989). A common front factor of f = 1.112 reconciles the pair with D = +0.257, but with two measurements and two unknowns that is a solve rather than a test. Section 9 (`--pair-npz`, `--pair-sigma`) added to `scripts/analyze_linear_flux_constancy.py`, verified inert without a pair (35 keys, 0 changed) and exact when a run is paired with itself. `docs/linear_flux_constancy.md` sections 6, 12, 13 and 14 revised. Next: a dx = 2 km repeat at sigma = 50 km is now the decisive check, since it does not depend on the fill-in model. |
| 2026-09-02 | Full-record wavelet panels in the Stage 1 deep dive redrawn as two logarithmic views with the cone of influence marked; the minute-165 burst is qualified rather than removed | The panels used `power /= power.max()` on a linear jet map, and the startup transient is two orders above everything after it. Measured at S3, z = 49.75 km in the sigma = 50 km control: record peak 47.0 in power units, strongest feature after 6 h 2.015, only 13.7 dB down, which renders at 0.04 of full scale and reads as empty; late median 53.7 dB down. Replaced with a 2x1 figure at the same filenames. Top panel is 10 log10(power/peak) on a fixed -50 to 0 dB scale, so amplitudes stay comparable across time and between the two station heights. Bottom panel divides each frequency by its own median over the last WAVELET_BG_FRAC = 0.5 of the record, which surfaces a weak late feature at its own level but discards amplitude, so it is blanked below WAVELET_REL_FLOOR_DB = -40 dB absolute; without that blanking the amplified noise floor produced tens of dB of banding out of nothing, which the first draft showed. The cone of influence, already returned by `wavelet()` and previously discarded, is now shaded on both. Consequences: in the linear control the late record is empty in the wave band and not elsewhere. Above 1 mHz, where the startup transient sat, the strongest reliable power after 6 h is 48.5 dB below the record peak, an amplitude bound of 0.020 m/s against the 5.377 m/s transient. The only surviving late power is at 0.065 mHz, a 4.3 h period at -22.4 dB, which is not a wave but the slow approach of the still-arriving long modes; the time series confirms it without the transform, being 97.7% a linear ramp at -0.066 m/s per h over t > 6 h with a residual sd of 0.018 m/s and lag-1 autocorrelation 0.998. That is the predicted signature: a stationary wave has zero phase speed and so gives a constant w at a fixed station, and the residual drift is the 13.8% of the flux that has not arrived at that height by 12 h. In the U=57/h=3 run the minute-165 burst survives outside the cone at z = 49.75 km near minute 160 at 3-5 mHz but lies almost entirely inside it at z = 57.25 km, so its apparent growth to t_end is not separable from the transform edge and the upper-station strength should not be read from these panels. Recorded in Section 6 of `docs/arsia_stage1_U57_h3_deep_dive.md`. Deep dive rerun 17.5 s (U57/h3) and 26.8 s (the sigma = 50 km control); `deepdive_metrics.json` unchanged. |
| 2026-09-02 | Linear control twin at sigma = 50 km extended from 12 to 15 h to score the two available accounts of the flux drift's time dependence out of sample; both are excluded, and the drift falls inside the reporting thresholds for the first time | Motivation: sections 10 and 13 both rested on the fill-in account, which had been fitted to the record it was built from rather than tested by it, and a three-parameter relaxation a exp(-t/tau) + D_inf with no wave physics in it describes the same twelve hourly drifts about as closely (rms 0.026-0.030%/km). The two agree over the record and separate beyond it. Both were fitted to t <= 12 h and projected forward before launching, as new section 10 (`--forecast-from`, `--forecast-to`, `--forecast-fit-lo`) of `scripts/analyze_linear_flux_constancy.py`, with `relax_fit()` solving the two linear parameters exactly on a scan in tau and rejecting any window whose tau exceeds its own length. Two corrections to the framing came out of that prediction. The arrived fraction being frozen at 0.8619 over 9-13.4 h holds at the evaluation top only: mode n = 2 enters the band at 2.84 h and reaches 55 km at 13.39 h, so its front crosses one dz = 500 m level about every 0.1 h throughout, and T moves from -0.401 to -0.283%/km over 9-12 h. Consequently the step in T when n = 2 clears the band top is +0.0145%/km, not the +0.2 a step in the arrived fraction at one height would give, because the band-fitted slope responds to fronts inside 10-55 km rather than to the fraction arrived at its top; the separation is carried by the 12-15 h increment, +0.150 against +0.067%/km, and not by a step. Run: `julia -t auto --project=. scripts/sweep_arsia_mons_stage2.jl linear u=40 h=0.1 hours=15 sigma=50`, wall-clock 3694.87 s (61.6 min) for 73266 steps to t = 54000 s, converged, dt flat at 0.737-0.740 s with no reduction, 733 snapshots, 2.5 GB into `output/arsia_mons_linear15h_LX2000_sig50/`; extraction 32.5 s, analysis 1 s. Cost per step is about 14% above the 12 h twin's 2587.66 s / 58609 steps, consistent with analysis work sharing the machine. It is an extension and not a repetition: all 113 snapshots shared with the 12 h twin are bit-identical on every stored field, so hours 0-12 are the same integration. Controls unchanged: min_x u below the sponge +14.529 m/s (0.363 U, no overturning), max \|w\| below the sponge 8.529 m/s, mean flow within 0.320 m/s of U at 55 km, and static instability at t_end confined to 58.75-61.75 km with zero unstable area fraction anywhere in 10-55 km. Outcome, forecasts against measurement for the drift of the hour-mean profile over 10-55 km: at 13/14/15 h the fill-in account gave -0.0435/+0.0317/+0.0367 and the relaxation family -0.0887 to -0.0302 / -0.0579 to +0.0301 / -0.0379 to +0.0762, against measured -0.0644/-0.0383/-0.0332. The relaxation is refused by the hour-on-hour increments, which need no anchor for D and no choice of fit window: at the fitted tau of 2.3-2.5 h it requires +0.020 to +0.046%/km over 14-15 h, and the measured increment is +0.0052 against a per-snapshot spread of 0.0030, with the half-hourly trace reading -0.0331, -0.0314, -0.0349. The drift has stopped falling, and at about -0.035 rather than at the D_inf of -0.000 to +0.223 the family projects. The fill-in account is refused by the level: its out-of-sample error is -0.070%/km at both 14 and 15 h, about 23 times the per-snapshot spread, with an rms of 0.0584 over the three scored hours against 0.0123 for the family median, and applied to the 15 h record its correction now moves the drift from -0.0332 to +0.1293, accounting for -290% of it against -76% at 12 h. The nulling front factor is 0.567 at 15 h against 0.534 for the same twin at 12 h and 0.989 at sigma = 100 km, so it is behaving as a fitted parameter; the f = 0.989 null should no longer be read as support for the model, the f = 1.112 pair reconciliation stands as the solve it was described as, and the statistic that 97.0% of the drift is arrival transient is not quotable alone, its companions being -76% and -290%. What the extension establishes without either account: over 10-55 km in the final hour the drift moved from -0.1138 to -0.0332%/km, peak-to-peak from 8.66 to 4.20%, accumulated departure from -5.12 to -1.49%, and pointwise drag band RMS from 0.1427 to 0.0815 m/s/sol, so the record is flat at the 5% and 0.1%/km thresholds for the first time. Since spectrum, grid, column, wind and sponge did not move, a numerical flux divergence fixed in height would not have fallen with the horizon, and -0.0332%/km is a model-free upper bound on it. Separately, the flux divergence over 10-55 km, entirely below the 60 km sponge base, is matched by the measured acceleration of the x-mean flow at every hour from 9 to 15 h: band RMS ratios 1.013/1.054/1.081/1.123/1.058 with correlations 0.9964-0.9992 and residuals 6.6-13.5%, and over the final hour the residual is 0.0058 m/s/sol rms with mean -0.0013, the acceleration slightly exceeding the divergence rather than falling short. That bounds loss and not misplacement, since a scheme transporting momentum to the wrong height conserves it and deposits into du/dt exactly as physical deposition does, so the drift in MF(z) remains the only handle and remains an upper bound; but it does mean the band-integrated divergence should not be called spurious drag while the column is still filling. Still open: the run is not steady at 15 h, with the flux median falling about 1.5%/h and the mean flow accelerating at 0.086 m/s/sol band RMS, so whether -0.035%/km is a numerical floor or a slower transient is not settled, and neither account survives to extrapolate it. Operational note: the analysis environment was rebuilt because Homebrew had moved to Python 3.14 and the 3.13 interpreter carrying numpy was gone; `~/.venvs/magnusp` (uv, CPython 3.13.15, numpy 2.5.2, scipy 1.18.1, matplotlib, netCDF4 1.7.4) reproduces the recorded sigma = 100 km analysis exactly, 36 keys and 0 differences, and section 10 is additive, adding only `horizon`. `docs/linear_flux_constancy.md` section 14 added and section 15 revised. Next: the dx = 2 km repeat, now at sigma = 50 km and 15 h. |
| 2026-09-02 | Momentum flux of the 15 h linear control split into upward- and downward-propagating parts, on the record already in hand; a returning wave carrying 4.6% of the flux is measured, growing, and consistent with a top-boundary origin, and the "converged" reading of the previous entry does not survive | Motivation: the drift d ln\|MF\|/dz is one projection of the flux profile and the level is the other, and a height-uniform loss moves the level while leaving the drift untouched, so section 2 cannot see it. Over 12.0-14.7 h the band flux falls at 1.42%/h while the drift stalls at -0.035%/km, which is that combination, and below the sponge the geometry admits one channel for such a loss. Method, added as section 11 of `scripts/analyze_linear_flux_constancy.py` (`--updown-hours`, `--updown-window`, `--updown-max-snaps`, `--updown-msens`, `--updown-synthetic`, `--sponge-k0`): with rho_0 the x-mean density, Psi = FFT_x[sqrt(rho_0) w'] and Xi = FFT_x[sqrt(rho_0) u'] satisfy steady anelastic continuity dPsi/dz = -ikXi + Psi/(2H), so the vertical derivative comes from the stored u rather than a difference in z; writing Psi = A exp(imz) + B exp(-imz) gives Re(Xi Psi*) = -(m/k)(\|A\|^2 - \|B\|^2) with the AB* cross terms cancelling identically in the real part. Three consequences: the split is pointwise in z and needs no vertical transform or window, which matters because the 45 km band holds only about 1.6 vertical wavelengths at the band-median lambda_z of 30.7 km and a windowed FFT split would leak the upward peak into the downward half; the 1/(2H) term drops out of the flux exactly; and the flux carries no cross term at all, so a partial standing wave produces no modulation in MF(z) and the interference has to be looked for in the variance rather than in the flux residual. Floor fixed before use on a synthetic WKB wave built from the same ridge spectrum, density and N(z): 8.4e-05 in flux for a purely upward wave, set by the WKB amplitude gradient, with imposed reflected fractions of 0.05 and 0.10 returned as 0.05004 and 0.10002. On the run the dependence on m is quadratically small as predicted, R = 0.0486 with m from the mean state, 0.0491 with background N and U = 40, and 0.0527/0.0499/0.0486/0.0497 at m x 0.90/0.95/1.05/1.10; the net flux does not depend on m and its drift is -0.0348%/km in every row, which is a check on the partition. Closure: the two parts sum to <rho_0 u'w'>_x identically, so the test is against <rho u'w'>_x, where the largest level difference is 1.14% of the flux, the triple correlation <rho' u'w'>_x; evanescent wavenumbers carry 1.9e-04. Result: \|MF_dn\|/\|MF_up\| over 10-55 km is 0.0248 at 8.00 h, 0.0410 at 11.99 h and 0.0486 at 14.74 h (down shares 2.42/3.94/4.64%), with the half-hourly trace running 0.0187 to 0.0494 over 6.0-14.99 h at a linear 0.0037 per hour and no sign of settling. The height structure points above the band: 0.037 at 9.75 km, 0.079 at 54.75 km, 0.18 at 59.75 km just below the sponge base, tending to 1.00 above 70 km where the net flux is 0.16 uPa and only a standing residual remains. Two independent corroborations: the variance cross term peaks at a vertical scale of 15.0 km against the predicted lambda_z/2 = 15.3 km at 8 and 12 h, with predicted depths 2r/(1+r^2) of 0.31 and 0.39 against measured peak-to-peak 0.47 and 0.75 (at 14.74 h the dominant scale moves to the 45 km band length, so two of three); and the sigma = 100 km twin gives R = 0.0138 at 11.74 h with its flux at 0.644 of the steady linear value, against 0.0410 and 0.837 for the narrower ridge, less of its spectrum having reached the sponge. Consequences for the previous entry. The drift of the upward component alone is +0.0445%/km against -0.0348 for the net, both inside the 0.1%/km threshold but opposite in sign, so the residual is at least as consistent with the height structure of the returning wave as with a divergence in the incident one and the two are not separable here. The sub-band slope is not a measurement at the level it has been read: over 40-55 km it is -0.195%/km weighting the covariance by rho and -0.337 weighting by its x-mean, a 0.142%/km spread on a definitional choice with no content at this order, against 0.0005 for the full band, and removing the returning wave does not remove the non-monotonicity (-0.083/+0.060/-0.270 across the thirds). Reflection is not the main channel: over 12.0-14.7 h the incident flux falls at 1.123%/h and the net at 1.417, so the returning wave accounts for 16.9% of the loss, and it accounts for 21.4% of the 19.47% shortfall against the steady linear reference (net 19.47%, upward alone 15.31%), the remainder still carrying the factor 2.08 N_0 ambiguity of item (g). The other 83% is not a weakening source: the x-mean wind at 0.25, 1.75, 4.75 and 9.75 km changes by at most 0.0010 m/s across 8, 12 and 15 h while the incident flux falls at about 1.1%/h. Sponge scales, reported as scale arguments and not as a reflection calculation: K(z) = K0 sin^2(pi/2 (z-60)/30) from `src/sponge.jl` with K0 = 0.05 1/s set at `scripts/sweep_arsia_mons_stage2.jl:129`, an order of magnitude above the 0.005 default recorded at `src/parameters.jl:87` as the Srivastava optimal; K0/omega_hat = 33-199 over the flux-carrying modes, and the damping rate overtakes the intrinsic frequency 1.4-3.3 km above the base, 0.04-0.11 lambda_z, so on the scale the wave responds to the layer is closer to a step than to its nominal 0.98 lambda_z thickness (0.14-0.36 lambda_z at K0 = 0.005). Verification: section 11 is additive, the sigma = 100 km analysis reproducing 447 recorded leaves with 0 changed at rel_tol 1e-9 and only `horizon` and `updown` added. New figure `linear_flux_updown.png`. `docs/linear_flux_constancy.md` section 15 added and section 16 reordered. Next: the sponge variation at fixed grid, ahead of the dx = 2 km repeat, because a reflected component of this size would differ between two runs sharing a boundary treatment but not a grid and would be attributed to truncation. |
| 2026-09-03 | Up/down decomposition of the linear control resolved per horizontal wavenumber, on the records already in hand; the attribution of the downward component to the sponge is not supported by the timing on either run, and the previous entry's twin corroboration is reinterpreted as composition rather than origin | Motivation: the previous entry localised the downward component at the sponge base on geometry and on scale arguments about the damping layer, neither of which was tested, and the decisive test was already latent in the data: `updown_split()` forms fup and fdn per horizontal wavenumber and only their sum over k was being reported, while its rfft grid is the grid `ridge_spectrum()` builds, so ridge mode n is column n of both. Method, added as `--updown-modes`, `--updown-mode-heights` and `--updown-mode-share-floor` of `scripts/analyze_linear_flux_constancy.py` with the helper `mode_timing()`: each mode has its own clock, reaching z at t_up = I(z)/(k U^2) and receiving anything reflected at the sponge base only at t_rt = (2 I(z_spg) - I(z))/(k U^2), and since cg_z scales with k the record ends inside that gap for the long modes, the round trips of n = 1, 2, 3 to 9.75 km being 53.15, 26.57 and 17.72 h against 14.99 h. Predictions fixed before the table was produced: a sponge origin needs R_k at the floor in those three and orders R_k with K0/omega_hat, which runs 33.2 to 397.9 and increases with lambda_x; an interior origin needs R_k nonzero wherever t_up has elapsed. Modes are classed on their own clock at each hour, counted only once the sample time exceeds 1.5 times their transit, and modes below 1e-3 of the steady flux are excluded from the aggregates because their ratio is set by the noise of the partition; on the sigma = 100 km ridge that floor changes the returned-class error bar from +-0.086 to +-0.002 and it leaves the sigma = 50 km numbers untouched, its smallest mode carrying 0.017. Two floors established first. The per-mode method floor, remeasured on the synthetic purely upward WKB wave because the 8.4e-05 of the previous entry is a share-weighted sum bounding no single column, is 2e-05. That synthetic is steady and so does not exercise the partition's one approximation, dPsi/dz from steady continuity dropping d(rho')/dt in a compressible solver, so a second partition was built on dPsi differenced in z, accurate to about 0.2% at dz = 500 m against lambda_z = 31 km; with eps the relative difference between the two estimates, an error of that size puts B = -eps Psi/2 and manufactures about eps^2/4 in every mode, the same quadratic law the m_scale sensitivity measures. Measured eps at 9.75 km is 0.100-0.169, a bound of 0.0029, and it falls rather than grows over the record while the signal triples. Result: at 9.75 km the modes that cannot have received a sponge reflection carry share-weighted R_k = 0.0073 at 8 h, 0.0162 at 12 h and 0.0202 +- 0.004 at 14.74 h, a factor of 7 above the closure bound and three orders above the floor, with the differenced-derivative partition giving 0.0160 for the same class and the three modes reading 0.0263, 0.0117 and 0.0242 individually. The scale ordering also runs the wrong way for a damping layer, the share-weighted d r_n / d log10 lambda_x over the returned modes being -0.107 at 9.75 km, +0.008 at 29.75 and +0.065 at 54.75, the last at a height where the round trip is only 1.2 times the upward transit and the classes are not separated. Replicated on the sigma = 100 km twin, which shares the k grid but weights n = 1 at 0.182 against 0.048: its unreturned class carries 0.0079 at 8 h and 0.0193 +- 0.008 at 11.74 h, its two derivative estimates agree to 1%, and its slopes are negative at all three heights. Three corrections to the previous entry follow. The step at the round trip does not replicate: at 9.75 km the returned class exceeds the unreturned one on sigma = 50 km (0.0288 against 0.0162) and the ordering reverses on sigma = 100 km (0.0137 against 0.0193), so the 41% a face-value reading would put at the round trip is an upper bound and not a measurement. The twin comparison quoted as "the direction a top-boundary origin predicts" measures composition rather than origin, the unreturned class being the same on both runs while 86.6% of the twin's flux sits in modes that could not have returned by 11.74 h against 39.1% for the narrower ridge at 11.99 h. And R is filling rather than growing, at a lower level than the aggregate correction gives: because R_k is already a per-mode ratio a returned mode carries no filling of its own, the returned-class mean is close to constant in time at every height (0.0321/0.0288/0.0345 at 9.75 km, 0.0310/0.0330/0.0385 at 29.75, 0.0520/0.0523/0.0535 at 54.75 across 8, 12 and 14.74 h), and dividing the band R by the share that has had time to return implies r = 0.0454/0.0538/0.0585 against 0.0384/0.0381/0.0422 measured directly, about 40% lower. The R(z) rise of 0.037 to 0.079 is also partly composition, modes not yet arrived contributing to the numerator of the band ratio and almost nothing to its denominator, n = 1 at 54.75 km not having arrived at all (t_up = 26.77 h) yet reading R_k = 2.30; the returned-class asymptote runs 0.0345/0.0385/0.0535, a factor of 1.55 across the band against 2.12, still leaving a height dependence a single reflector at 60 km would not give, about 0.98%/km in r or 0.49%/km of symmetric attenuation per leg. A spectral shape test was run alongside and does not settle the question: fitting \|MF_dn,k\| = C \|MF_up,k\|^alpha gives alpha = 0.31 +- 0.21, 0.64 +- 0.34 and 0.57 +- 0.44 on sigma = 50 km and 0.85 +- 0.20, 0.74 +- 0.21 and 0.91 +- 0.43 on the twin, so five of six are within 1.5 standard errors of the alpha = 1 a reflection requires, and the twin's fits have the longer lever arm, 139-fold against 10-fold in the denominator. What survives from the previous entry: a downward component is present, well above the floor of the method, and the band ratio grows because the column is filling. What does not: the localisation of its source, which is now open between a distributed reflector, the lower boundary, the start-up adjustment and the interior scheme. Verification: purely additive, the 15 h analysis reproducing all 1039 recorded leaves with 0 changed and only `updown/per_mode` added. Wall clock about 13 s for the whole analysis against 4.7 s, the per-mode reporting itself free (13.1 s with it, 13.9 s with `--updown-modes 0`), the increase being the second partition formed for every snapshot; no simulation time. New figure `linear_flux_permode.png`. `docs/linear_flux_constancy.md` section 16 added and the old section 16 renumbered to 17 with its items reordered. Next: the dx = 2 km repeat is now ahead of the sponge variation, the reason for deferring it having been that a reflected component would differ between runs sharing a boundary treatment but not a grid, which the timing removes by making the interior scheme a live candidate; prediction stated first, that an interior origin should drop the unreturned-class 0.0202 by about a factor of four to the closure bound while a boundary or start-up origin leaves it unchanged. |
| 2026-09-03 | Stationarity of the downward component tested on the records already in hand; the per-mode classification of the previous entry is not overturned by the frequencies actually present, its central result stands at about half the level, and the growth it left open is not resolved as a trend | Motivation: the previous entry classes a mode by t_rt = (2 I(z_spg) - I(z))/(k U^2), which is the round trip of a stationary mountain component with ground-based frequency zero and omega_hat = -kU, and the split of section 11 applies m = sqrt(N^2/U^2 - k^2), which is the same assumption. A component at the same k and another frequency travels at omega_hat^2/(Nk) rather than kU^2/N, so its round trip is the tabulated one times (kU/omega_hat)^2, and the split reports a downward share ((1-r)/(1+r))^2 with r = kU/\|omega_hat\| from it even when it is purely upward, that share rising toward 1 as r moves away from 1 in either direction. The unreturned class was therefore unreturned only for waves at the mountain frequency, a structural limit of the classification rather than a defect in its measurement, and the one place the previous entry could fail without any of its numbers being wrong. Method, added as section 12 of `scripts/analyze_linear_flux_constancy.py` (`--transient-stride`, `--transient-window`, `--transient-bins`, helpers `transient_pass()` and `tr_fig()`; `updown_split()` now also returns A, B and mk): the x-transform of exp(i(kx - omega t)) carries exp(-i omega t), so a stationary field holds A_k and B_k fixed and anything else rotates them at omega = -d arg/dt, measured as -arg<B_{i+1} B_i*>/dt over the trailing window with the modulus normalised by the geometric mean of the two energies so a growing amplitude cannot push it past 1. The flux being quadratic, <\|B\|^2> = \|<B>\|^2 + var(B) splits it exactly into a part a fixed amplitude would carry, which the stationary clock governs, and a part it could not; the same split on A is the control. Cost: the whole analysis 27 s against 15 s with the section off, 245 further snapshot reads at stride 3, the only part of the script whose cost scales with the length of the record; no simulation time. Result 1, the loophole does not open. At 9.75 km on sigma = 50 km the unreturned modes rotate at omega/kU of -0.287, 1.470 and 0.692, giving \|omega_hat\|/kU of 1.279, 0.470 and 0.301 and corrected round trips of 32.1, 120.5 and 187.3 h against a 14.99 h record, so two of the three descend more slowly than assumed rather than faster; 3 of 4 unreturned entries are slowed. On the sigma = 100 km twin the shortest corrected round trip is 21.4 h against 11.99 h, with 4 of 6 slowed. Result 2, about a third of the ratio is not steady and is not returning flux. Over the last 3 h the unreturned class at 9.75 km reads R = 0.0146 on sigma = 50 km and 0.0139 on the twin, of which 0.0097 and 0.0095 are steady, against a per-mode floor of 1.7e-05 and closure bounds of 0.0029 and 0.0045, so the steady part is 3.3 and 2.1 times the closure bound and about 570 times the floor. Across the ten class entries on the two runs the fluctuating part of the downward branch is 0.21 to 0.97 of the fluctuating part of the upward branch, median 0.42, while the steady parts stand in a ratio of 0.0027 to 0.0504: entry by entry the fluctuating component is divided between the two branches 8 to 125 times more evenly, median 39, which is what a partition does with something whose vertical wavenumber is not the one it assumes and is not what a reflection does. Result 3, the growth is scatter. Binning the whole record at one hour with class membership recomputed in each bin gives, for the share-weighted ratio at 9.75 km, mean 0.01419, sd 0.00545, range 0.00783-0.02281 and a trend of +0.00103 +- 0.00052 per hour (2.0 se) on sigma = 50 km, and mean 0.01686, sd 0.00614, range 0.00878-0.02569 and -0.00061 +- 0.00068 per hour (0.9 se, wrong sign) on the twin; the sigma = 50 km sequence runs 0.0123, 0.0095, 0.0078, 0.0097, 0.0121, 0.0220, 0.0228, 0.0114, 0.0143, 0.0200 and the two hours the previous entry quotes, 0.0073 at 8 h and 0.0202 at 14.74 h, are its second-lowest and its highest. The one monotone quantity is the fluctuating part, falling 0.00304 to 0.00014 and 0.01110 to 0.00106 across the ten bins while tracking the upward branch's fluctuating part, which is the start-up settling. Two negatives alongside. No descent: the fluctuating class-summed downward flux is simultaneous within one 221 s sample interval over 39.25-29.25 km at coherences 0.96-0.98, bounding propagation from below at 45 m/s against the 1.16-1.27 m/s stationary group speed of the one mode in that class there, and below 29 km the coherence falls to 0.33-0.86 with lags that change sign between adjacent levels on the twin. No association with the overturning layer: with everything slower than 55 min removed the correlation of the downward flux at 9.75 km with the unstable area fraction is 0.052 and 0.112 for the class sum and -0.385 to 0.186 for individual modes, while the raw coefficients 0.334 and 0.293 are exceeded by the upward branch at 0.705 and 0.889, which is filling and cannot be driven by the layer, so what raw correlation exists is a shared trend; the unstable fraction peaks at 13.27 h and falls 44% by the end of the sigma = 50 km record while the ratio does not follow. What is left to account for: a downward component steady in the ground frame, at about the same level on two runs whose weighting of the affected modes differs by 3.8, not growing, not descending and not tracking the instability, which are the properties of something in the steady solution the solver settles into rather than of a propagating return. Also fixed: a local named `share` in the level-decay loop of section 11 shadowed the ridge flux-share array of section 8, harmless until section 12 read it; renamed `refl`. Verification: additive, the 15 h analysis reproducing all 1760 recorded leaves with none removed and none differing by more than 1.7e-11 relative, only `transient` added. New figure `linear_flux_transient.png`. `docs/linear_flux_constancy.md` section 17 added and the old section 17 renumbered to 18 with its items reordered. Next: the amplitude control, h = 30 m at sigma = 50 km and 15 h, launched the same day, because R is dimensionless and every remaining candidate except one is a linear operator on the wave, so only a source nonlinear in the amplitude can move it; prediction fixed before it finishes, that the steady unreturned-class 0.0097 either holds within the +-0.0017 the bin scatter supports, making the origin linear, or falls toward the 0.0029 closure bound, making it the overturning layer. dx = 2 km moves last, horizontal truncation entering as (k dx)^2 and so requiring the downward spectrum to be 144 times more spread than the upward one across n = 1-12 against a largest measured excess of 1.7. |
| 2026-09-03 | Amplitude control run at h = 30 m; the steady downward component in the modes that cannot have received a sponge reflection is not amplitude-independent, which removes every linear candidate for it, and the floor of the split is corrected upward by a factor of 14 | Motivation: R is dimensionless, and every candidate left after the previous entry except one is a linear operator acting on the wave. Sponge reflection, the lower boundary condition, the residue of the start-up adjustment and truncation error on either axis all scale the numerator and the denominator together and leave R unchanged when the forcing amplitude changes; a source nonlinear in the wave amplitude does not. That makes the amplitude the one variable separating the two classes without a second grid, a second sponge or a longer record. Run: `julia --project=. scripts/sweep_arsia_mons_stage2.jl linear sigma=50 hours=15 h=0.03`, everything else identical, 3684.6 s wall clock over 73055 steps against 3695.1 s and 73266 steps at h = 100 m, which confirms dt is set by the vertical acoustic limit and not by the wave. Analysed with `--h-max 30`. Prediction fixed before it finished: the steady unreturned-class ratio at 9.75 km either holds at 0.0097 within the +-0.0017 the bin scatter supports, making the origin linear, or falls toward the 0.0029 closure bound. Result, writing each quantity as h^p with the incident flux at h^2 so that p = 2 means the ratio is unchanged: the steady unreturned class goes 0.00969 -> 0.00199 at 9.75 km (p = 3.31) and 0.00707 -> 0.00100 at 29.75 km (p = 3.62), both h = 30 m values below the closure bounds of 0.0029 and 0.0038 and so no longer measurable. The controls establish that the comparison has the sensitivity: the fluctuating parts of both branches give p = 1.96 to 2.26, five of six within 0.07 of 2, so the component the previous entry identified as one the partition divides evenly between its branches is amplitude-independent, as that reading requires. The returned class behaves differently, p = 2.14 at 9.75 km, amplitude-independent to 14%, and 2.53 and 2.39 at 29.75 and 54.75 km. Two amplitudes fix a two-parameter reading, an amplitude-independent part plus one quadratic in h, which is what a source quadratic in the wave field gives since its radiated flux would go as h^4 against the incident h^2: at 9.75 km a = 0.00123 with 0.00846 quadratic at h = 100 m, at 29.75 km a = 0.00040 with 0.00667, so 87% and 94% of what was measured at h = 100 m is amplitude-dependent and both offsets lie below the closure bound. Consequence: a downward component whose steady part scales as h^3.3 to h^3.6 cannot be produced by the sponge, the lower boundary, the start-up adjustment or truncation on either axis, which is the whole of the candidate list the previous entry left open, and what remains is a nonlinearity that is not yet located. Supporting numbers: the steady linear reference scales as h^2 exactly (ratio 11.111) while the measured incident flux scales as 10.28, so the shortfall against that reference falls from 15.3% to 8.5% and about half of it is itself amplitude-dependent rather than a numerical loss; band R 0.0486 -> 0.0269; R(z) at 10/30/55 km 0.0371/0.0480/0.0786 -> 0.0257/0.0219/0.0739, so the height dependence inside the band largely disappears while the value just below the sponge barely moves; peak mean-flow deviation -2.150 -> -0.206 m/s against -0.194 predicted as h^2; min u/U 0.363 -> 0.762 against 0.809 predicted by linear scaling, the 6% shortfall being the sense a saturating h = 100 m run gives; the 40-55 km rho versus rho_0 weighting spread, which is the direct measure of the triple correlation, 0.142 -> 0.057 %/km against 0.043 predicted as h. One prediction fails and withdraws an interpretation: the N^2 < 0 area fraction was expected to go to zero and instead rose, 0.00207 -> 0.00750. Resolved in height, those points sit in both runs in a narrow 58.75-61.75 km layer straddling the sponge base at 60 km, seven levels at h = 100 m and five at h = 30 m, with none below, so that diagnostic measures a static instability generated at the top of the damping layer rather than wave overturning in the free column and does not scale with the wave; the previous entry's null result is unaffected but the layer should not be described as the overturning of the mountain wave and it is not the nonlinearity this run detects. Separate finding, the floor of the method. `synthetic_wkb()` builds its field from the same m = sqrt(N^2/U^2 - k^2) the split assumes, so it exercises the algebra, the conjugate-pair weighting and the WKB amplitude gradient but cannot see an error in m, and there is one: the anelastic vertical structure equation for Psi = sqrt(rho_0) w' at constant N and constant scale height carries m^2 = N^2/U^2 - k^2 - 1/(4H^2), and with the band-mean N = 0.00849 1/s and H = 9.61 km that last term is 6.00% of the first, so the assumed m is low by 3.0%. Passing the exact solution through the split returns 2.438e-04 over the band and 2.394e-04 to 2.556e-04 across modes 1-12, matching ((1-r)/(1+r))^2 with r = m_exact/m to five digits, while that field's own flux is height-independent to 0.0000%, which checks the construction. It is a systematic bias of one sign in every mode, so it does not average down, and it rather than the 1.7e-05 is the floor a per-mode ratio should be read against; added as `synthetic_exact()` and reported as `updown/exact_check` next to the WKB check. Nothing reverses: the margins become 84 times the floor for the 0.0202 and 40 times for the 0.0097, with the 0.0029 closure bound binding in either case, and the timing arguments are untouched since they turn on when a mode can have received anything rather than on how small the floor is. Verification: additive, the 15 h analysis reproducing all 1760 previously recorded leaves with none removed and none differing by more than 1.7e-11 relative, only `transient` and `updown/exact_check` added. `docs/linear_flux_constancy.md` section 18 added and the old section 18 renumbered to 19 with its items 1 to 3 rewritten. Next: locate the nonlinearity. A third amplitude at h = 10 m separates the offset-plus-quadratic reading (0.00131 predicted at 9.75 km) from a single exponent R going as h^1.32 (0.00047), a factor of 2.8, both above the 2.4e-04 floor and below the 0.0029 closure bound; going up is not available since min u/U reaches about 0.04 near h = 150 m. Alongside it, a wave-overturning diagnostic restricted to below 55 km reporting min N^2 rather than the fraction below zero, and the triple correlation resolved by horizontal wavenumber, both at no simulation cost. The sponge variation follows, now aimed at the returned class and run on the h = 30 m configuration, since the returned class is the one with the amplitude signature of a linear channel. |
| 2026-09-03 | Taylor-Goldstein theory side checked against the diagnostic; the amplitude reference's factor of 2.08 is identified exactly, partial reflection off the represented N(z) is quantified, and the amplitude-independent offset of the previous entry is shown to be unresolved | A companion notebook `scripts/math_mountain_waves.ipynb` with figures in `output/math_mountain_waves/` carries the Taylor-Goldstein derivations behind sections 15 to 18 and reproduces Nappo Figures 3.11 and 3.12. Two of its results bear on the record and both were checked independently before being quoted; the rest has not been verified. First, the amplitude reference. `linear_steady_flux()` normalises the WKB amplitude as (rho_0 N)^-1/2 at the launch level, so the reference carries N at the first grid point, where the MCD column has a +11.23 K/km inversion against a -4.66 K/km dry adiabat. Evaluating that function with the observed N(z) and with a constant band-mean N gives 2.0759, and N(0)/N_band = 0.017631/0.008493 = 2.0759 to four decimals, so the factor of 2.08 carried since section 9 is entirely the choice of launch N and not the WKB step, which narrows that open item to fixing the choice rather than replacing the method. Second, partial reflection off the profile itself, which is the distributed interior reflector the timing of the previous entries does not exclude. Integrating Phi'' + [N^2/U^2 - k^2 - s^2 + s'] Phi = 0 for Phi = sqrt(rho_0) w' with s = -d ln(rho_0)/dz / 2, downward from the sponge base with a pure upgoing radiation condition and decomposing at the base of the band, gives a flux-weighted reflected share of 4.87e-04 over the three unreturned modes at 9.75 km, with no solver and no sponge in the calculation. Two limits on it. It is not cumulative: at a base of 29.75 km it is 7.89e-04, larger than over the longer path, because contributions from different parts of the column interfere, whereas the fitted offsets run the other way (0.00123 at 9.75 km against 0.00040 at 29.75 km), so the two do not track in height. And it depends on how N(z) is represented rather than only on the data: 4.87e-04 uses the centred-difference N on the 500 m grid, which is the N the split itself forms, while linearly interpolating the same 500 m MCD column onto a finer grid makes dT/dz a staircase and raises the share to 1.84e-03, converged from 250 m to 62.5 m. The subagent's 6.3e-04 is the 500 m value and its convergence claim does not survive that refinement; the 500 m value is nonetheless the right one to compare against the run. Third and consequently, an uncertainty the previous entry did not carry. The amplitude-independent offset is 1.0989 R(30) - 0.0989 R(100) and so is borne almost entirely by the smaller-amplitude run; propagating the bin-to-bin standard errors of the steady part over ten one-hour bins, 0.00101 at h = 30 m and 0.00156 at h = 100 m, gives 0.00123 +- 0.00112, that is 1.1 standard errors from zero. The offset is therefore not resolved on two amplitudes, the 7.31e-04 of computable linear contributions (2.44e-04 partition bias plus 4.87e-04 reflection) sits inside its one standard error so it is consistent with being entirely linear, and the height mismatch above is a reason not to read that as an accounting. Unaffected: the amplitude-dependent part, 0.00846 of 0.00969 at h = 100 m, which is well resolved and is the finding of the previous entry. `docs/linear_flux_constancy.md` section 18 gains a subsection and an error bar on its two-parameter table, and section 19 item 6 records the N(0)/N_band identity. Note that `scripts/` is gitignored, so the notebook and its figures are untracked. No simulation time; these checks run in seconds. |

---

## Arsia Mons Stage 1 sweep results

Analysis of the 9-run Stage 1 sweep (configuration in the 2026-06-30 baseline row above). Per-run figures, `metrics.json`, the comparative `sweep_composite_w.png`, and `amec_response.csv/.md` are in `output/arsia_mons_stage1_analysis/`. Vertical-velocity amplitudes are sampled in the 45-53 km AMEC band at the final time (t ~ 10780 s). Hill metrics restrict to the hill window \|x_rel\| <= 500 km (the composite view); the global peak is over the full 2000 km width and is reported only to document a boundary feature.

| U_bar (m/s) | h_max (km) | Nh/U | hill peak \|w\| (m/s) | hill mean \|w\| (m/s) | late CoV | global peak \|w\| (m/s) | global peak x (km) | wall-clock (s) |
|---|---|---|---|---|---|---|---|---|
| 20 | 0.1 | 0.09 | 0.04 | 0.009 | 0.21 | 0.27 | -998 | 1214 |
| 20 | 1.0 | 0.90 | 0.24 | 0.056 | 0.05 | 0.37 | -998 | 1230 |
| 20 | 3.0 | 2.70 | 0.69 | 0.161 | 0.06 | 0.69 | -306 | 1253 |
| 40 | 0.1 | 0.04 | 0.17 | 0.038 | 0.19 | 13.50 | -946 | 1309 |
| 40 | 1.0 | 0.45 | 0.96 | 0.266 | 0.14 | 13.29 | -886 | 1327 |
| 40 | 3.0 | 1.35 | 3.36 | 0.791 | 0.11 | 13.86 | -890 | 1393 |
| 57 | 0.1 | 0.03 | 0.41 | 0.097 | 0.11 | 52.93 | -910 | 1509 |
| 57 | 1.0 | 0.32 | 2.59 | 0.742 | 0.20 | 50.76 | -910 | 1509 |
| 57 | 3.0 | 0.95 | 49.62 | 3.987 | 0.37 | 49.62 | -178 | 1622 |

The late CoV is the coefficient of variation of the hill-region AMEC mean \|w\| over the final third of the run, used here as a stationarity indicator.

### Observed patterns

- The hill-region AMEC mean \|w\| increases monotonically with both U_bar and h_max, spanning about 0.009 to 4.0 m/s across the grid. At fixed U_bar the dependence on h_max is approximately linear (e.g. U=40: 0.038, 0.266, 0.791 for h=0.1/1/3 km), consistent with the linear flat-bottom forcing w_bottom = U_bar*dh/dx scaling with h_max. The dependence on U_bar is steeper than linear (e.g. h=1 km: 0.056, 0.266, 0.742 for U=20/40/57). This is consistent with the vertical wavelength lambda_z = 2*pi*U/N growing with U_bar, allowing the wave to reach the AMEC band with less attenuation.
- The composite field shows phase lines tilting upstream (toward negative x) with height, the expected structure of a vertically propagating stationary orographic gravity wave. The AMEC-band amplitude maximum accordingly sits about 180-370 km upwind of the summit rather than directly above it.
- Wave activity is concentrated below the 60 km sponge base; the region above 60 km is comparatively quiet in all panels, consistent with sponge absorption rather than reflection into the AMEC band.
- Eight of the nine runs reach a quasi-steady hill-region amplitude after an initial gravity-wave front passes the band (the front peaks near 20 min, then the amplitude relaxes to a lower level), with late-window CoV between 0.05 and 0.21. The U=57/h=3 km corner (Nh/U=0.95) does not settle within 3 h; its hill-region amplitude continues to rise through t_end (late CoV 0.37) and its peak \|w\| of about 50 m/s is comparable to the steering wind. This pattern is consistent with nonlinear steepening near the overturning threshold (\|w\| ~ U). The reported ~50 m/s should be read as a late-time, still-increasing amplitude rather than a converged steady state; the run completed without timestep collapse but is the marginal case of the grid.

### Boundary feature

The unrestricted peak \|w\| is dominated by a persistent feature near the left inflow edge (x ~ -886 to -998 km, about 900 km upwind of the summit), reaching roughly 10-50 m/s. Its amplitude scales with U_bar and is nearly independent of h_max, and it is present from early in each run. This is consistent with a small stationary perturbation at the Dirichlet inflow boundary amplifying as rho^-1/2 with height, a structural limitation of the outflow x-BC under a uniform inflow. The feature is spatially separated from the hill (outside the \|x_rel\| <= 500 km window) and does not contaminate the hill-region AMEC metrics. It is an additional motivation for the Stage 2 periodic-x configuration.

### Figures and data

- `output/arsia_mons_stage1_analysis/sweep_composite_w.png` - 3x3 final-time w*sqrt(rho/rho0), shared color scale, AMEC band shaded, sponge base dashed.
- `output/arsia_mons_stage1_analysis/run_U<U>_h<h>/` - per-run snapshot, AMEC-level (z ~ 50 km) station time series, wavelet, mean-flow evolution, and metrics.json.
- `output/arsia_mons_stage1_analysis/amec_response.csv` and `.md` - source of the response table (hill and global metrics per run).

---

## Arsia Mons Stage 2 plan

Stage 2 launches the Stage 1 orographic forcing into the real reversing zonal wind to test whether a stationary mountain wave deposits momentum near the 45-53 km AMEC band. Code is in `scripts/sweep_arsia_mons_stage2.jl` (runner) and `scripts/analyze_arsia_stage2.py` (analysis). Status 2026-07-13: the 3x3 sweep has been executed (`output/arsia_mons_stage2/`, LX=4000 km, 9/9 converged at t=10800 s, wall-clock 2814-2971 s per run) plus a smoke run (`output/arsia_mons_stage2_smoke/`); the analyzer has been run over it (results section below). The width sweep (`output/arsia_mons_stage2_width/`, LX in {2000,3000,4000} km at U=40/h=1) has also been run and analyzed (width subsection below); hill-region results are width-independent to ~5-12%, so the Stage 2 configuration axis is closed.

### Configuration

- Background wind. Real `u_zonal(z)` from `external/mcd_profiles/mars_arsia_mons_my31.csv` (column 6), installed as a `u_init` Function closure interpolated onto `grid.z_full` and clamped at the endpoints. U(0)=0; the profile is a westward jet over ~4.5-53 km (magnitude up to ~57 m/s near 22 km), an eastward jet over ~53-78 km (up to ~47 m/s near 65 km), then westward again above. Zero-crossings at ~4.5, 53.3, 77.7 km are critical levels for a stationary (c=0) wave; the 53.3 km level coincides with the top of the AMEC band.
- Forcing. Orographic `w = U_bar dh/dx` with `U_bar` held at the Stage 1 eastward values {20, 40, 57} m/s. Because U(0)=0, this `U_bar` is an effective launching amplitude, deliberately decoupled from the background u(z) the wave then propagates through; holding it identical to Stage 1 isolates the effect of the shear. Gaussian hill sigma=100 km, h_max in {0.1, 1, 3} km.
- Horizontal BC. Periodic, which removes the Stage 1 Dirichlet inflow edge whose stationary perturbation amplified as rho^-1/2 with height. The momentum flux is reported as a full-domain zonal mean and the orographic amplitude as a hill-windowed mean of |w|.
- Vertical column. Unchanged from Stage 1 (90 km, dz=500 m, nz=180; sponge 60-90 km, K0=0.05). The 53.3 km critical level lies below the sponge base and stays in the clean region; the 77.7 km level lies inside the sponge and is damped, which is acceptable because it is above the AMEC band.
- Sweeps. The main grid is 3x3 (U_bar x h_max) at a baseline width (default LX=4000 km; override with `lx=NNNN`). Domain width is a separate axis {2000, 3000, 4000} km run at a representative (U=40, h=1 km) via the `width` mode. t_end=10800 s (3 h), as in Stage 1.

### Lessons carried from Stage 1

- Report the orographic response as a hill-windowed mean of |w| rather than an unrestricted domain peak (the Stage 1 peak was contaminated by the inflow-edge artifact).
- Use periodic-x to remove the inflow-edge artifact at its source.
- Track and report per-run wall-clock duration (recorded in `sweep_results.csv`, joined into the response table).

### Expected outcomes

- For a freely propagating stationary wave the momentum-flux profile MF(z) is approximately constant with height. A reduction in |MF| approaching the 53.3 km critical level would be consistent with momentum deposition in or just below the AMEC band; the analyzer reports a deposition fraction 1 - |MF(~52 km)| / |MF(~10 km)|.
- The wave launched in the westward jet is expected to encounter the 53.3 km critical level from below. If absorption is strong, the AMEC-band amplitude and the field above ~53 km may differ from the uniform-wind Stage 1 baseline at matched forcing.
- The lower critical level at ~4.5 km may limit how much of the surface-launched disturbance reaches the jet, so the effective source may be weaker than the nominal U_bar suggests.

### Potential challenges

- Periodic wrap-around. With the reversing wind, low-level disturbances advect westward and mid-level disturbances eastward, and both wrap eventually. At LX=2000 km the nearer x-boundary is ~1000 km from the hill, about 4.9 h of advection at the jet speed, longer than the 3 h integration, so wrap is not expected to reach the hill region within the run; the width axis tests this directly. A vertical sponge does not absorb horizontal advection, so width rather than the sponge is the control.
- Critical-level resolution. Near u(z)=0 the vertical wavelength lambda_z = 2 pi |U| / N shrinks toward zero, so with dz=500 m the wave becomes underresolved just below 53 km and numerical or viscous absorption may substitute for the physical critical-level absorption. The deposition diagnostic should be read as a combined critical-level-and-dissipation signal rather than a pure linear-theory result.
- Effective forcing amplitude. As in Stage 1, the flat-bottom linear BC cannot represent the literal ~20 km obstacle (Nh/U well above 1, blocked regime), so U_bar and h_max remain effective parameters rather than the physical mountain.
- Background maintenance. The strong shear above ~78 km (about -80 m/s at 90 km, steeper aloft) sits in the sponge and is relaxed to the windy background. The inviscid background is an exact discrete steady state and the diffusion preserves the jet, but any slow drift of u(z) should be checked against `u_init` during analysis.

### Outputs (when run)

- `output/arsia_mons_stage2/run_U<U>_h<h>/` and `output/arsia_mons_stage2_width/run_LX<lx>_U<U>_h<h>/` snapshots plus `sweep_results.csv` (adds an `LX_km` column).
- `scripts/analyze_arsia_stage2.py` produces `output/arsia_mons_stage2_analysis/` with per-run figures (w field with critical levels alongside the MF(z) and u(z) profiles), `stage2_composite_w.png`, `stage2_mf_overlay.png`, and `stage2_response.csv` / `.md`. Use `--width-mode` for the width sweep.

---

## Arsia Mons Stage 2 results

Analysis of the 9-run Stage 2 main sweep (real reversing u_zonal(z), periodic-x, LX=4000 km; configuration in the plan section above). Analyzer executed 2026-07-13 over `output/arsia_mons_stage2/` (miniforge base python; analysis wall-clock 4.3 s). Figures and tables in `output/arsia_mons_stage2_analysis/`. All amplitudes are final-time (t ~ 10760-10800 s); hill metrics use the \|x_rel\| <= 500 km window; MF(z) = <rho u' w'>_x is the full-domain zonal mean.

| U_bar (m/s) | h_max (km) | Nh/U | hill peak \|w\| (m/s) | hill mean \|w\| (m/s) | MF@30km (Pa) | MF@56km (Pa) | dep. frac 30->56 km | wall-clock (s) |
|---|---|---|---|---|---|---|---|---|
| 20 | 0.1 | 0.09 | 0.04 | 0.009 | +4.8e-08 | +8.1e-09 | 0.83 | 2814 |
| 20 | 1.0 | 0.90 | 0.37 | 0.081 | +4.7e-06 | +7.5e-07 | 0.84 | 2829 |
| 20 | 3.0 | 2.70 | 1.86 | 0.301 | +3.4e-05 | +5.2e-06 | 0.85 | 2841 |
| 40 | 0.1 | 0.04 | 0.08 | 0.017 | +1.9e-07 | +3.2e-08 | 0.83 | 2819 |
| 40 | 1.0 | 0.45 | 1.06 | 0.184 | +1.8e-05 | +2.7e-06 | 0.85 | 2833 |
| 40 | 3.0 | 1.35 | 3.87 | 0.564 | +8.5e-05 | +1.2e-05 | 0.86 | 2888 |
| 57 | 0.1 | 0.03 | 0.11 | 0.023 | +3.9e-07 | +6.5e-08 | 0.84 | 2819 |
| 57 | 1.0 | 0.32 | 1.76 | 0.287 | +3.2e-05 | +4.8e-06 | 0.85 | 2839 |
| 57 | 3.0 | 0.95 | 16.39 | 0.858 | +8.1e-05 | +1.2e-05 | 0.85 | 2971 |

The deposition fraction is 1 - MF(56 km)/MF(30 km), i.e. the fraction of the jet-level upward flux absorbed across the 45-56 km layer containing the AMEC band and the 53.3 km critical level. The 30 km reference (mid westward jet) replaces the analyzer's built-in 10 km reference, which sits too close to the 4.5 km critical level: below ~5 km most of the launched flux is already absorbed, and MF at 10 km is near zero or sign-flipped in the h >= 1 km runs, making the shipped `mf_deposited_frac` column of `stage2_response.csv` unstable (values from -6.3 to +0.9). The multi-level MF values above were computed from the final snapshots directly.

### Observed patterns

- Wave activity terminates at the 53.3 km critical level in every forced panel of `stage2_composite_w.png`; the region above it is quiet. This contrasts with Stage 1 (uniform wind), where waves filled the column up to the 60 km sponge base. The deposition fraction across 45-56 km is nearly uniform at 0.83-0.86 over the full forcing grid, including the linear h=0.1 km runs, and MF is reduced by more than 98% by 60 km. This uniformity across two decades of forcing amplitude is consistent with critical-level absorption acting as an approximately linear filter at this resolution (with dz=500 m the vertical wavelength is underresolved just below the critical level, so the diagnostic combines critical-level absorption with numerical/viscous dissipation, as anticipated in the plan).
- MF in the jet (10-53 km) is positive in all runs, the sign that decelerates the westward jet toward zero when deposited; the mean-flow response is therefore an eastward acceleration near the AMEC-band top, opposite in sign to the Stage 1 deceleration of the eastward uniform wind.
- Most of the launched flux never reaches the jet: MF drops sharply across the lowest ~5 km (the 4.5 km critical level), with strong small-scale structure trapped below it (visible below 5 km in the per-run w fields). The effective source entering the jet is much weaker than the nominal U_bar suggests, as anticipated in the plan.
- Hill-region AMEC mean \|w\| still increases monotonically with U_bar and h_max but is reduced relative to Stage 1 at matched forcing for U >= 40 (U=57/h=3: 0.86 vs 3.99 m/s; U=40/h=1: 0.18 vs 0.27 m/s), while the weakest-wind runs are comparable or slightly larger (U=20/h=1: 0.08 vs 0.06 m/s). The U=57/h=3 corner no longer shows the Stage 1 runaway growth: its hill peak \|w\| is 16.4 m/s (vs 49.6 m/s in Stage 1) and the final field is confined below 53.3 km.
- The periodic-x BC removed the Stage 1 inflow-edge artifact: the global peak \|w\| equals the hill peak in all nine runs (no boundary feature column is needed in the response table).

### Width sensitivity

Three additional runs at the representative (U=40, h=1 km) with LX in {2000, 3000, 4000} km (`sweep_arsia_mons_stage2.jl width`, executed 2026-07-13; all converged at t=10800 s; wall-clock 770/1080/1498 s, approximately linear in LX; total elapsed 1506 s with the three runs concurrent). Analyzer run with `--width-mode` (12.5 s) into `output/arsia_mons_stage2_width_analysis/`.

| LX (km) | hill RMS w (m/s, 5-55 km) | AMEC hill mean \|w\| (m/s) | MF@30km (Pa) | MF*LX (Pa km) | dep. frac 30->56 km | far/hill RMS ratio | wall-clock (s) |
|---|---|---|---|---|---|---|---|
| 2000 | 0.217 | 0.206 | +3.9e-05 | 78 | 0.90 | 0.28 | 770 |
| 3000 | 0.221 | 0.198 | +2.6e-05 | 77 | 0.82 | 0.28 | 1080 |
| 4000 | 0.211 | 0.184 | +1.8e-05 | 70 | 0.85 | 0.21 | 1498 |

- Hill-window metrics are width-independent within ~5% (RMS w) to ~11% (AMEC mean \|w\|, weakly decreasing with LX); the composite panels are visually indistinguishable in the \|x_rel\| <= 500 km window. The deposition fraction spread (0.82-0.90) brackets the main-sweep value of ~0.85.
- MF as a full-domain zonal mean dilutes as 1/LX (the wave occupies a fixed physical area near the hill), so the raw MF values are not comparable across widths; the width-integrated flux MF*LX agrees to ~10%. Cross-width MF comparisons should use MF*LX or restrict the average to a fixed window.
- No wrap-around contamination is detected: the far-field (\|x_rel\| > 700 km, 5-55 km) RMS w is 3-5x below the hill-window value at all widths, and the LX=2000 far field matches LX=3000. This is consistent with the plan's advection estimate (nearest boundary ~4.9 h away at the jet speed, longer than the 3 h integration).
- Practical consequence: the LX=4000 km baseline is width-robust, and LX=2000 km reproduces the hill-region physics within ~10% at half the wall-clock cost, making it a reasonable grid for longer integrations (e.g. the planned 6 h extension).

### Figures and data

- `output/arsia_mons_stage2_analysis/stage2_composite_w.png` - 3x3 final-time density-scaled w, critical levels dashed, AMEC band shaded, sponge base dotted.
- `output/arsia_mons_stage2_analysis/stage2_run_U<U>_h<h>.png` - per-run w field with the u(z) and MF(z) profiles alongside.
- `output/arsia_mons_stage2_analysis/stage2_mf_overlay.png` - all nine MF(z) profiles overlaid (shows the dominant low-level absorption).
- `output/arsia_mons_stage2_analysis/stage2_response.csv` / `.md` - response table (note the `mf_deposited_frac` caveat above).
- `output/arsia_mons_stage2_width_analysis/` - width-mode products (1x3 composite, MF overlay, response table); same `mf_deposited_frac` caveat applies.

---

## Arsia Mons Stage 1 U=40/h=3 time extension

A 6 h Stage 1 integration of the U_bar = 40 m/s, h_max = 3 km mode (Nh/U = 1.35), run to test whether the momentum budget approaches a quasi-steady state beyond the 3 h sweep t_end. Configuration is identical to the Stage 1 sweep (uniform eastward u_init = 40 m/s, outflow x-BC, LX = 2000 km, 90 km domain with the sponge over 60-90 km, MY31 MCD background, no radiative cooling) except for t_end. Launched with `julia -t auto --project=. scripts/sweep_arsia_mons_stage1.jl ext u=40 h=3 hours=6`; 32253 steps to t = 21600 s, wall-clock 2434.1 s.

### Why this mode was extended

The deep-dive next-steps list proposed extending U=57/h=3. The target was changed to U=40/h=3 because that column produces the cleaner wave field at h = 1 and 3 km. Baselines were already available in both stages, so no new 3 h runs were required:

| baseline | final t | AMEC hill mean \|w\| | MF at 30 km | last-hour normalized trend, MF at 30 km | last-hour normalized trend, AMEC \|w\| |
|---|---|---|---|---|---|
| `output/arsia_mons_stage1/run_U40.0_h3.0` | 10753 s | 0.79 m/s | -2.16 mPa | -0.79 per hour | +0.16 per hour |
| `output/arsia_mons_stage2/run_U40.0_h3.0` | 10784 s | 0.56 m/s | +0.085 mPa | +0.63 per hour | +0.54 per hour |

Both were still evolving at t_end, so the motivation that applied to U=57/h=3 applies here as well. The normalized trend is the least-squares slope over the window divided by the window mean magnitude, so values of order unity indicate a quantity changing on the run's own timescale. The Stage 2 MF at 30 km carries the opposite sign to Stage 1 because the real zonal wind reverses within the column, so the two MF values are not directly comparable in magnitude; the trends are the relevant comparison here.

### Response between 3 h and 6 h

Hill-window metrics use \|x_rel\| <= 500 km; all vertical metrics are restricted to z < 60 km, below the sponge base. The AMEC band is 45-53 km. Static instability is evaluated as N^2 = (g/theta) dtheta/dz < 0 with theta = T (P_surface/P)^kappa, using the solver's own Mars constants (gamma = 1.31, g = 3.72 m/s^2, R_gas = 188.9 J/(kg K), so kappa = 0.2366) rather than assumed values; the reported fraction is the fraction of hill-window grid points below 60 km satisfying that condition.

| quantity | t = 3 h | t = 6 h |
|---|---|---|
| AMEC-band hill mean \|w\| | 0.79 m/s | 4.06 m/s |
| hill-window max \|w\|, z < 60 km | 6.9 m/s | 51.1 m/s |
| height of that maximum | 59.75 km | 53.75 km |
| peak <u>_x - U_bar | -8.6 m/s at 59.75 km | -26.1 m/s at 55.25 km |
| AMEC-band mean <u>_x | 40.8 m/s | 22.2 m/s |
| minimum <u>_x below 60 km | 31.3 m/s | 13.9 m/s |
| minimum hill-window u | -21.8 m/s | -64.2 m/s |
| static-instability fraction (N^2 < 0, hill window, z < 60 km) | 0.004 | 0.086 |
| MF at 29.8 km | -2.16 mPa | -4.04 mPa |
| MF at 44.8 km | -0.78 mPa | -1.33 mPa |
| MF at 51.8 km | -0.53 mPa | -0.64 mPa |

### Observed patterns

- The extension reproduces the 3 h baseline over the shared interval. Sampled quantities agree to within ~1% at matched times (for example the 2-3 h window mean of the AMEC amplitude, 0.750 against 0.748 m/s, and of MF at 30 km, -1.51 against -1.50 mPa), so the two records can be read as one continuous integration.
- The momentum flux approaches a plateau after about 4 h. Over the last hour the normalized trends are +0.16, +0.16 and +0.29 per hour at 30, 45 and 52 km with a coefficient of variation of 0.06-0.16, against -0.79, -0.69 and -0.67 per hour over the 2-3 h window. The sign change indicates the flux magnitudes stop growing and fluctuate about a level rather than continuing to increase.
- The mean flow does not plateau over the same interval. The peak deviation deepens from -8.6 to -26.1 m/s, and the rate falls only from about -4.6 m/s per hour over 2-3 h to about -3.2 m/s per hour over 5-6 h. The decelerated layer also descends. Taking its base as the lowest level at which the mean flow has slowed by more than 5 m/s, it moves from 57.3 km at 3 h to 53.8 km at 4 h, 50.8 km at 4.5 h, 41.3 km at 5 h and 35.3 km at 6 h, so the AMEC band lies within it from about 5 h onward. The metric excludes z < 2 km, which is necessary: the 0.75 km level carries a persistent deviation reaching -5.7 m/s by 6 h that is a separate near-surface feature rather than part of the descending layer, and including it collapses the reported base to 0.75 km from 4 h onward. A weaker 1 m/s threshold reaches lower, near 23 km at 6 h, but at that threshold the descent is not monotonic and the layer base is not cleanly separable from the same near-surface signal.
- The AMEC-band amplitude continues to grow, and its normalized trend increases rather than decreases (+0.25 per hour over the last hour against +0.16 per hour over 2-3 h). The band amplitude and the local mean wind are coupled here: as <u>_x in the band falls from 40.8 to 22.2 m/s the local Nh/U rises, which is consistent with the amplitude growth continuing while the flux is saturating.
- Overturning indicators grow monotonically through the record. The static-instability fraction rises from 0.000 at 1 h to 0.004 at 3 h, 0.023 at 4 h, 0.053 at 5 h and 0.086 at 6 h, concentrated near 53-57 km, and the hill-window minimum u becomes increasingly negative. The final hill-window max \|w\| of 51.1 m/s exceeds U_bar.
- Taken together, the data do not support a quasi-steady interpretation of the 6 h state. They are more consistent with the mode entering a saturated, breaking regime in which the flux is amplitude-limited while the deposited momentum continues to decelerate the mean flow, and in which the wave is progressing toward a self-induced critical level where <u>_x approaches zero.
- The 6 h U=40/h=3 state resembles the 3 h U=57/h=3 state in these indicators (that run had a static-instability fraction of 0.063, a hill-window minimum u of -103 m/s and a hill-window max \|w\| of 57.3 m/s at 3 h). The U=40 case carries the larger Nh/U (1.35 against 0.95) but develops more slowly, which is consistent with both the advection time and the vertical group velocity scaling with U.
- The self-induced deceleration in Stage 1 develops near 53-55 km, close to the 53.3 km critical level that the real reversing wind imposes in Stage 2. Whether that coincidence is more than incidental has not been tested.

### Caveats

- The breaking layer at 53-57 km sits 3-7 km below the sponge base at 60 km, and the height of the peak deceleration moves from 59.75 km at 3 h to 55.25 km at 6 h. An influence of the sponge on the late-time deceleration maximum cannot be excluded from this run alone.
- At dx = 4 km and dz = 500 m the overturning eddies are not resolved. The post-4 h state should be read as an indication of the onset of a breaking regime rather than as a quantitatively resolved breaking solution.
- Radiative cooling is disabled and molecular viscosity is capped, so no thermal damping limits the amplitude growth. The outstanding Eckermann et al. (2011) CO2 radiative damping item bears directly on this, since scale-dependent damping would be expected to reduce the growth rate above the AMEC band.
- Only the Stage 1 configuration was extended. The Stage 2 U=40/h=3 run was also shown to be trending at its 3 h t_end and has not been extended.

### Figures and data

- `output/arsia_mons_stage1_ext6h/run_U40.0_h3.0/` - 323 snapshots at 74.0 s cadence, 1.17 GB; `sweep_results.csv` alongside carries the wall-clock duration.
- `output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_convergence.png` - AMEC-band amplitude, MF at three levels and the peak mean-flow deviation against time, with the 3 h baseline overlaid dashed and its t_end marked.
- `.../deepdive_trendfit.png` - the same five series with the least-squares fits over the 2-3 h and 5-6 h windows drawn over their windows and annotated with slope, normalized trend and CoV. This is the figure that supports the trend values quoted above.
- `.../deepdive_hovmoller.png` - time-height evolution of the mean-flow deviation, momentum flux and hill-window RMS w over the full 6 h below 60 km, with contours of the deceleration depth. This is the continuous form of the six-profile figure and supports the descent statement.
- `.../deepdive_breaking.png` - the full 323-snapshot record of the unstable area fraction, minimum hill-window u, maximum hill-window \|w\| against U_bar and minimum zonal-mean u, plus a time-height map of where instability appears. This replaces an earlier six-point hand sampling.
- `.../deepdive_series.npz` - cached time-height arrays (`ubar_zt`, `mf_zt`, `unstable_frac_zt`, `wrms_hill_zt` and the derived scalar series) written by the deep-dive script, so the figure scripts do not re-read the snapshots.
- `output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_mf_ubar_profiles.png` - MF(z) and <u>_x(z) - U_bar at 0.5, 1, 2, 3, 4.5 and 6 h, showing the descending decelerated layer.
- `output/arsia_mons_stage1_ext6h_analysis/run_U40.0_h3.0/deepdive_w_profiles.png`, `deepdive_station_timeseries_w_z57.png`, `deepdive_wavelet_w_full_z50.png`, `deepdive_wavelet_w_full_z57.png`, `deepdive_metrics.json`.
- `scripts/analyze_ext6h_mechanism.py` - mechanism attribution against analytical results
  (WKB growth, overturning onset, momentum-budget closure, Lindzen saturation level and
  kinematic descent speed, deposition partition, sponge and resolution tests). Reads the same
  `deepdive_series.npz` plus the MCD column; 1.5 s wall-clock; writes
  `mechanism_metrics.json` alongside. Results and interpretation in
  `docs/ext6h_instability_mechanism.md`.
- Reproduce with `python3 scripts/analyze_arsia_stage1_deepdive.py --run-dir output/arsia_mons_stage1_ext6h/run_U40.0_h3.0 --compare-run output/arsia_mons_stage1/run_U40.0_h3.0` (193.3 s wall-clock, dominated by reading 478 snapshots). That pass also writes `deepdive_series.npz`.
- The three evidence figures are then produced from that cache by `scripts/viz_ext6h_trendfit.py`, `scripts/viz_ext6h_hovmoller.py` and `scripts/viz_ext6h_breaking.py`, which take `--npz` and `--out` and run in about 5 s for all three combined. They are kept as separate per-figure scripts to match the existing one-script-per-product convention in `scripts/`.
