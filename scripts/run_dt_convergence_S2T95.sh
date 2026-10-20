#!/usr/bin/env bash
# =============================================================================
# Section 14 suite of docs/reflection_convergence_tests.md: time-step convergence of
# the Stage 2 flux on the 95 km top, and the launch-layer hypothesis.
#
#   scripts/run_dt_convergence_S2T95.sh dry      resolved configuration of every case
#   scripts/run_dt_convergence_S2T95.sh check    0.1 h of S04 and C2R, bit for bit
#   scripts/run_dt_convergence_S2T95.sh launch   Z5 blend gate, then all cases in parallel
#   scripts/run_dt_convergence_S2T95.sh status   simulated time reached by each case
#   scripts/run_dt_convergence_S2T95.sh archive  copy to output/, verify, print the rm
#   scripts/run_dt_convergence_S2T95.sh monitor  per-step cost and CPU clocks every 60 s; keeps the
#                                                machine from idle sleep while cases run
#
# A second argument `g` selects the ground-forcing cases of the Section 14 addendum
# (S14G1-G3, C2R configuration; stamps in S14G.start), `a2` the cases of addendum 2
# (S14M1-M2 MCD face pair, S14K1 split sponge, S14X1 z-then-x order, S14R1 dz = 250 m,
# S14P1-P3 ground-face probes; stamps in S14A2.start), and `a3` the cases of addendum 3
# (V1a-b face-path ridge F/F_lin, V2u/V2m H = 0 equilibrium, V3g/V3f Queney ghost and
# face, S1 30 h stationarity, S2 MCD dz = 250 face, S3E-F launch-consistent face pair,
# S4 MCD z-then-x; stamps in S14A3.start), and `s15` the cases of Section 15 (S15F1 and
# S15F2, the h = 1 km dz pair on the face path, S15F1 doubling as PR1-U40-h1.0; the
# other eight PR1 face-path 24 h production cases; S15T1, the CFL 0.2 twin of
# PR1-U40-h0.1; S14V3p, the Queney face path with pbar=1; stamps in S15.start), and `s16`
# the cases of Section 16 (S16<scheme><level>: the 2 h face-path ridge under the Lie
# splitting L, Strang sweeps and operators S, alternating sweeps with Strang operators A
# and Strang sweeps with Lie operators X, each at CFL 0.4, 0.2, 0.1 and 0.05 for levels
# 1 to 4; stamps in S16.start) for dry, launch, status, archive and monitor; `check` is
# shared. Logs are named by run ID
# (case_id): S14<case> for the earlier sets. At launch the set records the git commit,
# the dirty flag and `git diff` (<stamps>.git.diff), writes a provenance header at the
# top of each log, and each case appends one row to docs/run_matrix.csv when it ends. `monitor` appends to $LOGS/<stamps>.cost.tsv until every case of
# the set has an end stamp.
#
# Bash reads a script as it runs, so a launch runs a copy of this file kept in $LOGS:
#   cp scripts/run_dt_convergence_S2T95.sh $LOGS/S14A3.launcher.sh
#   S14_REPO=$PWD bash $LOGS/S14A3.launcher.sh launch a3   (and monitor a3)
# S14_REPO points the copy at the repository; without it the script uses its own parent.
#
# Runs write to $MAGNUSP_OUTPUT (default ~/MAGNUSP_runs), one process per case with
# julia -t 1 and OPENBLAS_NUM_THREADS=1. Logs: $MAGNUSP_OUTPUT/logs/S14<id>.log, start
# and end stamps in $MAGNUSP_OUTPUT/logs/S14.start. `launch` refuses to start unless
# `check` has passed with the current driver (marker file holds the driver's hash).
# =============================================================================
set -u
cd "${S14_REPO:-$(dirname "$0")/..}"
export OPENBLAS_NUM_THREADS=1
export MAGNUSP_OUTPUT="${MAGNUSP_OUTPUT:-$HOME/MAGNUSP_runs}"
DRIVER=scripts/sweep_arsia_mons_stage2.jl
QUENEY=scripts/run_mountain_benchmark_strict.jl
ANALYZE=scripts/analyze_dt_convergence.py
PY="${PYTHON:-python3}"
LOGS="$MAGNUSP_OUTPUT/logs"
STAMPS="$LOGS/S14.start"
MARKER="$LOGS/S14.check_ok"
ARCHIVE="${S14_ARCHIVE:-output}"   # S14_ARCHIVE: the archive of a checkout run from a worktree
SET="${2:-s14}"
CASES="A B C D E F Z"
[ "$SET" = g ] && { CASES="G1 G2 G3"; STAMPS="$LOGS/S14G.start"; }
[ "$SET" = a2 ] && { CASES="M1 M2 K1 X1 R1 P1 P2 P3"; STAMPS="$LOGS/S14A2.start"; }
[ "$SET" = a3 ] && { CASES="S2 V3g V3f V1b S3F S1 V1a S3E S4 V2m V2u"; STAMPS="$LOGS/S14A3.start"; }
[ "$SET" = s15 ] && { CASES="F2 T1 F1 h01u20 h01u40 h01u57 h1u20 h1u57 h3u20 h3u40 h3u57 V3p"
                      STAMPS="$LOGS/S15.start"; }
[ "$SET" = s16 ] && { CASES="s16S4 s16X4 s16A4 s16L4 s16S3 s16X3 s16A3 s16L3 s16S2 s16X2 s16A2 s16L2 s16S1 s16X1 s16A1 s16L1"
                      STAMPS="$LOGS/S16.start"; }
LIN="linear mode=6 h=0.005 hours=18"
COMMON="pbar=1 ramp=2 k0=0.005 lz=95"
PROBE="linear mode=6 h=0.005 hours=1 pbar=1 k0=0.005 lz=95"
RIDGE="linear sigma=50 h=0.03 hours=18"
PR24="hours=24 $COMMON bcface=1"
SO="linear sigma=50 h=0.03 hours=2 $COMMON bcface=1 final=1"

# Section 16: the scheme (fourth character of the case) and the CFL level (fifth) of
# an s16 case, as driver arguments and as directory tags.
s16_scheme() {
    case "${1:3:1}" in
        L) echo "" ;; S) echo "strang=1 opstrang=1" ;; A) echo "alt=1 opstrang=1" ;; X) echo "strang=1" ;;
    esac
}
s16_stag() {
    case "${1:3:1}" in
        L) echo "" ;; S) echo "_strang_opS" ;; A) echo "_alt_opS" ;; X) echo "_strang" ;;
    esac
}
s16_level() {
    case "${1:4:1}" in
        1) echo "" ;; 2) echo "cfl=0.2 outint=200" ;; 3) echo "cfl=0.1 outint=400" ;; 4) echo "cfl=0.05 outint=800" ;;
    esac
}
s16_ltag() {
    case "${1:4:1}" in
        1) echo "" ;; 2) echo "_cfl0.2" ;; 3) echo "_cfl0.1" ;; 4) echo "_cfl0.05" ;;
    esac
}
RUN_MATRIX=docs/run_matrix.csv
QOUT="$MAGNUSP_OUTPUT/mountain_benchmark_strict_10h_dx250"
mkdir -p "$LOGS"

# Blend height: the smallest of 8, 10 and 12 km that passes the Z5 gate (Section 14).
# Only the background column enters the rule. Empty if none passes.
ublend_height() {
    local z
    for z in 8 10 12; do
        if "$PY" "$ANALYZE" zblend --zb-km "$z" | grep -q '"gate_pass": true'; then
            echo "$z"; return
        fi
    done
}

# Driver arguments of each case (after the driver path); UB is the blend height.
case_args() {
    case "$1" in
        A) echo "hours=9 $COMMON u=40 h=0.1 cfl=0.2 outint=200" ;;
        B) echo "hours=9 $COMMON u=40 h=0.1 cfl=0.1 outint=400" ;;
        C) echo "linear mode=6 h=0.005 hours=18 $COMMON cfl=0.2 outint=200" ;;
        D) echo "hours=9 $COMMON u=40 h=0.1 dz=250 outint=200" ;;
        E) echo "hours=9 $COMMON u=-40 h=0.1 ublend=$UB" ;;
        F) echo "hours=9 $COMMON u=-40 h=0.1 ublend=$UB cfl=0.2 outint=200" ;;
        Z) echo "hours=1 $COMMON u=40 h=0 dz=250 outint=200" ;;
        G1) echo "$LIN $COMMON bcface=1" ;;
        G2) echo "$LIN $COMMON bcface=1 cfl=0.2 outint=200" ;;
        G3) echo "$LIN $COMMON cfl=0.1 outint=400" ;;
        M1) echo "hours=9 $COMMON u=40 h=0.1 bcface=1" ;;
        M2) echo "hours=9 $COMMON u=40 h=0.1 bcface=1 cfl=0.2 outint=200" ;;
        K1) echo "$LIN $COMMON bcface=1 spsplit=1" ;;
        X1) echo "$LIN $COMMON bcface=1 zx=1" ;;
        R1) echo "$LIN $COMMON bcface=1 dz=250 outint=200" ;;
        P1) echo "$PROBE probe=25" ;;
        P2) echo "$PROBE cfl=0.2 probe=50" ;;
        P3) echo "$PROBE cfl=0.1 probe=100" ;;
        V1a) echo "$RIDGE $COMMON bcface=1" ;;
        V1b) echo "$RIDGE $COMMON bcface=1 cfl=0.2 outint=200" ;;
        V2u) echo "linear mode=6 h=0 hours=1 pbar=1 bcface=1" ;;
        V2m) echo "hours=1 $COMMON u=40 h=0 bcface=1" ;;
        V3g) echo "--dx=250 --tend_hr=10 --outdir=$QOUT/run_U18_b10" ;;
        V3f) echo "--dx=250 --tend_hr=10 --bcface=1 --outdir=${QOUT}_bcface/run_U18_b10" ;;
        S1) echo "linear mode=6 h=0.005 hours=30 $COMMON bcface=1" ;;
        S2) echo "hours=9 $COMMON u=40 h=0.1 bcface=1 dz=250 outint=200" ;;
        S3E) echo "hours=9 $COMMON u=-40 h=0.1 ublend=8 bcface=1" ;;
        S3F) echo "hours=9 $COMMON u=-40 h=0.1 ublend=8 bcface=1 cfl=0.2 outint=200" ;;
        S4) echo "hours=9 $COMMON u=40 h=0.1 bcface=1 zx=1" ;;
        F1) echo "$PR24 u=40 h=1" ;;
        F2) echo "$PR24 u=40 h=1 dz=250 outint=200" ;;
        T1) echo "$PR24 u=40 h=0.1 cfl=0.2 outint=200" ;;
        h01u20) echo "$PR24 u=20 h=0.1" ;;
        h01u40) echo "$PR24 u=40 h=0.1" ;;
        h01u57) echo "$PR24 u=57 h=0.1" ;;
        h1u20) echo "$PR24 u=20 h=1" ;;
        h1u57) echo "$PR24 u=57 h=1" ;;
        h3u20) echo "$PR24 u=20 h=3" ;;
        h3u40) echo "$PR24 u=40 h=3" ;;
        h3u57) echo "$PR24 u=57 h=3" ;;
        V3p) echo "--dx=250 --tend_hr=10 --bcface=1 --pbar=1 --outdir=${QOUT}_bcface_pbarP/run_U18_b10" ;;
        s16*) echo "$SO $(s16_scheme "$1") $(s16_level "$1")" ;;
    esac
}

# Driver of each case: the Queney cases (V3) use the Srivastava strict benchmark driver.
case_driver() {
    case "$1" in
        V3g|V3f|V3p) echo "$QUENEY" ;;
        *) echo "$DRIVER" ;;
    esac
}

# Sweep directory name and run subdirectory of each case.
case_dir() {
    case "$1" in
        A) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2/run_U40.0_h0.1" ;;
        B) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_cfl0.1/run_U40.0_h0.1" ;;
        C) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_cfl0.2/run_U40.0_h0.005" ;;
        D) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250/run_U40.0_h0.1" ;;
        E) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_ub$UB/run_U-40.0_h0.1" ;;
        F) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_ub${UB}_cfl0.2/run_U-40.0_h0.1" ;;
        Z) echo "arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h_dz250/run_U40.0_h0.0" ;;
        G1) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.005" ;;
        G2) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface_cfl0.2/run_U40.0_h0.005" ;;
        G3) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_cfl0.1/run_U40.0_h0.005" ;;
        M1) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1" ;;
        M2) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface_cfl0.2/run_U40.0_h0.1" ;;
        K1) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface_spsplit/run_U40.0_h0.005" ;;
        X1) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface_zx/run_U40.0_h0.005" ;;
        R1) echo "arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h0.005" ;;
        P1) echo "arsia_mons_linear1h_LX2000_mode6_K0_0.005_LZ95_pbarP/run_U40.0_h0.005" ;;
        P2) echo "arsia_mons_linear1h_LX2000_mode6_K0_0.005_LZ95_pbarP_cfl0.2/run_U40.0_h0.005" ;;
        P3) echo "arsia_mons_linear1h_LX2000_mode6_K0_0.005_LZ95_pbarP_cfl0.1/run_U40.0_h0.005" ;;
        V1a) echo "arsia_mons_linear18h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.03" ;;
        V1b) echo "arsia_mons_linear18h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface_cfl0.2/run_U40.0_h0.03" ;;
        V2u) echo "arsia_mons_linear1h_LX2000_mode6_pbarP_bcface/run_U40.0_h0.0" ;;
        V2m) echo "arsia_mons_stage2_1h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.0" ;;
        V3g) echo "mountain_benchmark_strict_10h_dx250/run_U18_b10" ;;
        V3f) echo "mountain_benchmark_strict_10h_dx250_bcface/run_U18_b10" ;;
        S1) echo "arsia_mons_linear30h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.005" ;;
        S2) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h0.1" ;;
        S3E) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_ub8_bcface/run_U-40.0_h0.1" ;;
        S3F) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_ub8_bcface_cfl0.2/run_U-40.0_h0.1" ;;
        S4) echo "arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface_zx/run_U40.0_h0.1" ;;
        F1) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h1.0" ;;
        F2) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_dz250_bcface/run_U40.0_h1.0" ;;
        T1) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface_cfl0.2/run_U40.0_h0.1" ;;
        h01u20) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U20.0_h0.1" ;;
        h01u40) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1" ;;
        h01u57) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U57.0_h0.1" ;;
        h1u20) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U20.0_h1.0" ;;
        h1u57) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U57.0_h1.0" ;;
        h3u20) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U20.0_h3.0" ;;
        h3u40) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h3.0" ;;
        h3u57) echo "arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U57.0_h3.0" ;;
        V3p) echo "mountain_benchmark_strict_10h_dx250_bcface_pbarP/run_U18_b10" ;;
        s16*) echo "arsia_mons_linear2h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface$(s16_stag "$1")$(s16_ltag "$1")/run_U40.0_h0.03" ;;
    esac
}

# Run ID of each case, which names its log and its stamp and run-matrix rows.
case_id() {
    case "$1" in
        F1|F2|T1) echo "S15$1" ;;
        V3p) echo "S14V3p" ;;
        s16*) echo "S16${1#s16}" ;;
        h01u*) echo "PR1-U${1#h01u}-h0.1" ;;
        h1u*) echo "PR1-U${1#h1u}-h1.0" ;;
        h3u*) echo "PR1-U${1#h3u}-h3.0" ;;
        *) echo "S14$1" ;;
    esac
}

driver_hash() { shasum -a 256 "$DRIVER" "$QUENEY" src/*.jl | shasum -a 256 | cut -c1-16; }

mode="${1:-}"
case "$mode" in
dry)
    UB="$(ublend_height)"; UB="${UB:-8}"
    for c in $CASES; do
        echo "### $(case_id "$c"): $(case_driver "$c") $(case_args "$c")"
        [ "$(case_driver "$c")" = "$DRIVER" ] || continue   # the Queney driver has no dryrun
        # shellcheck disable=SC2046
        julia -t 1 --project=. "$DRIVER" dryrun $(case_args "$c") 2>&1 |
            grep -E "vertical|grid|ublend|widths|pbar_mode|output|ERROR"
    done
    ;;
check)
    # 0.1 h of the S04 and C2R commands, written to a scratch root and compared with the
    # archived records snapshot by snapshot. Every earlier record uses the defaults of
    # dx, dz, outint and ublend, so this checks that the patched driver leaves them intact.
    # Two further runs check the face path against S14G1 and that the ground-face probe
    # leaves the C2R solution unchanged.
    CK="$MAGNUSP_OUTPUT/s14_check"
    rm -rf "$CK"; mkdir -p "$CK/probe" "$CK/queney_new" "$CK/queney_ref" "$CK/queney_face" "$CK/queney_pbar"
    MAGNUSP_OUTPUT="$CK" julia -t 1 --project=. "$DRIVER" \
        hours=0.1 $COMMON u=40 h=0.1 > "$LOGS/S14check_S04.log" 2>&1 &
    MAGNUSP_OUTPUT="$CK" julia -t 1 --project=. "$DRIVER" \
        linear mode=6 h=0.005 hours=0.1 $COMMON > "$LOGS/S14check_C2R.log" 2>&1 &
    MAGNUSP_OUTPUT="$CK" julia -t 1 --project=. "$DRIVER" \
        linear mode=6 h=0.005 hours=0.1 $COMMON bcface=1 > "$LOGS/S14check_G1.log" 2>&1 &
    MAGNUSP_OUTPUT="$CK/probe" julia -t 1 --project=. "$DRIVER" \
        linear mode=6 h=0.005 hours=0.1 $COMMON probe=5 > "$LOGS/S14check_probe.log" 2>&1 &
    # Addendum 3: 0.1 h of the face-path MCD run S14M1, and the Queney driver with its
    # default path against the copy kept before the bcface option was added (dx = 1000 m,
    # 0.7 h, snapshots at steps 0 and 2000), plus a run of its face path that must complete.
    MAGNUSP_OUTPUT="$CK" julia -t 1 --project=. "$DRIVER" \
        hours=0.1 $COMMON u=40 h=0.1 bcface=1 > "$LOGS/S14check_M1.log" 2>&1 &
    julia -t 1 --project=. "$QUENEY" --dx=1000 --tend_hr=0.7 --outdir="$CK/queney_new" \
        > "$LOGS/S14check_queney_new.log" 2>&1 &
    julia -t 1 --project=. "$LOGS/a3_ref/scripts/run_mountain_benchmark_strict.pre_a3.jl" --dx=1000 --tend_hr=0.7 \
        --outdir="$CK/queney_ref" > "$LOGS/S14check_queney_ref.log" 2>&1 &
    julia -t 1 --project=. "$QUENEY" --dx=1000 --tend_hr=0.7 --bcface=1 --outdir="$CK/queney_face" \
        > "$LOGS/S14check_queney_face.log" 2>&1 &
    # Section 15: the Queney face path with pbar=1 (S14V3p) must complete on :pressure.
    julia -t 1 --project=. "$QUENEY" --dx=1000 --tend_hr=0.7 --bcface=1 --pbar=1 --outdir="$CK/queney_pbar" \
        > "$LOGS/S14check_queney_pbar.log" 2>&1 &
    # Section 16: 0.1 h of the face-path ridge S14V1a on the default splitting, and the
    # Strang and alternating paths with Strang operators, which must reach t_end and
    # write its snapshot (final=1).
    mkdir -p "$CK/s16"
    MAGNUSP_OUTPUT="$CK" julia -t 1 --project=. "$DRIVER" \
        linear sigma=50 h=0.03 hours=0.1 $COMMON bcface=1 > "$LOGS/S14check_V1a.log" 2>&1 &
    MAGNUSP_OUTPUT="$CK/s16" julia -t 1 --project=. "$DRIVER" \
        linear sigma=50 h=0.03 hours=0.1 $COMMON bcface=1 final=1 strang=1 opstrang=1 > "$LOGS/S16check_S.log" 2>&1 &
    MAGNUSP_OUTPUT="$CK/s16" julia -t 1 --project=. "$DRIVER" \
        linear sigma=50 h=0.03 hours=0.1 $COMMON bcface=1 final=1 alt=1 opstrang=1 > "$LOGS/S16check_A.log" 2>&1 &
    wait
    ok=1
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/arsia_mons_linear0.1h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.005" \
        --ref "$ARCHIVE/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.005" || ok=0
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/probe/arsia_mons_linear0.1h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005" \
        --ref "$ARCHIVE/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005" || ok=0
    [ -f "$CK/probe/arsia_mons_linear0.1h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005/face_probe.nc" ] ||
        { echo "probe file missing"; ok=0; }
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/arsia_mons_stage2_0.1h_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.1" \
        --ref "$ARCHIVE/arsia_mons_stage2_24h_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.1" || ok=0
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/arsia_mons_linear0.1h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005" \
        --ref "$ARCHIVE/arsia_mons_linear18h_LX2000_mode6_K0_0.005_LZ95_pbarP_ramp2h/run_U40.0_h0.005" || ok=0
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/arsia_mons_stage2_0.1h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1" \
        --ref "$ARCHIVE/arsia_mons_stage2_9h_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.1" || ok=0
    "$PY" "$ANALYZE" bitcheck --new "$CK/queney_new" --ref "$CK/queney_ref" || ok=0
    [ -f "$CK/queney_face/snapshot_0002000.nc" ] && grep -q 'ground w: face' "$LOGS/S14check_queney_face.log" ||
        { echo "Queney face-path run did not complete"; ok=0; }
    [ -f "$CK/queney_pbar/snapshot_0002000.nc" ] && grep -q 'pbar_mode=pressure' "$LOGS/S14check_queney_pbar.log" ||
        { echo "Queney pbar=1 face-path run did not complete"; ok=0; }
    "$PY" "$ANALYZE" bitcheck \
        --new "$CK/arsia_mons_linear0.1h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.03" \
        --ref "$ARCHIVE/arsia_mons_linear18h_LX2000_sig50_K0_0.005_LZ95_pbarP_ramp2h_bcface/run_U40.0_h0.03" || ok=0
    for s in S A; do
        grep -q 't=360/360 s steps=[0-9]* OK' "$LOGS/S16check_$s.log" ||
            { echo "Section 16 scheme $s did not reach t_end"; ok=0; }
    done
    if [ "$ok" = 1 ]; then
        driver_hash > "$MARKER"; echo "check PASSED (driver hash $(cat "$MARKER"))"
    else
        rm -f "$MARKER"; echo "check FAILED: launch is blocked"; exit 1
    fi
    ;;
launch)
    [ -f "$MARKER" ] && [ "$(cat "$MARKER")" = "$(driver_hash)" ] ||
        { echo "run 'check' first (missing or stale $MARKER)"; exit 1; }
    UB="$(ublend_height)"
    run_cases="$CASES"
    [ "$SET" = g ] || [ "$SET" = a2 ] || [ "$SET" = a3 ] || [ "$SET" = s15 ] && UB=n/a
    if [ -z "$UB" ]; then
        echo "Z5 gate fails at 8, 10 and 12 km: S14E and S14F are withdrawn (Section 14)"
        run_cases="A B C D Z"; UB=none
    fi
    GIT_COMMIT=$(git rev-parse --short=12 HEAD 2>/dev/null || echo none)
    GIT_DIRTY=$([ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ] && echo dirty || echo clean)
    git diff > "${STAMPS%.start}.git.diff" 2>/dev/null
    DH=$(driver_hash); LH=$(shasum -a 256 "$0" | cut -c1-16)
    NCONC=$(echo $run_cases | wc -w | tr -d ' ')
    echo "$(date '+%F %T') S14 launch set $SET (ublend=$UB km) driver $DH launcher $LH ($0) git $GIT_COMMIT $GIT_DIRTY concurrency $NCONC" >> "$STAMPS"
    [ -f "$RUN_MATRIX" ] || echo "run_id,set,start,end,exit,real_s,steps,t_end_s,mean_dt_s,s_per_step,concurrency,git_commit,git_dirty,driver_src_hash,launcher_hash,command,local_dir,archive_dir" > "$RUN_MATRIX"
    for c in $run_cases; do
        mkdir -p "$MAGNUSP_OUTPUT/$(case_dir "$c")"   # the Queney driver does not create its outdir
        (
            id=$(case_id "$c"); log="$LOGS/$id.log"; s=$(date +%s)
            cmd="julia -t 1 --project=. $(case_driver "$c") $(case_args "$c")"
            {
                echo "# run_id $id | set $SET | start $(date '+%F %T')"
                echo "# command: OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS $cmd"
                echo "# git $GIT_COMMIT ($GIT_DIRTY; diff in ${STAMPS%.start}.git.diff) | driver+src hash $DH | launcher $LH ($0)"
                echo "# julia threads 1 | BLAS threads $OPENBLAS_NUM_THREADS | concurrency $NCONC cases"
                echo "# local $MAGNUSP_OUTPUT/$(case_dir "$c") | archive $ARCHIVE/$(case_dir "$c")"
            } > "$log"
            # shellcheck disable=SC2046
            /usr/bin/time -p julia -t 1 --project=. "$(case_driver "$c")" $(case_args "$c") >> "$log" 2>&1
            rc=$?; e=$(date +%s)
            echo "$id exit=$rc start=$s end=$e" >> "$STAMPS"
            real=$(awk '/^real /{print $2}' "$log" | tail -1)
            read -r steps tend <<< "$(sed -nE 's/^Simulation complete: ([0-9]+) steps, t=([0-9.]+) s.*/\1 \2/p' "$log" | tail -1)"
            mdt=$( [ -n "${steps:-}" ] && echo "scale=6; $tend / $steps" | bc || echo NA)
            sps=$( [ -n "${steps:-}" ] && [ -n "$real" ] && echo "scale=4; $real / $steps" | bc || echo NA)
            echo "$id,$SET,$(date -r "$s" '+%F %T'),$(date -r "$e" '+%F %T'),$rc,${real:-NA},${steps:-NA},${tend:-NA},$mdt,$sps,$NCONC,$GIT_COMMIT,$GIT_DIRTY,$DH,$LH,\"$cmd\",$MAGNUSP_OUTPUT/$(case_dir "$c"),$ARCHIVE/$(case_dir "$c")" >> "$RUN_MATRIX"
        ) &
        sleep 2
    done
    echo "launched: $run_cases (ublend=$UB km); stamps in $STAMPS"
    ;;
status)
    UB="$(grep -o 'ublend=[0-9a-z]*' "$STAMPS" 2>/dev/null | tail -1 | cut -d= -f2)"; UB="${UB:-8}"
    for c in $CASES; do
        d="$MAGNUSP_OUTPUT/$(case_dir "$c")"
        last=$(ls "$d"/snapshot_*.nc 2>/dev/null | tail -1)
        if [ -z "$last" ]; then echo "$(case_id "$c")  no snapshots"; continue; fi
        t=$(ncdump -h "$last" | grep -o 'time_s = [0-9.e+]*' | awk '{print $3}')
        end=$(grep "^$(case_id "$c") " "$STAMPS" 2>/dev/null | tail -1)
        printf "%s  t = %8.1f s (%5.2f h)  %s\n" "$(case_id "$c")" "$t" "$(echo "$t/3600" | bc -l)" "${end:-running}"
    done
    ;;
archive)
    # Copies each sweep directory to output/, verifies it by checksum and prints the
    # command that removes the local copy; it does not remove anything itself.
    UB="$(grep -o 'ublend=[0-9a-z]*' "$STAMPS" 2>/dev/null | tail -1 | cut -d= -f2)"; UB="${UB:-8}"
    rmlist=""; done_sweeps=" "
    for c in $CASES; do
        sweep="$(dirname "$(case_dir "$c")")"
        case "$done_sweeps" in *" $sweep "*) continue ;; esac   # cases that share a sweep
        done_sweeps="$done_sweeps$sweep "
        [ -d "$MAGNUSP_OUTPUT/$sweep" ] || { echo "$(case_id "$c"): $sweep not found"; continue; }
        rsync -a "$MAGNUSP_OUTPUT/$sweep" "$ARCHIVE/"
        diffs=$(rsync -a --dry-run --checksum --itemize-changes "$MAGNUSP_OUTPUT/$sweep/" "$ARCHIVE/$sweep/" | wc -l)
        n1=$(find "$MAGNUSP_OUTPUT/$sweep" -type f | wc -l); n2=$(find "$ARCHIVE/$sweep" -type f | wc -l)
        echo "$(case_id "$c"): $sweep files $n1/$n2, checksum differences $diffs"
        [ "$diffs" -eq 0 ] && [ "$n1" -eq "$n2" ] && rmlist="$rmlist \"$MAGNUSP_OUTPUT/$sweep\""
    done
    mkdir -p "$ARCHIVE/reflection_tests/logs"
    if [ "$SET" = s15 ]; then
        for c in $CASES; do cp "$LOGS/$(case_id "$c").log" "$ARCHIVE/reflection_tests/logs/"; done
        cp "${STAMPS%.start}".* "$LOGS"/S14check_* "$ARCHIVE/reflection_tests/logs/"
    else
        cp "$LOGS"/S14* "$ARCHIVE/reflection_tests/logs/"
    fi
    echo "logs copied to $ARCHIVE/reflection_tests/logs/"
    [ -n "$rmlist" ] && echo "verified; remove the local copies with:" && echo "rm -rf$rmlist"
    ;;
monitor)
    # One row per minute: P- and E-cluster clock [MHz] and usage, CPU power [W] and
    # temperature from macmon, then for each case the last logged step and the wall
    # seconds per step since the previous row (logs are flushed every 500 steps).
    # The machine idle-sleeps after 1 min (pmset sleep 1), and the runs then advance only
    # in DarkWake windows at a reduced clock (Section 14 addendum 2). -i alone is not
    # honoured once the machine is in DarkWake; -s (PreventSystemSleep, on AC) is. The
    # monitor holds both for as long as it runs, i.e. until every case has ended.
    caffeinate -s -i -w $$ &
    COST="${STAMPS%.start}.cost.tsv"
    [ -f "$COST" ] || printf "time\tP_MHz\tP_use\tE_MHz\tcpu_W\ttemp_C%s\n" \
        "$(for c in $CASES; do printf '\t%s_step\t%s_s_per_step' "$(case_id "$c")" "$(case_id "$c")"; done)" > "$COST"
    while :; do
        now=$(date +%s)
        hw=$(macmon pipe -s 1 -i 1000 2>/dev/null | "$PY" -c 'import json,sys
d=json.loads(sys.stdin.readline())
print("%d\t%.2f\t%d\t%.2f\t%.1f" % (d["pcpu_usage"][0], d["pcpu_usage"][1], d["ecpu_usage"][0], d["cpu_power"], d["temp"]["cpu_temp_avg"]))')
        row="$(date '+%F %T')\t$hw"; running=0
        for c in $CASES; do
            st=$(grep -a '^Step' "$LOGS/$(case_id "$c").log" 2>/dev/null | tail -1 | awk '{print $2}')
            st="${st:-0}"; cost=NA
            eval "ps=\${prev_step_$c:-}; pt=\${prev_t_$c:-}"
            if [ -n "$ps" ] && [ "$st" -gt "$ps" ]; then
                cost=$(echo "scale=4; ($now - $pt) / ($st - $ps)" | bc)
            fi
            if [ -z "$ps" ] || [ "$st" -gt "$ps" ]; then
                eval "prev_step_$c=$st; prev_t_$c=$now"
            fi
            grep -q "^$(case_id "$c") " "$STAMPS" 2>/dev/null || running=1
            row="$row\t$st\t$cost"
        done
        printf "%b\n" "$row" >> "$COST"
        [ "$running" = 1 ] || break
        sleep 60
    done
    echo "monitor finished; $COST"
    ;;
*)
    sed -n '2,30p' "$0"; exit 1 ;;
esac
