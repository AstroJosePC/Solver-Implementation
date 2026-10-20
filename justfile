# ==============================================================================
# MAGNUS-P Solver command runner (justfile)
# Front end of scripts/run_dt_convergence_S2T95.sh that queues whole batches in the
# Task Spooler queue shared with MarsAtmos, so DSMC jobs and solver batches take
# turns on the 16 performance cores.
# ==============================================================================

set shell := ["zsh", "-cu"]

# Same socket as MarsAtmos/justfile: one queue, one slot, for every heavy job
export TS_SOCKET := env_var('HOME') + "/.dsmc_ts.sock"

ts_bin := `command -v ts || command -v tsp || echo "ts"`
repo := justfile_directory()
launcher := "scripts/run_dt_convergence_S2T95.sh"
py := env_var_or_default("PYTHON", env_var('HOME') + "/.venvs/magnusp/bin/python")

# Runs of set <set> write to $MAGNUSP_RUNS/<set> (logs in its logs/), as d1 and d2 do
runs := env_var_or_default("MAGNUSP_RUNS", env_var('HOME') + "/MAGNUSP_runs")

default:
    @just --list

# ------------------------------------------------------------------------------
# 1. Launcher modes, run in the foreground
# ------------------------------------------------------------------------------

# Resolve the cases, commands and directories of a set
dry set:
    @S14_REPO="{{repo}}" MAGNUSP_OUTPUT="{{runs}}/{{set}}" PYTHON="{{py}}" bash {{launcher}} dry {{set}}

# Reproduce the reference records bit for bit and write the launch marker
check set:
    @S14_REPO="{{repo}}" MAGNUSP_OUTPUT="{{runs}}/{{set}}" PYTHON="{{py}}" bash {{launcher}} check {{set}}

# Simulated time and end stamp of every case of a set
status set:
    @S14_REPO="{{repo}}" MAGNUSP_OUTPUT="{{runs}}/{{set}}" PYTHON="{{py}}" bash "{{runs}}/{{set}}/logs/{{set}}.launcher.sh" status {{set}}

# Copy a finished set to the archive with checksums (prints the rm command, runs none)
archive set:
    @S14_REPO="{{repo}}" MAGNUSP_OUTPUT="{{runs}}/{{set}}" PYTHON="{{py}}" bash "{{runs}}/{{set}}/logs/{{set}}.launcher.sh" archive {{set}}

# ------------------------------------------------------------------------------
# 2. Batches
# ------------------------------------------------------------------------------

# Run a whole set in the foreground: guards, optional check, launch, monitor until every case ends
batch set check="1" freeze="1":
    #!/usr/bin/env zsh
    set -eu -o pipefail
    cd "{{repo}}"
    logs="{{runs}}/{{set}}/logs"; frozen="$logs/{{set}}.launcher.sh"
    [ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "refused: tracked changes in the tree"; exit 1; }
    git tag --points-at HEAD | grep -q "^run/{{set}}-" || { echo "refused: no run/{{set}}-<YYYYMMDD> tag on HEAD $(git rev-parse --short HEAD)"; exit 1; }
    mkdir -p "$logs"
    if [ "{{freeze}}" = 1 ]; then
        [ ! -e "$frozen" ] || { echo "refused: $frozen exists (set already launched?)"; exit 1; }
        cp {{launcher}} "$frozen"
    fi
    cmp -s {{launcher}} "$frozen" || { echo "refused: $frozen differs from {{launcher}} at HEAD"; exit 1; }
    export S14_REPO="{{repo}}" MAGNUSP_OUTPUT="{{runs}}/{{set}}" PYTHON="{{py}}"
    if [ "{{check}}" = 1 ]; then bash {{launcher}} check {{set}} | tee "$logs/{{set}}.check.out"; fi
    bash "$frozen" launch {{set}} | tee "$logs/{{set}}.launch.out"
    bash "$frozen" monitor {{set}} > "$logs/{{set}}.monitor.out" 2>&1
    bash "$frozen" status {{set}}

# Queue a whole set behind the jobs already in the shared spooler (freezes the launcher now)
queue set check="1":
    #!/usr/bin/env zsh
    set -eu
    cd "{{repo}}"
    logs="{{runs}}/{{set}}/logs"; frozen="$logs/{{set}}.launcher.sh"
    git tag --points-at HEAD | grep -q "^run/{{set}}-" || { echo "refused: no run/{{set}}-<YYYYMMDD> tag on HEAD $(git rev-parse --short HEAD)"; exit 1; }
    [ ! -e "$frozen" ] || { echo "refused: $frozen exists (set already launched?)"; exit 1; }
    mkdir -p "$logs"; cp {{launcher}} "$frozen"
    {{ts_bin}} -L "MAGNUSP_{{set}}" "{{just_executable()}}" --justfile "{{repo}}/justfile" batch {{set}} {{check}} 0
    echo "[+] set {{set}} queued at $(git rev-parse --short HEAD) (socket $TS_SOCKET); the job refuses to start if HEAD moves"

# Hold the queue while a set launched outside the spooler is running (waits on its monitor)
hold set:
    {{ts_bin}} -L "MAGNUSP_{{set}}_hold" zsh -c 'while pgrep -f "{{set}}.launcher.sh monitor {{set}}" > /dev/null; do sleep 60; done'
    @echo "[+] queue held until the monitor of set {{set}} exits"

# ------------------------------------------------------------------------------
# 3. Shared spooler (also holds the MarsAtmos DSMC jobs)
# ------------------------------------------------------------------------------

# Active and queued jobs of the shared spooler
q:
    @{{ts_bin}}

# Follow the output of a job (Ctrl+C leaves the job running)
tail id="":
    @if [ -z "{{id}}" ]; then {{ts_bin}} -t; else {{ts_bin}} -t {{id}}; fi

# Full output of a finished or running job
log id:
    @{{ts_bin}} -c {{id}}

# Remove a queued job that has not started
remove id:
    @{{ts_bin}} -r {{id}}
