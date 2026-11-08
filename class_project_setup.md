# ELMFIRE class project setup

A personal copy of ELMFIRE and the 2018 Camp Fire case in `$HOME`, with every
run and all output on `$SCRATCH`, readable by classmates.

| What | Where |
|---|---|
| ELMFIRE (code and prebuilt binaries) | `~/elmfire` |
| Camp Fire case | `~/elmfire-class/camp_fire` |
| Job script and run instructions (`README.md`) | `~/elmfire-class/` |
| Run directories and Slurm logs | `$SCRATCH/elmfire-class/` |

Sources, all copied (nothing is moved):

| Copy | From |
|---|---|
| `~/elmfire` | `$GROUP_HOME/codes/elmfire` (commit `cbf924a`, ELMFIRE 1.1, built with gcc 12.4.0 + OpenMPI 5.0.5) |
| `~/elmfire-class/README.md`, `camp_fire.slurm` | `$OAK/jgoell/Wildfire/elmfire-class/` |
| `~/elmfire-class/camp_fire` | `$OAK/jgoell/ELMFIRE_VnV_Suite/cases/Validation/landscape_scale/camp_fire` |

## 1. Copy everything

Run from `sh_dev` (about 790 MB in total: 703 MB of ELMFIRE, 364 MB of which
is git history, and 84 MB for the case).

```bash
sh_dev

# ELMFIRE
rsync -a $GROUP_HOME/codes/elmfire/ $HOME/elmfire/
cd $HOME/elmfire/build/linux/bin
for l in $(find . -maxdepth 1 -type l); do ln -sfn "$(basename "$(readlink "$l")")" "$l"; done
ls -l          # links should read e.g. elmfire -> elmfire_1.1

# Project folder (README + job script), then the Camp Fire case into it
rsync -a $OAK/jgoell/Wildfire/elmfire-class/ $HOME/elmfire-class/
SRC=$OAK/jgoell/ELMFIRE_VnV_Suite/cases/Validation/landscape_scale/camp_fire/
rsync -a --exclude={/data/fuels_and_topography,/data/weather,/outputs,/figures,/report,/run_case.sh,/compile_case.sh} $SRC $HOME/elmfire-class/camp_fire/

mkdir -p $SCRATCH/elmfire-class
du -sh $HOME/elmfire $HOME/elmfire-class
exit
```

Why the extra steps:

- The unversioned names in `build/linux/bin` (`elmfire`, `elmfire_post`, ...)
  are symlinks with absolute paths back to the group copy; the loop repoints
  them at the files in your copy. The versioned binaries (`elmfire_1.1`, ...)
  are real files and work as-is with the same modules loaded.
- Left out of the case copy:
  - `data/fuels_and_topography` and `data/weather`: not used by anything in
    the case, and on a different grid (UTM 11N);
  - `outputs`, `figures`: results of the earlier run;
  - `report`, `run_case.sh`, `compile_case.sh`: the suite's LaTeX report,
    which doesn't build on Sherlock (missing font). The job script replaces
    `run_case.sh`.
- `scripts/preprocess.py` is still copied but won't run here: it checksums the
  original `Camp_VAL.rar`, which isn't copied. It isn't needed to run the case.

## 2. Running the case

See `~/elmfire-class/README.md`. In short:

```bash
cd $SCRATCH/elmfire-class
sbatch ~/elmfire-class/camp_fire.slurm
```

## 3. Sharing with classmates

Home and scratch are closed to everyone but the owner by default. Each
classmate needs a Sherlock account. Run once per person, replacing `SUNETID`
(syntax from https://www.sherlock.stanford.edu/docs/storage/data-sharing/;
home uses NFSv4 ACLs, scratch uses POSIX ACLs).

```bash
# Home: pass through (no listing), read + run the code, case and scripts
nfs4_setfacl -a A::SUNETID@sherlock:X $HOME
nfs4_setfacl -R -a A:fd:SUNETID@sherlock:RX $HOME/elmfire $HOME/elmfire-class

# Scratch: pass through, read the project folder and anything created in it later
setfacl -m u:SUNETID:x $SCRATCH
setfacl -R -m u:SUNETID:rX $SCRATCH/elmfire-class
setfacl -R -d -m u:SUNETID:rX $SCRATCH/elmfire-class

# Check (the home filesystem can silently downgrade an entry)
nfs4_getfacl $HOME/elmfire
getfacl $SCRATCH/elmfire-class
```

Use `rwX` instead of `rX` on the scratch lines if classmates should be able
to write into the project folder.

A classmate runs the case from your copy like this; the job copies the case
into their own `$SCRATCH`, so read access is enough:

```bash
export ELMFIRE=~jgoell/elmfire
export CASE=~jgoell/elmfire-class/camp_fire
export VENV=$GROUP_HOME/venvs/elmfire-vnv     # their own venv; see README
mkdir -p $SCRATCH/elmfire-class && cd $SCRATCH/elmfire-class
sbatch ~jgoell/elmfire-class/camp_fire.slurm
```

Limits for classmates:

- The default venv (`$OAK/jgoell/venvs/elmfire-vnv`) is on gorle's Oak space,
  which they can't read. They need their own; `README.md` has the commands.
- `serc` is the SERC partition. Classmates outside it should submit with
  `sbatch -p normal ...`.

## Scratch purge

Files on `$SCRATCH` that are not modified for 90 days are deleted. Copy
anything needed for the final report to `$HOME` or Oak before then.
