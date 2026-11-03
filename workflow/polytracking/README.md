# Build the shared polytracking simulation sources

This branch contains the core source snapshot for the shared Fun4All pion
simulation and polytracking workflow. It uses frozen release `ana.576` on
an sPHENIX AlmaLinux 9 host with CVMFS.

Build into a new installation prefix, outside the source checkout:

```bash
bash workflow/polytracking/build.sh "$PWD/../install-polytracking-ana576"
```

The script builds the required packages in dependency order, installs their
headers, libraries and ROOT dictionaries, and runs their package checks.
It uses two build processes by default; set `TPC_BUILD_JOBS` to change this.
An interrupted build can be continued by repeating the same command with
`TPC_BUILD_RESUME=1`.

The matching macros and launchers are in
[the macros source branch](https://github.com/997873828mcx/macros/tree/share/polytracking-simulation-ana576-20261102/detectors/sPHENIX/polytracking_simulation_share).
They provide single-sample and Condor entrypoints. To use libraries rebuilt
from this branch, set `POLYTRACKING_CORE_INSTALL` to the new installation
prefix before running those entrypoints.

The source includes the drift/LUT implementation, 55 ns readout response,
diffusion and ADC baseline handling, polytracking and truth evaluation used
by the shared workflow. The complete package sources are retained so their
libraries and dictionaries can be rebuilt consistently. A rebuild can have
different binary hashes because compiler options and build paths affect the
resulting files.

LUTs, field/gas conditions and generated data remain external inputs. The
macros branch records the matching inputs and uses the existing shared LUT.
Changing the physical transport configuration requires corresponding inputs.
