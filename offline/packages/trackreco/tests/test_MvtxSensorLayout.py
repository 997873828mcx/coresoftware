#!/usr/bin/env python3
"""Compile exact production MVTX layout helpers without ROOT/ACTS dependencies."""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def main():
    here = Path(__file__).resolve().parent
    source = (here.parent / "MakeActsGeometry.cc").read_text()
    begin = "// MVTX_SENSOR_LAYOUT_HELPERS_BEGIN"
    end = "// MVTX_SENSOR_LAYOUT_HELPERS_END"
    if source.count(begin) != 1 or source.count(end) != 1:
        raise RuntimeError("expected one production MVTX layout helper block")
    helpers = source.split(begin, 1)[1].split(end, 1)[0]
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    with tempfile.TemporaryDirectory(prefix="mvtx-sensor-layout-") as directory:
        work = Path(directory)
        (work / "production_helpers.inc").write_text(helpers)
        executable = work / "test_MvtxSensorLayout"
        subprocess.run(compiler + ["-std=c++20", "-O2", "-Wall", "-Wextra",
                       "-Wshadow", "-Werror", "-I", str(work),
                       str(here / "test_MvtxSensorLayout.cc"), "-o", str(executable)],
                       check=True)
        subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    main()
