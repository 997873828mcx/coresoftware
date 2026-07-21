#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
g4tpc_dir="$(cd "${script_dir}/.." && pwd)"
build_dir="${TMPDIR:-/tmp}/tpc_helix_looper_test.$$"
mkdir -p "${build_dir}"
trap 'rm -rf "${build_dir}"' EXIT

cxx="${CXX:-g++}"
eigen_include="${EIGEN_INCLUDE:-/usr/include/eigen3}"

"${cxx}" -std=c++17 -O2 \
  -I"${g4tpc_dir}" \
  -I"${eigen_include}" \
  "${script_dir}/tpc_helix_looper_test.cc" \
  "${g4tpc_dir}/TpcTrackHelixFitter.cc" \
  -o "${build_dir}/tpc_helix_looper_test"

"${build_dir}/tpc_helix_looper_test"
