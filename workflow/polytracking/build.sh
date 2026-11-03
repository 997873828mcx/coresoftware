#!/usr/bin/env bash
set -eo pipefail
if [[ $# != 1 || ${1:-} == --help ]]; then
  echo "Usage: bash $0 NEW_INSTALL_PREFIX" >&2
  echo "Optional: TPC_BUILD_JOBS=2; TPC_BUILD_RESUME=1 to continue the same source build." >&2
  exit 2
fi
tpc_share_source=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
tpc_share_prefix=$(realpath -m -- "$1")
tpc_share_build="${tpc_share_prefix}-build"
tpc_share_revision=$(git -C "$tpc_share_source" rev-parse HEAD)
if [[ -e $tpc_share_prefix || -e $tpc_share_build ]]; then
  if [[ ${TPC_BUILD_RESUME:-0} != 1 || ! -r $tpc_share_build/source-revision ||
        $(cat "$tpc_share_build/source-revision") != "$tpc_share_revision" ]]; then
    echo "Choose a new prefix, or set TPC_BUILD_RESUME=1 for the same source revision." >&2
    exit 2
  fi
fi
case "$tpc_share_prefix/" in
  "$tpc_share_source/"*) echo "Choose an installation prefix outside the source checkout." >&2; exit 2 ;;
esac
mkdir -p "$tpc_share_prefix" "$tpc_share_build"
printf '%s\n' "$tpc_share_revision" > "$tpc_share_build/source-revision"
unset MYINSTALL LD_LIBRARY_PATH ROOT_INCLUDE_PATH PYTHONPATH LD_PRELOAD ROOTSYS OFFLINE_MAIN
export USER="$(id -un)"
export LOGNAME="$USER"
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.576
export LD_LIBRARY_PATH="$tpc_share_prefix/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$tpc_share_prefix/include:${ROOT_INCLUDE_PATH:-}"
export PATH="$tpc_share_prefix/bin:$PATH"
export CPPFLAGS="-I$tpc_share_prefix/include"
export CXXFLAGS="$(root-config --cflags) -O2 -g0 -std=c++20"
export LDFLAGS="-L$tpc_share_prefix/lib -L$tpc_share_prefix/lib64"
tpc_share_jobs=${TPC_BUILD_JOBS:-2}
[[ $tpc_share_jobs =~ ^[1-9][0-9]*$ ]] || { echo "TPC_BUILD_JOBS must be positive." >&2; exit 2; }
for tpc_share_package in \
  offline/framework/ffamodules \
  offline/packages/trackbase_historic \
  simulation/g4simulation/g4decayer \
  simulation/g4simulation/g4main \
  offline/packages/PHGarfield \
  offline/packages/tpc \
  offline/packages/tpccalib \
  offline/packages/trackreco \
  offline/packages/tpctrackreco \
  simulation/g4simulation/g4tpc \
  offline/packages/TrackingDiagnostics
do
  tpc_share_name=${tpc_share_package##*/}
  if [[ -f $tpc_share_build/$tpc_share_name.complete ]]; then
    echo "Reusing completed $tpc_share_name"
    continue
  fi
  echo "Building $tpc_share_name"
  mkdir -p "$tpc_share_build/$tpc_share_name"
  (
    cd "$tpc_share_build/$tpc_share_name"
    "$tpc_share_source/$tpc_share_package/autogen.sh" --prefix="$tpc_share_prefix"
    make -j"$tpc_share_jobs"
    make -j"$tpc_share_jobs" check
    make install
  ) > "$tpc_share_build/$tpc_share_name.log" 2>&1 || {
    echo "Build/check failed; inspect $tpc_share_build/$tpc_share_name.log" >&2
    exit 1
  }
  touch "$tpc_share_build/$tpc_share_name.complete"
done
echo "Installed: $tpc_share_prefix"
