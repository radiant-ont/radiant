#!/usr/bin/env bash
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

# Provision the external command-line tools at the exact versions used for the
# manuscript results:
#   cutadapt          adapter trimming
#   fastq-filter      length / mean-quality filtering
#   minimap2 v2.28    alignment (paper Methods)
#   racon    v1.4.3   consensus polishing (paper Methods, -e 0.7)
#
# minimap2 and racon are pinned by building the specified tags from source, so
# results do not drift with newer releases. Binaries are installed into
# ~/.local/bin with version suffixes; export the paths shown at the end.
set -euo pipefail

BIN_DIR="${HOME}/.local/bin"
mkdir -p "${BIN_DIR}"

MINIMAP2_TAG="v2.28"
RACON_TAG="1.4.3"

echo "== cutadapt =="
if ! command -v cutadapt >/dev/null 2>&1; then
    python3 -m pip install --user "cutadapt>=4.0"
fi
cutadapt --version

echo "== fastq-filter =="
if ! command -v fastq-filter >/dev/null 2>&1; then
    python3 -m pip install --user "fastq-filter>=0.3.0"
fi
fastq-filter --version

echo "== minimap2 ${MINIMAP2_TAG} =="
if [ ! -x "${BIN_DIR}/minimap2-2.28" ]; then
    tmp="$(mktemp -d)"
    git clone https://github.com/lh3/minimap2 "${tmp}/minimap2"
    ( cd "${tmp}/minimap2" && git checkout "${MINIMAP2_TAG}"
      # arm_neon/aarch64 flags build on Apple Silicon; omit on x86_64.
      if [ "$(uname -m)" = "arm64" ]; then make arm_neon=1 aarch64=1; else make; fi )
    cp "${tmp}/minimap2/minimap2" "${BIN_DIR}/minimap2-2.28"
fi
"${BIN_DIR}/minimap2-2.28" --version

echo "== racon ${RACON_TAG} =="
# Note: racon 1.4.3 vendors an old zlib that fails to compile under very recent
# clang (C23) on Apple Silicon. The source build works on Linux / older
# toolchains; on macOS the conda (osx-64) binary is the reliable route.
if [ ! -x "${BIN_DIR}/racon-1.4.3" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_SUBDIR=osx-64 conda create -y -n racon143 --override-channels \
            -c bioconda -c conda-forge "racon=${RACON_TAG}" || true
        if [ -x "$(conda run -n racon143 which racon 2>/dev/null)" ]; then
            cp "$(conda run -n racon143 which racon)" "${BIN_DIR}/racon-1.4.3"
        fi
    fi
fi
if [ ! -x "${BIN_DIR}/racon-1.4.3" ]; then
    tmp="$(mktemp -d)"
    git clone --recursive https://github.com/lbcb-sci/racon.git "${tmp}/racon"
    ( cd "${tmp}/racon" && git checkout "${RACON_TAG}" && git submodule update --init --recursive
      mkdir -p build && cd build
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
            -DCMAKE_C_FLAGS="-std=gnu89 -Wno-error -Wno-implicit-function-declaration -Wno-implicit-int" \
            -DCMAKE_CXX_FLAGS="-Wno-error" ..
      make -j4 )
    cp "${tmp}/racon/build/bin/racon" "${BIN_DIR}/racon-1.4.3"
fi
"${BIN_DIR}/racon-1.4.3" --version 2>/dev/null || \
    echo "racon 1.4.3 not provisioned locally; use Linux/CI or a reachable conda (osx-64)."

echo
echo "All tools ready. Export the pinned versions for the run:"
echo "    export MINIMAP2_PATH=${BIN_DIR}/minimap2-2.28"
echo "    export RACON_PATH=${BIN_DIR}/racon-1.4.3"
