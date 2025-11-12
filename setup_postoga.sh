#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

VENV_DIR=".venv"

if ! command -v uv >/dev/null 2>&1; then
    echo "Error: uv is required but not found in PATH." >&2
    exit 1
fi

if ! command -v maturin >/dev/null 2>&1; then
    echo "Error: maturin is required but not found in PATH." >&2
    exit 1
fi

echo "Creating virtual environment with uv..."
uv venv

if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
    echo "Error: expected ${VENV_DIR}/bin/activate but it was not created." >&2
    exit 1
fi

echo "Installing Python dependencies via uv pip..."
source "${VENV_DIR}/bin/activate"
uv pip install "."

echo "Building rustools extension with maturin..."
(
    cd rustools
    maturin develop --release
)

echo "Environment ready. Activate it anytime with 'source ${VENV_DIR}/bin/activate' and run postoga normally."
