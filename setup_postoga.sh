#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

if ! command -v uv >/dev/null 2>&1; then
    echo "Error: uv is required but not found in PATH." >&2
    exit 1
fi

if ! command -v maturin >/dev/null 2>&1; then
    echo "Error: maturin is required but not found in PATH." >&2
    exit 1
fi

echo "Installing postoga with uv..."
uv pip install "."

echo "Building rustools extension with maturin..."
(
    cd rustools
    maturin develop --release
)

echo "Installation complete. You can now run postoga directly."
