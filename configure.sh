#!/usr/bin/env bash


# author: Alejandro Gonzales-Irribarren  
# email: jose.gonzalesdezavala1@unmsm.edu.pe
# github: https://github.com/alejandrogzi


GET_BINARIES=true
INSTALL_PYTHON=false

if $INSTALL_PYTHON && [[ -f "./modules/requirements.py" ]]; then
    echo "Installing Python dependencies..."
    ./modules/requirements.py
fi


if command -v cargo &> /dev/null; then
    if $GET_BINARIES; then
        echo "Installing Rust binaries..."
        cat ./modules/rust-binaries.txt | xargs -n 1 cargo install
    else
        if [[ ! -f "./bed2gtf/Cargo.toml" ]] || [[ ! -f "./bed2gff/Cargo.toml" ]] || [[ ! -f "./noel/Cargo.toml" ]]; then
            git submodule init bed2gtf bed2gff noel && git submodule update bed2gtf bed2gff noel
        fi
            echo "Building bed2gtf..."
            cd ./bed2gtf && cargo build --release && cargo install --path .

            echo "Building bed2gff..."
            cd ../bed2gff && cargo build --release && cargo install --path .

            echo "Building noel..."
            cd ../noel && cargo build --release && cargo install --path .
    fi
else
    echo "cargo not found, please install Rust"
fi
