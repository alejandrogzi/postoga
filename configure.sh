#!/usr/bin/env bash


# author: Alejandro Gonzales-Irribarren  
# email: jose.gonzalesdezavala1@unmsm.edu.pe
# github: https://github.com/alejandrogzi


GET_BINARIES=true
INSTALL_PYTHON=true
CONFIG_FILE="$HOME/.bashrc"


add_to_config() {
    local line="$1"
    local file="$2"
    if ! grep -qF "$line" "$file"; then
        echo "$line" >> "$file"
    fi
}


if $INSTALL_PYTHON && [[ -f "./modules/requirements.py" ]]; then
    echo "Installing Python dependencies..."
    ./modules/requirements.py
fi


if command -v cargo &> /dev/null; then
    if $GET_BINARIES; then
        add_to_config 'export PATH="$HOME/.cargo/bin:$PATH"' "$CONFIG_FILE"

        echo "Installing Rust binaries..."
        cat ./modules/rust-binaries.txt | xargs -n 1 cargo install
    else
        if [[ -f "./bed2gtf/Cargo.toml" ]] && [[ -f "./bed2gff/Cargo.toml" ]]; then
            echo "Building bed2gtf..."
            cd ./bed2gtf && cargo build --release && cargo install --path .

            echo "Building bed2gff..."
            cd ../bed2gff && cargo build --release && cargo install --path .
        else
            echo "bed2gtf and bed2gff not found"
            git submodule init bed2gtf bed2gff && git submodule update bed2gtf bed2gff
            echo "Building bed2gtf..."
            cd ./bed2gtf && cargo build --release && cargo install --path .
            cd ../bed2gff && cargo build --release && cargo install --path .
        fi
    fi
else
    echo "cargo not found, installing Rust...\n"
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env

    add_to_config 'export PATH="$HOME/.cargo/bin:$PATH"' "$CONFIG_FILE"

    echo "Installing Rust binaries..."
    cat ./modules/rust-binaries.txt | xargs -n 1 cargo install
fi