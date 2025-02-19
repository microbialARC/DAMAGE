#!/usr/bin/env bash
BIN_DIR="$CONDA_PREFIX/bin"
cp -r DAMAGE/workflow/ $CONDA_PREFIX
cp DAMAGE/bin/cli.py "$BIN_DIR/DAMAGE"

chmod +x "$BIN_DIR/DAMAGE"

if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
    echo "export PATH=\"$BIN_DIR:\$PATH\"" >> "$HOME/.bashrc"
    echo "Added $BIN_DIR to PATH. Please restart your shell or run 'source ~/.bashrc'"
fi

echo "DAMAGE installed successfully! Run 'DAMAGE -h' to start with help page."