#!/bin/bash

echo "Building Atom Simulation..."

# Try different SDL2 paths for various systems
if pkg-config --exists sdl2 gl glu; then
    echo "Using pkg-config for dependencies..."
    gcc -o atom_simulation main.c $(pkg-config --cflags --libs sdl2 gl glu) -lm -Wall -O2
elif [ -f /usr/include/SDL2/SDL.h ]; then
    echo "Using system SDL2..."
    gcc -o atom_simulation main.c -I/usr/include/SDL2 -lSDL2 -lGL -lGLU -lm -Wall -O2
else
    echo "Trying basic compilation..."
    gcc -o atom_simulation main.c -lSDL2 -lGL -lGLU -lm -Wall -O2 2>/dev/null
fi

if [ $? -eq 0 ]; then
    echo "Build successful! Run with: ./atom_simulation"
else
    echo "Build failed. Make sure SDL2 development libraries are installed:"
    echo "Ubuntu/Debian: sudo apt install libsdl2-dev libglu1-mesa-dev"
    echo "Fedora: sudo dnf install SDL2-devel mesa-libGLU-devel"
    echo "Arch: sudo pacman -S sdl2 glu"
fi
