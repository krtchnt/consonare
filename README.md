# consonare

Automatic Overtone Analysis and Dissonance Profiling for Single-Note Recordings

> [!IMPORTANT]
> This project is developed under the **01204496 Algorithmic-Oriented Digital Signal
  Processing for Computer Engineers** course of **Department of Computer Engineering**,
  **Faculity of Engineering**, **Kasetsart University**.

> **Project Developers**: *Kritchanat Thanapiphatsiri (6610501955)*

## Installation
1. Install [rust](https://www.rust-lang.org/tools/install).
2. Once finished, open the terminal, and run these commands:
```sh
git clone https://github.com/krtchnt/consonare
cd consonare
cargo build --release
```
Optionally, to enable dissonance graph plotting, use
`cargo build --release --features visualise` instead. This will make execution time
slower, with the benefit of also plotting and saving the computed dissonance graph to
disk.

## Usage
To use the program, type `target/release/consonare-cli` in the terminal first. Then:
- Type the path to the audio file then enter to run the automatic overtone analysis and
  dissonance profiling for that recording. The program will output the results in the
  terminal once it's done. Example sample files are available to use in
  [`samples/`](./samples).
- Before typing the audio path, type `--verbose` then run the program to show all
  procedure output from step 1 to 11/12 in the terminal, not only for step 10 and 11.
- To change the exported directory for CSVs and graphs, type `--out-dir` then the new path
  to the exporting directory. By default, any file exports are in the same directory as
  where you're running the command from.
- Type `--help` then enter to display all options available for this command to be used.
  Many of them can be used together at the same time, for example:
  `--verbose --out-dir ...`

Example usage:
```sh
# 1. display all options available and how to use them
target/release/consonare-cli --help

# 2. analyse an audio file
target/release/consonare-cli my-instrument.wav

# 3. like `2.` but show all output from step 1 to 11/12
target/release/consonare-cli my-instrument.wav --verbose

# 4. like `2.` but outputs CSVs and graphs in another directory instead of the current
target/release/consonare-cli my-instrument.wav --out-dir my-directory/my-sub-directory
```

Make sure the audio file is a mono WAV (`.wav`) file with negligible noise and consisting
of only single note in the recording.

## Attributions
- All sample recordings is owned by [MuseHub](https://www.musehub.com) from the
  (free MuseSounds)[https://www.musehub.com/free-musesounds].
