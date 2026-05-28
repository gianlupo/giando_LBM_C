#!/bin/bash

echo "Postprocessing..."

TEST_DIR=$(pwd)

python3 contours.py

cd ../../post

ffmpeg -framerate 12 -pattern_type glob -i 'vor_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 vor:gif
ffmpeg -framerate 12 -pattern_type glob -i 'vel_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 vel.gif
ffmpeg -framerate 12 -pattern_type glob -i 'rho_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 rho.gif
ffmpeg -framerate 12 -pattern_type glob -i 'u_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 u.gif
ffmpeg -framerate 12 -pattern_type glob -i 'v_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 v.gif
ffmpeg -framerate 12 -pattern_type glob -i 'w_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 w.gif
ffmpeg -framerate 12 -pattern_type glob -i 'K_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 K.gif

cd $TEST_DIR

python3 profiles.py

cd ../../post

ffmpeg -framerate 12 -pattern_type glob -i 'u_profile_*.png' -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 u_profile.gif

cd $TEST_DIR

echo "Done."
