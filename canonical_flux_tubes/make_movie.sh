# usage: make_movie.sh date prefix output_name frames_per_second 
ffmpeg -framerate $4 -i ../../output/canonical_flux_tubes/$1/$2%04d.png -r $4 -pix_fmt yuv420p \
 -vcodec libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -y ../../output/canonical_flux_tubes/$1/$3.mp4
