wget "https://www.digitalrocksportal.org/media/projects/47/archive.zip" --no-check-certificate
unzip -j archive.zip origin/311/images/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw -d input/.
#rm input/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw
rm archive.zip
