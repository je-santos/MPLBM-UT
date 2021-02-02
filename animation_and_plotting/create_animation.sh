start=`date +%s`

for i in {0..1}
do
    python pyvista_lbm_animation.py $i
done

python create_gif.py

end=`date +%s`
runtime=$((end-start))

echo 'Creating the animation took this long: ' $runtime ' seconds'
