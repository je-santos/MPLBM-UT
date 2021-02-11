start=`date +%s`

run_pyvista_animation() {
    j=$1
    
    #echo 'Process # ' $j
    #sleep 3
    python pyvista_lbm_animation.py $j
}


num_proc=4 #$(ulimit -u) this command get max number of processes, but limiting at 4 because memory 

for i in {0..18}
do
    j=$i #keep track of the for loop with j so i can be used for limiting processes
    
    ((i=i%num_proc)); ((i++==0)) && wait
    
    run_pyvista_animation "$j" &
    
    #echo 'Process # ' $j
done
wait

python create_gif.py

end=`date +%s`
runtime=$((end-start))

echo 'Creating the animation took this long: ' $runtime ' seconds'
