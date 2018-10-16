#!/bin/bash
#
# Will run the program on 1, 2, 3, and 4 parallel 
# images in sequence and save the output. It will
# then diff the output to catch for any differences.

# run the program in sequence
for n in {1..4}; do
  echo Running app on $n images..
  mpiexec --oversubscribe -n $n ./app > out_${n}.txt
done

# compare the output
for n in {2..4}; do
  echo diff between num_images = 1 and 2:
  diff out_[12].txt
done

# remove temporary output files
rm out_*.txt
