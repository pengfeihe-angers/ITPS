# ITPS
The source code for ITPS

1) For the input parameters:
-i {the name of instance} --seed {seed} -rep {the number to count the times} -alb {the low bound} -pls {the maximal nonimproving iterations for LS} -pth {the maximal nonimproving iterations for threshold phase} -prb {the probability for not choosing the best candidate solution} -prc {the probability for choosing NF2 in threshold phase} -tnf {the percentage of random taken nodes in NF2} -dep {the depth of the two phases}
e.g.
-i cycle100.rnd --seed 0 -rep 0  -alb 1 -pls 0.2 -pth 3 -prb 0.5 -prc 0.1 -tnf 0.05 -dep 5

2) The command for compile the cpp file: "g++ BCP_irace.cpp -O3 -lm -Wall -o BCP"

3) The instances are in "instance.zip"

4) There is an example script written in python to generate the execution list for the readers: script_list.py

If you have questions, please contact the author: rjtkxj@gmail.com
