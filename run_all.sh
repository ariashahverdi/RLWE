#!/bin/bash

#SBATCH --job-name="math_crypt"
#SBATCH -n 1 
#SBATCH -t 2-00:00:00
#SBATCH --ntasks=5
#SBATCH --mem=128gb
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

. /usr/share/Modules/init/bash
module load sage

# Location of Challenges File for example
#'/home/<username>/Desktop/rlwe-challenges-v1/rlwe-challenges-v1/challenges/chall-id0001-rlwed-m256-q769-l3-short-toy'
chall_dir=$1
echo $chall_dir
folder_name="${chall_dir##*/}"

# LeakPattern Supported = 1&7 mod 16, 1&15 mod 16, 1 mod 8, 1 mod 16 
leakpattern=$2
echo $leakpattern

# Either running "Challenges" or "NewHope" parameter set
mode=$3
echo $mode

# Either 'Binomial' or 'Gaussian'
# Doesn't matter if you are running challenges
sampling=$4
echo $sampling

filename='$chall_dir-$mode.txt'

if [ "$mode" == "Challenges" ]; then

	for i in $chall_dir/*.secret; do
	    filename=${i%.*}
	    
	    echo "************************************************"

	    name="${filename##*/}"
	    echo $name
	    
	    ./run.sh "$chall_dir/$folder_name.challenge" "$filename.instance" "$filename.secret" $leakpattern $mode $sampling '1' '1'

	    echo ""
	    echo ""
	    echo ""
	done
else
	for i in {1..200}; do
		
		for j in {1..1}; do
			
			echo "************************************************"
			echo $i $j

			#srun -N 1 --mem=32gb run.sh "1" "2" "3" $leakpattern $mode $sampling $i $j  
			./run.sh "1" "2" "3" $leakpattern $mode $sampling $i $j 
		
		done
	
		wait

	done
fi
