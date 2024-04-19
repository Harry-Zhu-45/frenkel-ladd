#!/bin/bash
#SBATCH --job-name="FLFE-cpp"
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mem=1g
#SBATCH --nodes=1
#SBATCH --time=0-30:00:00
#SBATCH --mail-type=END,FAIL

num_runs=10

values=(1.3218454836828071 7.429242462300721 20.448998156279607 44.35579276774375 85.59512159890338 152.4321293389935 251.20242463993492 378.0986975526874 510.1482256352902 606.495129089744)

# boxtype="cubic"
# echo $boxtype
# replicate=4
# echo $replicate

boxtype="xyz"
echo $boxtype
nx=3
ny=3
nz=6
echo $nx, $ny, $nz

echo

for value in "${values[@]}"; do
    echo "Running for value $value..."
    for ((i = 1; i <= $num_runs; i++)); do
        # result=$(/home/mzhu/Desktop/FrenkelLadd/fl-cpp/build/main.out -l "$value" -b "$boxtype" -r "$replicate")
        result=$(/home/mzhu/Desktop/FrenkelLadd/fl-cpp/build/main.out -l "$value" -b "$boxtype" -x "$nx" -y "$ny" -z "$nz")

        msd=$(echo "$result" | grep "MSD:" | awk '{print $2}')
        echo "$msd"
    done
    echo
done
