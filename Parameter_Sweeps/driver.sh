#Note: some runs (when controlling for tau or expression level for high clustering resolution for Sestan_DLPFC or Allen_MTG) did not finish in 7 days
#In that case, we reran starting at the first unfinished iteration and combined the two files to make a complete file

python make_scripts.py

for file in run1-*.sh;
do
    sbatch -p hbfraser $file
done
