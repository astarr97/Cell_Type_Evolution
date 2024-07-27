o = open("commands_parameter_sweeps.sh")
c = 0
for line in o:
    c += 1
    out = open("run1-" + str(c) + ".sh", 'w')
    if "GeneSet" in line:
        time="168:00:00"
    elif "1\n" in line:
        time="168:00:00"
    else:
        time="168:00:00"
    out.write("#!/bin/bash\n#SBATCH --time=" + time + "\n#SBATCH -p hbfraser\n#SBATCH --mem=16GB\n#SBATCH -c 1\n\n")
    out.write(line)
    out.close()

o.close()
