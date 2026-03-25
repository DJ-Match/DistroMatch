#!/bin/bash

# Check if the user provided the required argument
if [ "$#" -ne 1 ]; then
	echo "Usage: $0 <number_of_nodes>"
	exit 1
fi

# Get the number of nodes/tasks from the argument
NODES=$1

# Define the output file name
OUTPUT_FILE="job_p${NODES}.sh"

# Generate the SLURM batch script
cat <<EOF > $OUTPUT_FILE
#!/bin/bash

#SBATCH -J scaling
#SBATCH -A p72884
#SBATCH -t 02:00:00
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512
#SBATCH --nodes=$NODES               # Request nodes
#SBATCH --ntasks=$NODES              # MPI tasks total
#SBATCH --ntasks-per-node=1          # 1 task per node
#SBATCH --cpus-per-task=1            # 1 core per task
#SBATCH --exclusive 

# Execute your program
./run.sh run_p$NODES 1
EOF

echo "SLURM batch script '$OUTPUT_FILE' generated successfully."
