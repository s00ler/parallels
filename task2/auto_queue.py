from __future__ import print_function
import subprocess
import time
from functools import reduce

processes = [1, 2, 4, 8, 16, 32, 64, 128]
grid_sizes = [128, 256, 512]

tasks = []
for procs in processes:
    for grid_size in grid_sizes:
        tasks.append((procs, grid_size))
subprocess.check_output(['bash', '-c', "module add slurm"]).decode()
while tasks:
    bashCommand = "squeue -p test"
    output = subprocess.check_output(['bash', '-c', bashCommand]).decode()
    output = output.split('\n')
    output = list(map(lambda l: l.split(' '), output))
    output = [line for line in output if 'dagerasi' in line]
    free_spaces = 3 - len(output) + 1
    for i in range(free_spaces):
        task = tasks.pop()
        launchcmd = "sbatch -p test -n {} ./mpitask2 {} 0".format(task[0],
                                                                  task[1])
        output = subprocess.check_output(['bash', '-c', launchcmd]).decode()
        print(output)
        time.sleep(5)
    time.sleep(30)

print("All tasks launched!")
