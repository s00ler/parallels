from __future__ import print_function
import subprocess
import time
import sys


processes = [1, 2, 4, 8, 16, 32, 64, 128]
grid_sizes = [128, 256, 512]

tasks = []
for procs in processes:
    for grid_size in grid_sizes:
        tasks.append((procs, grid_size))
try:
    start = int(sys.argv[2])
except Exception:
    start = 0

while tasks:
    bashCommand = "squeue -p test"
    output = subprocess.Popen(bashCommand.split(' '),
                              stdout=subprocess.PIPE).communicate()[0]
    output = list(map(lambda l: l.split(' '), output.split('\n')))
    output = [line for line in output if 'dagerasi' in line]
    free_spaces = 3 - len(output)
    for i in range(free_spaces):
        task = tasks.pop()
        launchcmd = "sbatch -p test -n {0} impi ./{} {1} 0".format(task[0], sys.argv[1], task[1])
        output = subprocess.Popen(launchcmd.split(' '),
                                  stdout=subprocess.PIPE).communicate()[0]
        print(output)
        time.sleep(2)
    time.sleep(60)

print("All tasks launched!")
