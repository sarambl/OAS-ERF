from subprocess import run
file1 = open('archive_cases.txt', 'r')
Lines = file1.readlines()

count = 0
# Strips the newline character
for line in Lines:
    count += 1
    print(f'find {line[:-1]} -type d ')
    run(f'find {line[:-1]} -type d ', shell=True)
    #print("Line{}: {}".format(count, line.strip()))