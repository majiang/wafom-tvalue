from sys import stdin

x = 2
for line in stdin:
    if line.strip():
        print eval(line.strip().replace('^', '**'))
    else:
        print
