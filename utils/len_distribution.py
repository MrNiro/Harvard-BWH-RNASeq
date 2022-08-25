import sys

if __name__ == '__main__':
    filename = ""
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    f = open(filename)
    len_distribution = {}
    for i, line in enumerate(f.readlines()):
        if i % 4 == 1:
            reads_len = len(line.strip())
            if reads_len not in len_distribution:
                len_distribution[reads_len] = 1
            else:
                len_distribution[reads_len] += 1

    output = list(map(lambda x: (x, len_distribution[x]), len_distribution.keys()))
    output.sort(key=lambda x: x[0])
    print(output)
