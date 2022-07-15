import sys
import os


if __name__ == '__main__':
    data_folder = "/home/ubuntu/data/"

    if len(sys.argv) > 2:
        data_folder = sys.argv[1]

    input_sample = open("/home/ubuntu/rnaseq/samplesheet.csv", "w")
    input_sample.write("sample,fastq_1,fastq_2,strandedness\n")
    for _, _, filenames in os.walk(data_folder):
        filenames.sort(key=lambda x: int(x.split("_")[0]))
        for idx, f_name in enumerate(filenames):
            sample = "control_REP%d," % (idx + 1)
            fastq = data_folder + f_name + ",,"
            strandedness = "forward\n"

            line = sample + fastq + strandedness
            input_sample.write(line)

    input_sample.close()
