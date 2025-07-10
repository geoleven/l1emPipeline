import sys

def process_file(input_file, output_file):
    newlines = []
    with open(input_file, "r") as infile:
        for line in infile.readlines()[1:]:
            feature = line.strip().split('\t')
            txstart = int(feature[2]) + 399
            txend = int(feature[3]) -400
            estart = int(feature[5]) + 399
            eend = int(feature[6]) -400
            newlines.append(
                f"{feature[0]}\t{feature[1]}\t{txstart}\t{txend}\t{feature[4]}\t{estart}\t{eend}\t{feature[7]}\t{feature[8]}\t{feature[9]}\t{feature[10]}\n"
            )

    with open(output_file, "w") as outf:
        for line in newlines:
            outf.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python annotation400.py <input_file> <output_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file)
