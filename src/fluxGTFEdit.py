import re

def processGTF(inputFile, outputFile):
    with open(inputFile, 'r') as infile, open(outputFile, 'w') as outfile:
        tid = 0
        for line in infile:
            if line.startswith('#') or line.strip() == '':
                outfile.write(line)
                continue
            # Write original line
#            outfile.write(line)
            # Split line and process
            parts = line.strip().split('\t')
            tid += 1
            parts[0] = "transcript" + str(tid)
            if parts[6] == '+':
                parts[4] = str(int(parts[4]) - int(parts[3]) )
                parts[3] = '1'
            else:
                parts[4] = str(int(parts[4]) - int(parts[3]) )
                parts[3] = '0'
            outfile.write('\t'.join(parts) + '\n')

            if len(parts) >= 9:
                # Modify third column
                new_parts = parts.copy()
                new_parts[2] = 'exon'

                #new_parts[3] = '1'
                #new_parts[4] = str(int(parts[4]) - int(parts[3]))

                # Extract gene_id
                attributes = new_parts[-1]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                gene_id = gene_id_match.group(1) if gene_id_match else 'unknown'
                # Append exon_number and exon_id
                new_attributes = attributes.rstrip(';') + f'; exon_number "1"; exon_id "{gene_id}.1"'
                new_parts[-1] = new_attributes
                # Write modified line
                outfile.write('\t'.join(new_parts) + '\n')

processGTF('cusannot.sorted.gtf', 'cusannot.sorted.double.gtf')