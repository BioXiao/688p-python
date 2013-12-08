

blast_results = open('Pelobates_cultripes-v-nt.blastn', 'r')

sequences = {}



for line in blast_results:
    line = line.strip()
    
    if line[0]=='#':
        continue
    
    columns = line.split('\t')
    
    sequence_name = columns[0]
    description = columns[-1]
    
    genus = description.split(' ')[0]
    
    if not sequence_name in sequences:
        sequences[sequence_name] = set()
        
    sequences[sequence_name].add(genus)
    
out_file = open("blast_result_genera.csv", 'w')
    
for sequence, value in sequences.iteritems():
    out_file.write(sequence+'\n')
    out_file.write(', '.join(value)+'\n')
    
out_file.close()