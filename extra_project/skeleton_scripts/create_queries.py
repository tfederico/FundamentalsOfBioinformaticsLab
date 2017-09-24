# 18. Sept. 2017
# This script write the queries from the database

import sys
from datetime import datetime
start_time = datetime.now()

def write_query(file_name, fasta_string, fasta_line):
	with open('../queries/'+file_name+'.fasta','w') as query_file:
		query_file.write(fasta_line)
		query_file.write(fasta_string)
def read_db(db_file_name):
	identifier = ''
	seq = ''
	with open(db_file_name) as db_file:
		for line in db_file:
			if line[0] == '>':
				if not identifier == '': 
					write_query(identifier, seq, fasta_line)
					seq = ''
					fasta_line = line
					line = line[1:].strip()
					line_array = line.split('|')
					identifier = line_array[3] 
				else:
					pass
			else:
				seq += line.strip()

def main(db_file_name):
	read_db(db_file_name)
	print 'The script '+__file__+' ran successfully with an excecution time of '+str(datetime.now()-start_time)

if __name__ == '__main__':
	db_file_name = sys.argv[1]
	main(db_file_name)
