# 1. Backup, Saved on 5. Sept. [Federico]
# 2. Backup, Saved on 6. Sept. [Ben] - Assignment for Week 1 Done

#!/usr/bin/python

# This script downloads all the sequences in the list provided by the user, and puts them into a single file to be used
# as a database.
# It also stores the fetched fasta sequences into individual files which will be used as queries for (PSI-)BLAST search.

import urllib2
import argparse
import os.path


def fetch_one_fasta(uniprot_id):
    """
    Fetch the fasta formatted sequence of the uniProtID supplied in the argument.
    :param uniProtID: id of the protein to fetch in uniprot database.
    :return: sequence data in fasta format.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # Define the variable 'url' as a string with the URL to the Fasta formatted sequence of the uniProtID
    # in the uniprot website (www.uniprot.org).

    url = "http://www.uniprot.org/uniprot/"+uniprot_id+".fasta"
#   print url
    ##########################
    ###  END CODING HERE  ####
    ##########################


    fh = urllib2.urlopen(url)
    result = fh.read()
    fh.close()
    return result

def fetch_all_sequences(query_folder, uniprot_ids, database_file):
    """
    Fetch the fasta formatted sequence for each uniProtID in the list.
    :param query_folder: folder containing fasta files.
    :param uniprot_ids: File with the list of protein IDs
    :param database_file: database files to save the results of the fetching.
    """
    if not query_folder.endswith("/"):
        query_folder += "/"

    print("Processing the list of ids...")
    # uniprot_ids is file_object, gets already opened in skeleton, in "r" (read) mode
    # database_file is also file_object, gets already opened in skeleton, in "a" (append) mode
    for line in uniprot_ids:

        ##########################
        ### START CODING HERE ####
        ##########################
        # Fetch the fasta formatted sequence for each uniProtID.
        # Store the fasta sequences as individual fasta file in your query directory.
        # Store all the fasta sequences in one single fasta file as well. These individual files will be used
        # as (PSI-)BLAST queries later on.
	    line = line.rstrip() # see documentation, removes "\n"
	    my_fasta = fetch_one_fasta(line) # retrieves the fasta data by calling the above function
	    my_query = open(query_folder+"/"+line+".fasta", "w") 	# opens new fasta file object with file name
								                                # PDB-code.fasta in write mode in query_folder
	    my_query.write(my_fasta)
	    my_query.close()
	    database_file.write(my_fasta)
        ##########################
        ###  END CODING HERE  ####
        ##########################

    print("Processing finished.")
    uniprot_ids.close()
    database_file.close()

def check_query_folder(query_folder):
    if not os.path.exists(query_folder):
        return False
    return True


def main(query_folder, uniprot_ids, database_file):
    if check_query_folder(query_folder):
        fetch_all_sequences(query_folder, uniprot_ids, database_file)
    else:
        print("Query folder does not exist. "
              "Please make sure that the specified folder exists before you run this script.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script downloads all the sequences in the list provided'
                                                 'by the user, and puts them into a single file to be used'
                                                 'as a database. It also stores the fetched fasta sequences into'
                                                 'individual files which will be used as queries for (PSI-)BLAST search')

    parser.add_argument("-ids", "--uniprot_ids", type=argparse.FileType("r"),
                        help="File with the list of protein IDs", required=True)
    parser.add_argument("-db", "--dbfile", type=argparse.FileType("w"),
                        help="Output file containing all fetched sequences", required=True)
    parser.add_argument("-q", "--qfolder", help="Output folder for your queries", required=True)
    args = parser.parse_args()

    query_folder = args.qfolder
    uniprot_ids = args.uniprot_ids
    database_file = args.dbfile

    main(query_folder, uniprot_ids, database_file)
