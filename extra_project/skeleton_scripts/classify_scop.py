

#!/usr/bin/python
from datetime import datetime
import itertools
import argparse

start_time = datetime.now()

def retrieve_scop_data(scop_file): # done
    """
    Reads a databaset file from SCOP and returns a dictionary with the protein IDS mapping.
    :param scop_file: database file containing mapping of PDB's to SCOP ID's.
    :return: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    """

    scop_dic = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You can parse SCOP data in various ways. E.g. you can use dictionary of dictionaries
    # {proteinID: {"class": class, "fold": fold, "superfamily": superfamily, 'family': family}}
    with open(scop_file, "r") as scop_database:
        for line in scop_database:
            if line[0] != "#":
                line = line.strip()
                line_array = line.split("\t")
                pdb_id = line_array[1].upper()
                scop_data = line_array[5]
                scop_data_array = scop_data.split(",")[0:4]
                scop_class = scop_data_array[0].split("=")[1]
                scop_fold = scop_data_array[1].split("=")[1]
                scop_superfamily = scop_data_array[2].split("=")[1]
                scop_family = scop_data_array[3].split("=")[1]
                scop_dic[pdb_id] = [scop_class, scop_fold, scop_superfamily, scop_family]
            else:
                pass
    ########################
    ### END CODING HERE ####
    ########################

    return scop_dic


def compute_similarity_score(prot1_scop, prot2_scop): # done
    """
    Computes the score for two proteins on the basis of the data from SCOP database.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # You need to decide whether you need this function for SCOP database.
    score = 0
    index = 0
    go_on = True
    while go_on and index < 4:
        if prot1_scop[index] == prot2_scop[index]:
            score += 1
        else:
            go_on = False
        index += 1

    return score
    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(prot1_scop, prot2_scop): # done
    """
    Returns the similarity score between two proteins.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :param pair: a tuple with the UniProt IDs of the two proteins to compare. (Not present)
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    score = compute_similarity_score(prot1_scop, prot2_scop)

    if score >= 4:
        parameter = "similar"
    elif score <= 1:
        parameter = "different"
    else:
        parameter = "ambiguous"

    return parameter
    ########################
    ### END CODING HERE ####
    ########################

# If you will use the numeric score for SCOP (similar to GO), you may want to use check_similarity_for_protein_pair
# with other arguments. See the example below.
# def check_similarity_for_protein_pair(score, threshold):
#    pass

def generate_all_possible_protein_pairs(protein_ids): # done
    """
    Returns a list containing all unique protein pairs.
    :param protein_ids: list of all proteins IDs.
    :return: list of possible unique protein pairs.
    """
    pairs = list()
    uniprot_ids = protein_ids
    ##########################
    ### START CODING HERE ####
    ##########################
    for a in range(0, len(uniprot_ids)):
        for b in range(a+1, len(uniprot_ids)):
            pairs.append((uniprot_ids[a], uniprot_ids[b]))
    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def assign_homology(scop_dict, conv_dict, uniprot_pairs): # done
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param scop_dict: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    :param protein_ids_pdbs: dictionary with UniprotID as key and PDB ID as a value.
    :param pairs: list of all possible unique protein pairs.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar).
    """
    scop_homology = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.
    scop_pdbs = scop_dict.keys()
    for pair in uniprot_pairs:
        pdb1 = conv_dict[pair[0]]
        pdb2 = conv_dict[pair[1]]
        if pdb1 in scop_pdbs and pdb2 in scop_pdbs:
            similarity_score = check_similarity_for_protein_pair(scop_dict[pdb1], scop_dict[pdb2]) # returns one out of {"similar", "different", "ambiguous"}
            scop_homology[pair] = similarity_score
        else:
            pass
    ########################
    ### END CODING HERE ####
    ########################
    return scop_homology


def write_results(filename, scop_homology): # done
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param output_file: the name of the output file.
    :param scop_homology: dictionary (keys: protein pairs as tuples; values: one of the value - different/similar/ambiguous)
    """
    with open(filename, "w") as f:
        for (p1, p2), value in scop_homology.iteritems():
            f.write("\t".join([p1, p2, value]) + "\n")


def read_protein_ids_file(filename): # done
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    protein_ids = list()
    with open(filename, "r") as uni_file:
        for line in uni_file:
            line = line.strip()
            protein_ids.append(line)
    #######################
    ### END CODING HERE ###
    #######################
    return protein_ids


def read_lookup_table(filename): # done
    """
    Reads the specified file and returns the dictionary with UniprotID as key and PDB ID as a value.
    :param filename: file with the mapping between Uniprot ids and PDB ids.
    :return: dictionary with UniprotID as key and PDB ID as a value.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    conversion_dictionary = {}
    with open(filename, "r") as lookup_table:
	    for line in lookup_table:
       		line = line.strip()
      		line_array = line.split("\t")
      	 	uniprot = line_array[0]
		pdb = line_array[1]
		conversion_dictionary[uniprot] = pdb
    #######################
    ### END CODING HERE ###
    #######################
    return conversion_dictionary


def main(uniprot_ids_file_name, output_file_name, lookup_file_name, scop_file_name):
    ##########################
    ### START CODING HERE ####
    ##########################
    conversion_dictionary = read_lookup_table(lookup_file_name) # dic, uniprot -> pdb
    #print(conversion_dictionary) # debugging
    uniprot_ids = read_protein_ids_file(uniprot_ids_file_name) # list, uniprot
    #print uniprot_ids +  "\n" +  "Length: " + len(uniprot_ids)) # debugging
    uniprot_pairs = generate_all_possible_protein_pairs(uniprot_ids) # list of 2-tuples, excl. (a, a), (b, b), if (a, b) then excl. (b, a)
    #print(uniprot_pairs)
    scop_dictionary = retrieve_scop_data(scop_file_name) # pdb -> scop_data (class, fold, superfamily, family)
    #print scop_dictionary
    scop_homology_dictionary = assign_homology(scop_dictionary, conversion_dictionary, uniprot_pairs)
    #print scop_homology_dictionary
    write_results(output_file_name, scop_homology_dictionary)
   
    print scop_homology_dictionary.values()    
    sim = 0
    amb = 0
    dif = 0
    for item in scop_homology_dictionary.values():
	if item == "similar":
		sim += 1
	elif item == "ambiguous":
		amb += 1
	elif item == "different":
		dif += 1
    print ("Similar: " +str(sim))
    print("different: " +str(dif))
    print("amb : " +str(amb))

    print "Script " + str(__file__) + " ran successful, with an execution time of "+ str(datetime.now()-start_time) + "."
    #######################
    ### END CODING HERE ###
    #######################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from SCOP database (from local file)'
                                                 ' and provides an output file'
                                                 ' with the strings "ProteinID   ProteinID   similarity", where'
                                                 ' similarity is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-o", "--output_file", help="Output file name")
    parser.add_argument("-ids", "--protein_ids_file", help="File with the protein Uniprot IDs")
    parser.add_argument("-pdb", "--pdb_id_file", help="File with the mapping between Uniprot ids and PDB ids")
    parser.add_argument("-s", "--scop_file", help="SCOP database file")

    args = parser.parse_args()

    input_file = args.protein_ids_file
    output_file = args.output_file
    pdb_id_file = args.pdb_id_file
    scop_file = args.scop_file

    main(input_file, output_file, pdb_id_file, scop_file)
