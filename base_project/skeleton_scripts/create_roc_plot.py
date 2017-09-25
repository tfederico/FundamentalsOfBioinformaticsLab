#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, and calculates the AUC.

import argparse
import numpy
import matplotlib
matplotlib.use('AGG')
import pylab
from datetime import datetime
start_time = datetime.now()

def parse_blast_results(filename):
    """
    Parse every protein pair's e-value out of a BLAST results file.
    :param filename: input file with BLAST results.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding e-value from BLAST (value).
    """

    blast_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) != 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])
                
            key = (arr[0], arr[1])
            if arr[2] == "NA":    # Substitute protein pairs whose e-value is
                value = 1e6       # not available with an e-value of 1 million
            else:
                value = float(arr[2])
            blast_results[key] = value

    return blast_results


def parse_benchmark_results(filename):
    """
    Parse every protein pair's classification out of the benchmark file.
    :param filename: input file with benchmark classifacations.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding call (value).
    """

    benchmark_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
    
            if len(arr) < 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])

            # Benchmark classifications are symmetric, so add both possible keys:                
            key1 = (arr[0],arr[1])
            key2 = (arr[1],arr[0])
            value = arr[2]
            benchmark_results[key1] = value
            benchmark_results[key2] = value

    return benchmark_results


def integrate(x, y):
    """
    Calculate the Area Under the Curve (AUC) for a given list of coordinates
    :param x: a list of x-coordinates
    :param y: a list of y-coordinates
    :return: a float with the surface area under the curve described by x and y
    """
    
    auc = 0.
    last_x = x[0]
    last_y = y[0]
    for cur_x, cur_y in zip(x, y)[1:]:
        #########################
        ### START CODING HERE ###
        #########################
	x_diff = cur_x - last_x
	y_mean =  (cur_y + last_y) / 2
	auc += x_diff * y_mean
        #########################
        ###  END CODING HERE  ###
        #########################
        last_x = cur_x
        last_y = cur_y
    return auc


def roc_plot(blast_evalues, benchmark_dict, png_filename):
    """
    Draw the ROC plot for a given set of e-values and corresponding benchmark classifications.

    :param blast_evalues: the dictionary produced by parse_blast_results()
    :param benchmark_dict: the dictionary produced by parse_benchmark_results()
    """

    ### Create the lists of coordinates

    x = [0] # array of the ROC plot's x-coordinates: False Positive Rate = FP/(FP/TN)
    y = [0] # array of the ROC plot's y-coordinates: True  Positive Rate = TP/(TP+FN)
    
    evalues = [(v, k) for k, v in blast_evalues.items()] # List of tuples consisting of (evalue, protein_pair)
    sorted_evalues = sorted(evalues)
    last_evalue = sorted_evalues[0][0]
    counter = 1
    sims = 0
    diffs = 0
    fpr = 0
    tpr = 0
    my_var = False
    for evalue, protein_pair in sorted_evalues:

        #########################
        ### START CODING HERE ###
        #########################
        # Iterate through the protein pairs, in order of ascending e-value
        # Treat every unique e-value as a homology threshold:
        #   Append one coordinate to x and y for each e-value, tracking how
        #   many "different" and "similar" pairs you've come across so far,
        #   (i.e. do not add 2 coordinates if 2 protein pairs have the same e-value.)
        # Ignore entries in the benchmark_dict classified as "ambiguous",
        #   as well as protein pairs for which you have no benchmark classification

	if protein_pair in benchmark_dict.keys():
		if not benchmark_dict[protein_pair] == "ambiguous":		
			if not evalue == last_evalue:
				x.append(fpr)
				y.append(tpr)
			scop = benchmark_dict[protein_pair]
			if scop == "different":
				if my_var == False:
					my_var = True
					print("Res: "+str(evalue))
				diffs += 1
			elif scop == "similar":
				sims += 1
			else:
				print "ERROR"
      			fpr = diffs/counter
			tpr = sims/counter
			counter + 1
		else:
			continue
	#########################
        ###  END CODING HERE  ###
        #########################
        last_evalue = evalue

    # Remember ROC plots start in (0,0) and end in (1,1)!
    x = numpy.array(x) / float(x[-1]) # At this point, x[-1] = sum(benchmark_dict.values() == "different")
    y = numpy.array(y) / float(y[-1]) #                y[-1] = sum(benchmark_dict.values() == "similar")

    ### Figure out the AUC
    auc = integrate(x, y)
    
    ### Draw the plot and write it to a file
    pylab.plot(x, y)
    pylab.plot([0,1],[0,1],'--k')
    pylab.xlabel('False Positive Rate')
    pylab.ylabel('True Positive Rate')
    pylab.title('AUC = %.3f' % auc)
    pylab.savefig(png_filename)


def main(blast_results_file, benchmark_results_file, png_file):
    # Parse the input files and retrieve every protein pair's e-value and benchmark classification.
    blast_evalues = parse_blast_results(blast_results_file)
    benchmark_results = parse_benchmark_results(benchmark_results_file)
    
    # Draw and save the ROC plot
    roc_plot(blast_evalues, benchmark_results, png_file)
    print "The script " + __file__ + " ran successfully with an excecution time of " + str(datetime.now()-start_time) + "."

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw and save a ROC plot to a file")
    parser.add_argument("-iblast","--input_blast_results", help="tab-separated BLAST results file", required=True)
    parser.add_argument("-ibench","--input_benchmark_results", help="tab-separated benchmark classification file", required=True)
    parser.add_argument("-o", "--output_png", help="output png file", required=True)

    args = parser.parse_args()
    
    blast_file = args.input_blast_results
    benchmark_file = args.input_benchmark_results
    png_file = args.output_png

    main(blast_file,benchmark_file, png_file)
