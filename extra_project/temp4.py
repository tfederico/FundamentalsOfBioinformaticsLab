import datetime
start_time = datetime.datetime.now()

in_gis_file = open("gis.txt", "r")
total = 0
mappings = 0

for a in in_gis_file:
	total += 1
	gis_1 = a.strip()
	in_mapping_file = open("idmapping_output_output.txt", "r")
	for b in in_mapping_file:
		b_list = b.strip().split("\t")
		gis_2 = b_list[0]
		pdb_2 = b_list[1]
		if gis_1 == gis_2:
			out_file = open("mapped_gis.txt", "a")
			out_file.write(gis_2 + "\t" + pdb_2 + "\n")
			out_file.close()
			mappings += 1
		else:
			pass
	in_mapping_file.close()

in_gis_file.close()
print("Mappings: " + str(mappings))
print("Total Counter: " + str(mappings))
print __file__ + " ran with a time of " + str(datetime.datetime.now() - start_time)
