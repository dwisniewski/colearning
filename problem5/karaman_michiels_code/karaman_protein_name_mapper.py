'''
	Because protein names in karaman are given in form of full names
	we need to map michiel crystal structures to karaman to run Figuren.r program correctly

	example run command: python karaman_protein_name_mapper.py ../helper_dataset/tmp_michiel_proteins.txt ../dataset/kaset_structures_mapping ../helper_dataset/karaman_protein_crystals_only.txt^C
'''
import sys

michiel_list_path = sys.argv[1]
mapping_path = sys.argv[2]
output_path  = sys.argv[3]


def loadMichielList(path):
	f = open(path)
	data = []
	for line in f:
		data.append(line.strip())
	return data


def loadMapping(path):
	mapping = dict()
	
	f = open(path)
	for line in f:
		line = line.strip().split("\t")
		fullname = line[3]
		crystal_str = line[1]
		mapping[crystal_str] = fullname
	return mapping

def writeResult(path, data, mapping):
	f = open(path, 'w+')
	for item in data:
		if item in mapping.keys():
			f.write(mapping[item] + "\n")
		else:
			f.write("NONE\n")
	f.close()

data = loadMichielList(michiel_list_path)
mapping = loadMapping(mapping_path)

print(data)
print "\n\n"
print mapping


writeResult(output_path, data, mapping)