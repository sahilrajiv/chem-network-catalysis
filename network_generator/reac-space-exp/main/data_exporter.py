import os
import csv

def mod_to_gephi(dg_obj):

	
	with open("data.csv", 'a') as data_file:
		csv_writer = csv.writer(data_file, delimiter='\t')
		csv_writer.writerow(['ID', 'Species', 'R or P', 'Reaction Name', 'Mass'])
		for e in dg_obj.edges:

			for i in range(len(e.rules)):
				# use 'i' as a "sub-index"
				for source in e.sources:
					l=[f'{e.id}_{i}',source.graph.smiles,'-1',list(e.rules)[i].name,str(source.graph.exactMass)]
					csv_writer.writerow(l)
				for target in e.targets:
					l=[f'{e.id}_{i}',target.graph.smiles,'1',list(e.rules)[i].name,str(target.graph.exactMass)]
					csv_writer.writerow(l)

	with open("data.csv", 'r') as r_file:
		with open("rxn.csv", 'w') as in_file:
			csv_reader=csv.reader(r_file, delimiter='\t')

			csv_writer = csv.writer(in_file, delimiter='\t')
			csv_writer.writerow(['Source','Target', 'Reaction Name', 'Mass'])

			for row in csv_reader:
				l=[]
				if row[2]=='-1':
					l.append(row[1])
					l.append(row[0])
					l.append(row[3])
					l.append(row[4])
				if row[2]=='1':
					l.append(row[0])
					l.append(row[1])
					l.append(row[3])
					l.append(row[4])
				
				csv_writer.writerow(l)


	with open("data.csv", 'r') as r_file:
		with open("type.csv", 'w') as in_file:
			csv_reader=csv.reader(r_file, delimiter='\t')
			csv_writer = csv.writer(in_file, delimiter='\t')
			next(csv_reader)
			csv_writer.writerow(['ID','Type'])
			ruleID=set()
			node=set()
			print(type(csv_reader))
			for row in csv_reader:
		
				ruleID.add(row[0])
				node.add(row[1])

			for m in node:
				csv_writer.writerow([m, 1])

			for m in ruleID:
				csv_writer.writerow([m, 2])



		

	


