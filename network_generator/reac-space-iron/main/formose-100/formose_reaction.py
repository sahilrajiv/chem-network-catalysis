
# This flag will be used to toggle one pair of extra Cannizarro 2 rules
with_formaldehyde = True

include("../main.py")
max_mass_limit = 100
include("../mod_to_neo4j_exporter.py")

postChapter("Formose Reaction")

formaldehyde = smiles("C=O", name="Formaldehyde")
glycoladehyde = smiles("OCC=O", name="Glycolaldehyde")
water = smiles("O", name="Water")
oh=smiles("[OH-]", name="Hydroxy ion")
hplus=smiles("[H+]", name="Proton")
ferric=smiles("[Fe+3]", name="Ferric Ion")
ferrous=smiles("[Fe+2]", name="Ferrous Ion")

'''dg = DG.load(inputGraphs, inputRules, "formose_6rounds_dec21.dg")
print("Finished loading from dump file")'''

# Number of generations we want to perform
generations = 50

dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []

# dump initial reactants as part of "G0"
write_gen_output(subset, generation=0, reaction_name="formose")

with dg.build() as b:
	for gen in range(generations):
		start_time = time.time()
		print(f"Starting round {gen+1}")
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=4)
		end_time = time.time()
		print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
		print(f'Products in generation {gen+1}:', len(res.subset))

		# The returned subset and universe do not contain redundant tautomers
		#subset, universe = clean_taut(dg, res, algorithm="CMI")
		subset, universe = res.subset, res.universe
		#print('Product set size after removal:', len(subset))
		# This step replaces the previous subset (containing tautomers) with the cleaned subset
		#res = b.execute(addSubset(subset) >> addUniverse(universe))
		
		export_to_neo4j(dg_obj = dg, generation_num = gen)
		write_gen_output(subset, gen+1, reaction_name="formose")

	mod_to_gephi(dg) #for ds libraries
	print('Completed')

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)
# dg.print()
# count_rules(dg)
