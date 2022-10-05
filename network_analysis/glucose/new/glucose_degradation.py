# There's no formaldehyde at the beginning, turn the Cannizarro rule for formaldehyde off
with_formaldehyde = False

include("../main.py")
# including these files using MOD's include() so that MOD's functions are callable in them
include("../mod_to_neo4j_exporter.py")
#include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

open_glucose = smiles("O=CC(O)C(O)C(O)C(O)C(O)", "Open Chain Glucose")
water = smiles("O", name="Water")
oh=smiles("[OH-]", name="Hydroxy ion")
hplus=smiles("[H+]", name="Proton")
ferric=smiles("[Fe+3]", name="Ferric Ion")
ferrous=smiles("[Fe+2]", name="Ferrous Ion")
# Number of generations we want to perform
generations = 4


dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []

# dump initial reactants as part of "G0"
write_gen_output(subset, generation=0, reaction_name="glucose_degradation")

# In the following block, apart from generating the reactions, we may print structures
# and reactions forming them that are not in the MS
#postSection("Structures not found in MS")
with dg.build() as b:
	for gen in range(generations):
		start_time = time.time()
		print(f"Starting round {gen+1}")
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=8)
		end_time = time.time()
		print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
		print('Original subset size:', len(res.subset))

		# The returned subset and universe do not contain redundant tautomers
		#subset, universe = clean_taut(dg, res, algorithm="CMI")
		subset, universe = res.subset, res.universe
		#print('Subset size after removal:', len(subset))
		# This step replaces the previous subset (containing tautomers) with the cleaned subset
		#res = b.execute(addSubset(subset) >> addUniverse(universe))
		export_to_neo4j(dg_obj = dg, generation_num = gen+1)
		write_gen_output(subset, gen+1, reaction_name="glucose_degradation")
	mod_to_gephi(dg)
	print('Completed')

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)

