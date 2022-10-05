include('iron_alpha_keto_redn.py')
include('iron_amino_formation.py')
include('iron_alkoxide.py')
include('iron_beta_elimination.py')
include('iron_aldalphaoxi.py')
# glyceraldehyde=smiles("O=C(C(=O)O)C", name="Glyceraldehyde")
# #glycolaldehyde=smiles("O=CCO", name="Glycolaldehyde")
# oh=smiles("[OH-]", name="Hydroxy ion")
# hplus=smiles("[H+]", name="Proton")
# ferric=smiles("[Fe+3]", name="Ferric Ion")
# ferrous=smiles("[Fe+2]", name="Ferrous Ion")
alphahydro_alphaketo_acid= [ruleGMLString("""rule [
	ruleID "Fe3 ---> Fe2 Conversion of Alpha Hydro Group of an ACID to Alpha Keto group" 
	labelType "term"
	left [
		node [ id 6 label "H" ]
        node [ id 8 label "Fe3+" ]
        node [ id 9 label "H" ]
		node [id 10 label "Fe3+"]
        edge [ source 5 target 6 label "-" ]
        edge [ source 4 target 5 label "-" ]
        edge [ source 4 target 9 label "-" ]
	]   
	context [
		
		edge [ source 1 target 2 label "=" ]
        edge [ source 1 target 3 label "-" ]
        edge [ source 1 target 4 label "-" ]
		
        edge [ source 4 target 7 label "-" ]
        edge [ source 3 target 11 label "-" ]

		node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "O" ]
        node [ id 4 label "C" ]
        node [ id 5 label "O" ]
        node [ id 7 label "*" ]
        node [ id 11 label "H" ]
		
	]	
	right [
		node [ id 6 label "H+" ]
        node [ id 8 label "Fe2+" ]
        node [ id 9 label "H+" ]
		node [id 10 label "Fe2+"]

        edge [ source 4 target 5 label "=" ]
	]   

	 


	  constrainAdj [ id 7 op "=" count 0
        nodeLabels [ label "N" label "S" label "O" ]
  
		edgeLabels [ label "=" ]
    ]



]""")]


alphaketo_alphahydro_acid= [ruleGMLString("""rule [
	ruleID "Fe2 ---> Fe3 Conversion of Alpha Keto Group of an ACID to Alpha Hydro group" 
	labelType "term"
	left [
        node [ id 6 label "H+" ]
        node [ id 8 label "Fe2+" ]
        node [ id 9 label "H+" ]
		node [id 10 label "Fe2+"]

        edge [ source 4 target 5 label "=" ]
		
	]   
	context [
		
		edge [ source 1 target 2 label "=" ]
        edge [ source 1 target 3 label "-" ]
        edge [ source 1 target 4 label "-" ]
	edge [ source 3 target 11 label "-" ]	
        edge [ source 4 target 7 label "-" ]

		node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "O" ]
        node [ id 4 label "C" ]
        node [ id 5 label "O" ]
        node [ id 7 label "*" ]
        node [ id 11 label "H" ]
		
	]	
	right [
		

        node [ id 6 label "H" ]
        node [ id 8 label "Fe3+" ]
        node [ id 9 label "H" ]
		node [id 10 label "Fe3+"]
        edge [ source 5 target 6 label "-" ]
        edge [ source 4 target 5 label "-" ]
        edge [ source 4 target 9 label "-" ]

	]   

	 


	  constrainAdj [ id 7 op "=" count 0
        nodeLabels [ label "N" label "S" label "O" ]
  
		edgeLabels [ label "=" ]
    ]

]""")]

k=GraphPrinter()
k.withIndex=True
for a in inputRules:
	a.print()



# for a in inputGraphs:
# 	a.print()

# subset = inputGraphs
# universe = []

# dg = DG(graphDatabase=inputGraphs, labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

# generation=0


# with dg.build() as b:
#      b.execute(
#          addSubset(inputGraphs)
#          >> repeat[2](
#              inputRules
#          ), verbosity=6
#     )
		



# n=DGPrinter()

# dg.print(n)

