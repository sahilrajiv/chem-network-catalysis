amino_rule = """rule [
	ruleID "Fe2 ---> Fe3 Alpha Keto Reduction, conversion to amino acid"
	left[
        node [id 15 label "H"]
        node [id 16 label "H+"]
        node [id 17 label "H+"]

        
        node [id 10 label "Fe2+"]
        node [id 8 label "Fe2+"]

        edge[ source 12 target 15 label "-"]
        edge[ source 6 target 5 label "="]


        ]
	context [
        node [id 1 label "C"]
        node [id 2 label "O"]
        node [id 3 label "H"]
        node [id 4 label "O"]
        node [id 5 label "C"]
        node [id 6 label "O"]
        node [id 7 label "*"]

        
   

        node [id 12 label "N"]
        node [id 13 label "H"]
        node [id 14 label "H"]


        edge[ source 1 target 2 label "-"]

        edge[ source 1 target 4 label "="]
        edge[ source 2 target 3 label "-"]
        
        edge[ source 1 target 5 label "-"]

        edge[ source 5 target 7 label "-"]
        edge[ source 12 target 13 label "-"]
        edge[ source 12 target 14 label "-"]


        
         
	]


	right [
        node [id 15 label "H"]
        node [id 16 label "H"]
        node [id 17 label "H"]

   
        node [id 10 label "Fe3+"]
        node [id 8 label "Fe3+"]

      

        edge[ source 5 target 12 label "-"]
        edge[ source 5 target 15 label "-"]

        edge[ source 6 target 16 label "-"]
        edge[ source 6 target 17 label "-"]

      
        	]
	
]
"""


amino_rule_pyrite= ruleGMLString(amino_rule)

# for p in inputRules:
#         p.print()


