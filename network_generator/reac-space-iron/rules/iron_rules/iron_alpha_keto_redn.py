

alpha_keto_redn = """rule [
	ruleID "Fe2 ---> Fe3 Alpha Keto Reduction of a Carboxylic Acid (Fe2+ to Fe3+)"
	left[
        node [id 7 label "H+"]
        node [id 8 label "Fe2+"]
        node [id 11 label "Fe2+"]

        edge[ source 4 target 5 label "="]
	]
	context [
        node [id 1 label "C"]
        node [id 2 label "O"]
        node [id 3 label "O"]

        node [id 10 label "H"]
        
        node [id 4 label "C"]
        node [id 5 label "O"]
        node [id 6 label "*"]

        edge[ source 3 target 10 label "-"]

        edge[ source 1 target 3 label "-"]
        edge[ source 1 target 2 label "="]

        edge[ source 1 target 4 label "-"]
        edge[ source 4 target 6 label "-"]
        
         
	]


	right [
        node [id 7 label "H"]
        node [id 8 label "Fe3+"]
        node [id 11 label "Fe3+"]

        edge[ source 4 target 5 label "-"]
        edge[ source 5 target 7 label "-"]
        	]
	
]
"""




alpha_keto_redn_fe= ruleGMLString(alpha_keto_redn)

