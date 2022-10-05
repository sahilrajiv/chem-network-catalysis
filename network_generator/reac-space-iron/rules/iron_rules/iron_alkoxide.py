

alkoxide = """rule [
	ruleID "Fe3 ---> Fe2 Aldehyde to acid"
	left[
        node [id 10 label "O-"] #hydroxy ion
        

        node [id 3 label "H"]
        node [id 9 label "Fe3+"]
        
        node [id 12 label "Fe3+"] 

        edge[ source 1 target 3 label "-"]
        
        
	]
	context [
        node [id 1 label "C"]
        node [id 2 label "O"]
        
        node [id 4 label "C"]
        node [id 5 label "O"]
        node [id 6 label "H"]
        node [id 7 label "H"]
        node [id 8 label "*"]
        node [id 11 label "H"] 

        # 

        edge[ source 1 target 2 label "="]
        
        edge[ source 4 target 1 label "-"]
        edge[ source 4 target 5 label "-"]
        edge[ source 5 target 6 label "-"]
        edge[ source 4 target 7 label "-"]
        edge[ source 4 target 8 label "-"]

        edge[ source 10 target 11 label "-"]

        # edge[ source 2 target 13 label "-"]
        # edge[ source 5 target 13 label "-"]

        
        
          
	]


	right [
                node [id 10 label "O"] #hydroxy ion
                

        node [id 3 label "H+"]
        node [id 9 label "Fe2+"]
        node [id 12 label "Fe2+"]
        
        edge[ source 10 target 1 label "-"]

       
	]

      
	
]
"""

alkoxide_form = ruleGMLString(alkoxide)

# p = GraphPrinter()

# for r in inputRules:
#     r.print(p)
