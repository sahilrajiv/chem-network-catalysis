alphahydroxytoalphaketo = [ruleGMLString("""rule [
	ruleID "Fe3 ---> Fe2 Conversion of Alpha Hydroxy Group of an aldehyde to Alpha Keto group" 
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

		node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        node [ id 5 label "O" ]
        node [ id 7 label "*" ]
        
		
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


alphahydroxytoalphaketo_inverse = [ruleGMLString("""rule [
	ruleID "Fe2 ---> Fe3 Conversion of Alpha Hydroxy Group of an aldehyde to Alpha Keto group, inverse" 
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
		
        edge [ source 4 target 7 label "-" ]

		node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        node [ id 5 label "O" ]
        node [ id 7 label "*" ]
        
		
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
