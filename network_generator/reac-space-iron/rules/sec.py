

secondary_alcohol_ox = ruleGMLString("""
rule [
    ruleID "Oxidation of Secondary Alcohols by Ferricyanide"
    #labelType "term"
    left [
        edge [ source 1 target 4 label "-" ]
        edge [ source 1 target 5 label "-" ]
        edge [ source 5 target 6 label "-" ]
        node [ id 7 label "Fe3+" ]  #Represent Iron like this
        node [ id 4 label "H" ]
        node [ id 6 label "H" ]
         ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "*" ]
        node [ id 3 label "*" ]
        node [ id 5 label "O" ]
        

        edge [ source 1 target 2 label "-" ]
        edge [ source 1 target 3 label "-" ]
       

   ]
    right [
        edge [ source 1 target 5 label "=" ]
        node [ id 4 label "H+" ]
        node [ id 6 label "H+" ]
        node [ id 7 label "Fe2+" ]
    ]
]
""")



Imine_Ox = ruleGMLString("""
rule [
    ruleID "Oxidation of an Imine by Ferricyanide"
    #labelType "term"
    left [
        edge [ source 1 target 3 label "-" ]
        edge [ source 1 target 4 label "=" ]
        edge [ source 4 target 5 label "-" ]

        node [ id 8 label "Fe3+" ]  #Represent Iron like this
        
        
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "*" ]
        node [ id 3 label "H" ]
        node [ id 4 label "N" ]
        node [ id 5 label "H" ]
        
        node [ id 10 label "C-" ]
        node [ id 11 label "N" ]
        node [ id 12 label "C-" ]
        node [ id 13 label "N" ]
        node [ id 14 label "C-" ]
        node [ id 15 label "N" ]
        node [ id 16 label "C-" ]
        node [ id 17 label "N" ]
        node [ id 18 label "C-" ]
        node [ id 19 label "N" ]
        node [ id 20 label "C-" ]
        node [ id 21 label "N" ]
         
        edge [ source 1 target 2 label "-" ]


        edge [ source 10 target 11 label "#" ]
        edge [ source 12 target 13 label "#" ]
        edge [ source 14 target 15 label "#" ]
        edge [ source 16 target 17 label "#" ]
        edge [ source 18 target 19 label "#" ]
        edge [ source 20 target 21 label "#" ]
        edge [ source 10 target 8 label "-" ]
        edge [ source 12 target 8 label "-" ]
        edge [ source 14 target 8 label "-" ]
        edge [ source 16 target 8 label "-" ]
        edge [ source 18 target 8 label "-" ]
        edge [ source 20 target 8 label "-" ]

    ]

    right [
        edge [ source 1 target 4 label "#" ]
        edge [ source 3 target 5 label "-" ]
        node [ id 8 label "Fe2+" ]
    ]
]
""")





k=GraphPrinter()
k.withIndex=True
for a in inputRules:
	a.print()




