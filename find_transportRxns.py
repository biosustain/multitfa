def find_transportRxns(model):
	
	""" Transport reaction definition
	if a metabolite other than ATP, H+ is present in two compartments at the same time (To account for various transport mechanisms) 
	
	or
	
	If its just one metabolite present in the reaction *To account for the passive diffusion of metabolite
	
	"""
	transport_rxns=[];
	
	for rxn in model.reactions:
		#rxn=model.reactions.get_by_id(rxn)
		
		met=list(set(i.name for i in rxn.metabolites))
		#print(rxn,met)
		
		transport_dict={'H+': 'C00080', 'ATP C10H12N5O13P3': 'C00002', 'ADP C10H12N5O10P2': 'C00008'};
		
		temp_count=0;
		for x in met:
			if ((transport_dict.get(x) != None) or (x in transport_dict.values())):
				continue
			else:
				temp_count=temp_count+1;
				#print(rxn, x)
				 
		if not temp_count > 1:
			transport_rxns.append(rxn)
	
	return transport_rxns
