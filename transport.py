def find_transportMets(cobra_model, reaction, Kegg_map):

	""" 
	transport metabolite definition:
		different id (same base id but different compartment)
		
		comparing the KEGG ids of products and reactants for reaction. Assuming transport metabolites have same KEGG id irrespective of compartment
	"""
	
	products_KEGG = {}
	
	for i in reaction.products:
		products_KEGG[Kegg_map[i.id]] = i.id
	
	for i in reaction.reactants:
		if Kegg_map[i.id] in products_KEGG.keys():
			transport_mets.append(i)
			transport_mets.append(products_KEGG[i])
			
		
	return transport_mets

def transport_constraints(reaction):
	"""
	pH differnce notation in to out
	"""
	transport_mets = find_transportMets(cobra_model,reaction,Kegg_map)
	
	n_H = len([met.id for met in transport_mets if Kegg_map[met.id] == 'C00080'])/2
	
	c = sum(met.charge for met in transport_mets)
	
	compartment_list = list(set([i.compartment for i in transport_mets]))
	
	pH_diff = test # Finish this
	
	del_psi = 33.33 * pH_diff - 143.33 # adopted from henry paper
	del_psi_G = c * F * del_psi
	
	del_pH_G = -2.3*h*R*T*pH_diff # Define constants
	
	transport_delG = del_psi_G + del_pH_G
	
	return transport_delG
