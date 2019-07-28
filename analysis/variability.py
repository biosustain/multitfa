def thermo_variability_analysis(model, reaction_List = None): # processes = None

    if reaction_List == None:
        reaction_ids = [r.id for r in model.reactions]
    else:
        reaction_ids = [r.id 
                        for r in model.reactions.get_by_any(reaction_List)]
    with model:
        model.optimize()
    
    