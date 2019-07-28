class solution():

    def __init__(self, objective_value, fluxes, reduced_costs, solver = None):

        self.solver = solver
        self.objective_value = objective_value
        self.fluxes = fluxes
        self.reduced_costs = reduced_costs

def get_solution(model, reactions = None, metabolites = None):

    pass