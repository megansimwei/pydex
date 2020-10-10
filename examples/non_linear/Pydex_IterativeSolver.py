import numpy as np
import math as math

# Maximum number of optimal candidates - Caratheodory's Theorem
def MaxCandidates(designer):
    num_parameters = len(designer.model_parameters)
    maxcandidates = (num_parameters * (num_parameters + 1)) / 2
    return maxcandidates


# Determine initial density of the grid
def GetInitialLevel(MaxCandidates, designer, N_Ti_variables):
    maxcandidates = MaxCandidates(designer)
    i = math.ceil(math.sqrt(maxcandidates))
    levels = [i] * N_Ti_variables
    return levels


# initialization of deltas
def GetInitialDelta(bounds, levels, N_Ti_variables):
    delta = np.zeros(N_Ti_variables)
    for x in range(N_Ti_variables):
        delta[x] = (bounds[x][1] - bounds[x][0]) / ((levels[x] - 1) * 2)
    return delta


# Get new partition and level around an optimal solution
def GetNewPartition(bounds, optimal_candidate, N_Ti_variables, delta, sub_levels):
    boundsAroundCandidates = []
    partitionlevels = []
    for x in range(N_Ti_variables):
        if optimal_candidate[x] == bounds[x][0]:
            newbounds = [optimal_candidate[x], optimal_candidate[x] + delta[x]]
            partitionlevels.append(sub_levels)
        elif optimal_candidate[x] == bounds[x][1]:
            newbounds = [optimal_candidate[x] - delta[x], optimal_candidate[x]]
            partitionlevels.append(sub_levels)
        else:
            newbounds = [optimal_candidate[x] - delta[x], optimal_candidate[x] + delta[x]]
            partitionlevels.append(sub_levels * 2 - 1)
        boundsAroundCandidates.append(newbounds)
    return boundsAroundCandidates, partitionlevels


# Set the termination range of variables with undefined termination range to be 0.01 of the total range
def GetAllTerminationRanges(termination_range, bounds):
    x = 0
    for range in termination_range:
        if range == 0:
            termination_range[x] = (bounds[x][1] - bounds[x][0]) * 0.01

        x += 1
    return termination_range


# Update Delta to be a quarter of the size, if delta is less than the termination range_given it gets the value of the termination range
def UpdateDelta(delta, termination_range, N_Ti_variables, sub_levels):
    division = 2 * (sub_levels - 1)
    for x in range(N_Ti_variables):
        if (delta[x] / division) > termination_range[x]:
            delta[x] = delta[x] / division
        else:
            delta[x] = termination_range[x]
    return delta


# Find optimal experimental design candidates iteratively from
def IterativeSolver(designer, bounds, package, optimizer, termination_range, showallplots = False, sub_levels=3):
    N_Ti_variables = len(bounds)
    termination_range = GetAllTerminationRanges(termination_range, bounds)
    levels = GetInitialLevel(MaxCandidates, designer, N_Ti_variables)
    designer.ti_controls_candidates = designer.enumerate_candidates(bounds, levels)
    designer.initialize(verbose=2)
    designer.design_experiment(criterion=designer.d_opt_criterion, package=package, optimizer=optimizer)
    designer.plot_optimal_controls(non_opt_candidates=True)
    delta = GetInitialDelta(bounds, levels, N_Ti_variables)
    count = 0

    while count < 10 and np.all(delta != termination_range):
        finalcandidates = []

        for each_candidate in designer.optimal_candidates:
            optimal_candidate = each_candidate[1]
            boundsAroundCandidates, partitionlevels = GetNewPartition(bounds, optimal_candidate, N_Ti_variables, delta,
                                                                      sub_levels)
            newcandidates = designer.enumerate_candidates(boundsAroundCandidates, partitionlevels)
            finalcandidates.append(newcandidates)
            mat = np.concatenate(finalcandidates)

        delta = UpdateDelta(delta, termination_range, N_Ti_variables, sub_levels)

        mat = np.unique(mat, axis=0)

        designer.ti_controls_candidates = mat

        designer.initialize(verbose=2)

        designer.design_experiment(criterion=designer.d_opt_criterion, package="scipy", optimizer="SLSQP", )

        designer.get_optimal_candidates()

        if showallplots == True:
            designer.plot_optimal_controls(non_opt_candidates=True)
            designer.print_optimal_candidates(write=False)
        count += 1

    if showallplots == False:
        designer.plot_optimal_controls(non_opt_candidates=True)
        designer.print_optimal_candidates(write=False)

    designer.show_plots()

    # hello tests

# testing branch
