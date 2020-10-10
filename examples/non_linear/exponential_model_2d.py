from pydex.core.designer import Designer
import numpy as np
import Pydex_IterativeSolver as Solve

def simulate(ti_controls, model_parameters):
    return np.array([
        model_parameters[0] +
        model_parameters[1] * np.exp(model_parameters[2] * ti_controls[0]) +
        model_parameters[3] * np.exp(model_parameters[4] * ti_controls[1])
    ])

designer_1 = Designer()
designer_1.simulate = simulate

mp = np.array([1, 2, 2, 10, 2])
designer_1.model_parameters = mp

bounds = [[-1, 1], [-1, 1], ]
#levels = [6, 6, ]
termination = [0, 0]

Solve.IterativeSolver(designer=designer_1, bounds=bounds, package="scipy", optimizer="SLSQP", showallplots=True, termination_range = termination)



#reso = 21j
#tic1, tic2 = np.mgrid[-1:1:reso, -1:1:reso]
#tic = np.array([tic1.flatten(), tic2.flatten()]).T
#designer.ti_controls_candidates = tic

#designer.initialize(verbose=2)

#criterion = designer.a_opt_criterion
#designer.design_experiment(criterion, write=False)
#designer.print_optimal_candidates()
#designer.plot_controls()
