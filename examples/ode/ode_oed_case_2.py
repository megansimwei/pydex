import numpy as np
from pyomo import dae as pod
from pyomo import environ as po

from pydex.core.designer import Designer


def simulate(model, simulator, ti_controls, sampling_times, model_parameters):
    """ fixing the control variables """
    # time-invariant
    model.theta_0.fix(model_parameters[0])
    model.theta_1.fix(model_parameters[1])
    model.alpha_a.fix(model_parameters[2])
    model.alpha_b.fix(0)
    model.nu.fix(model_parameters[3])
    # model.nu.fix(1)

    model.tau.fix(max(sampling_times))
    model.ca[0].fix(ti_controls[0])
    model.cb[0].fix(0)
    model.temp.fix(ti_controls[1])
    # no time-varying control for this example

    """ ensuring pyomo returns state values at given sampling times """
    model.t.initialize = np.array(sampling_times) / model.tau.value
    model.t.order_dict = {}  # to suppress pyomo warnings for duplicate elements
    model.t._constructed = False  # needed so we can re-initialize the continuous set
    model.t._data = {}
    model.t._fe = []
    model.t.value_list = []
    model.t.value = []
    model.t.construct()  # line that re-initializes the continuous set

    """ simulating """
    simulator.simulate(integrator='idas')
    simulator.initialize_model()

    """" extracting results and returning it in appropriate format """
    normalized_sampling_times = sampling_times / model.tau.value
    ca = np.array([model.ca[t].value for t in normalized_sampling_times])
    cb = np.array([model.cb[t].value for t in normalized_sampling_times])

    return np.array([ca, cb]).T


def create_model():
    """ defining the model """
    model = po.ConcreteModel()
    model.t = pod.ContinuousSet(bounds=(0, 1))
    model.tau = po.Var()

    model.temp = po.Var()

    model.ca = po.Var(model.t, bounds=(0, 50))
    model.cb = po.Var(model.t, bounds=(0, 50))
    model.dca_dt = pod.DerivativeVar(model.ca, wrt=model.t)
    model.dcb_dt = pod.DerivativeVar(model.cb, wrt=model.t)

    model.theta_0 = po.Var()  # model parameters
    model.theta_1 = po.Var()
    model.alpha_a = po.Var()
    model.alpha_b = po.Var()
    model.nu = po.Var()

    def _material_balance_a(m, t):
        k = po.exp(m.theta_0 + m.theta_1 * (m.temp - 273.15) / m.temp)
        return m.dca_dt[t] / m.tau == - k * (m.ca[t] ** model.alpha_a) * (
                model.cb[t] ** model.alpha_b)

    model.material_balance_a = po.Constraint(model.t, rule=_material_balance_a)

    def _material_balance_b(m, t):
        k = po.exp(m.theta_0 + m.theta_1 * (m.temp - 273.15) / m.temp)
        return m.dcb_dt[t] / m.tau == m.nu * k * (m.ca[t] ** model.alpha_a) * (
                model.cb[t] ** model.alpha_b)

    model.material_balance_b = po.Constraint(model.t, rule=_material_balance_b)

    return model


pre_exp_constant = 0.1
activ_energy = 5000
theta_0 = np.log(pre_exp_constant) - activ_energy / (8.314159 * 273.15)
theta_1 = activ_energy / (8.314159 * 273.15)

""" create a pyomo model """
model_1 = create_model()
simulator_1 = pod.Simulator(model_1, package='casadi')

""" create a designer """
designer_1 = Designer()

""" pass pyomo model and simulator to designer """
designer_1.model = model_1
designer_1.simulator = simulator_1

""" overwrite the designer's simulate function """
designer_1.simulate = simulate

""" specifying nominal model parameter """
theta_nom = np.array([theta_0, theta_1, 1, 0.5])  # value of theta_0, theta_1, alpha_a, nu
designer_1.model_parameters = theta_nom  # assigning it to the designer's theta

""" creating experimental candidates, here, it is generated as a grid """
n_s_times = 200  # number of equally-spaced sampling time candidates
n_c = 20 ** 2  # grid resolution of control candidates generated

# defining sampling time candidates
tau_upper = 200
tau_lower = 0
spt_candidates = np.array([np.linspace(tau_lower, tau_upper, n_s_times)
                           for _ in range(n_c)])

# specifying bounds for the grid
Ca0_lower = 1
Ca0_upper = 5
temp_lower = 273.15
temp_upper = 273.15 + 50
# creating the grid, just some numpy syntax for grid creation
Ca0_cand, temp_cand = np.mgrid[Ca0_lower:Ca0_upper:complex(0, np.sqrt(n_c)),
                      temp_lower:temp_upper:complex(0, np.sqrt(n_c))]
Ca0_cand = Ca0_cand.flatten()
temp_cand = temp_cand.flatten()
tic_candidates = np.array([Ca0_cand, temp_cand]).T

""" passing the experimental candidates to the designer """
designer_1.ti_controls_candidates = tic_candidates
designer_1.sampling_times_candidates = spt_candidates

"""
Specify measurable states:
A list or array with column numbers where the measurable states are returned in the simulate 
function. Optional, if un-specified assume all responses (from simulate function) measurable
"""
designer_1.measurable_responses = [0, 1]
# designer_1.measurable_responses = [0]
# designer_1.measurable_responses = [1]

""" optional information for plotting purposes, if unspecified empty axes titles """
designer_1.candidate_names = np.array(["Candidate {:d}".format(i + 1)
                                       for i, _ in enumerate(tic_candidates)])

""" optional information for estimability study """
designer_1.responses_scales = np.array([1, 1])

""" initializing designer """
designer_1.initialize(verbose=2)  # 0: silent, 1: overview, 2: detail

""" option to save current designer state """
# designer_1.save_state()

""" do an estimability study """
# designer_1.estimability_study(save_sensitivities=True)

""" (optional) plotting attributes """
designer_1.response_names = ["c_A", "c_B"]
designer_1.model_parameter_names = ["\\theta_0", "\\theta_1", "\\alpha", "\\nu"]

""" D-optimal design """
criterion = designer_1.d_opt_criterion
# criterion = designer_1.a_opt_criterion
# criterion = designer_1.e_opt_criterion

result = designer_1.design_experiment(criterion=criterion, optimize_sampling_times=False,
                                      write=False, fd_jac=False, package="cvxpy")

designer_1.print_optimal_candidates()
designer_1.plot_optimal_efforts()

designer_1.plot_optimal_predictions()
designer_1.plot_optimal_sensitivities()
