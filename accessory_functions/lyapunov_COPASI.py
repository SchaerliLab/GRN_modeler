# use the python interface for COPASI
# https://readthedocs.io/en/latest/API/html#module-task_timecourse
from basico import *

# load the model
model = load_model(model_name)

# get the original lyapunov settings
lyapunov_settings = get_task_settings(T.LYAPUNOV_EXPONENTS)
# change the problem and method
lyapunov_settings['problem'] = problem
lyapunov_settings['method'] = method

# use the settings
set_task_settings(T.LYAPUNOV_EXPONENTS,lyapunov_settings)

# calculate the exponents
exponents, _, _ = run_lyapunov(num_exponents=lyapunov_settings['problem']['ExponentNumber'],start_averaging_after=lyapunov_settings['problem']['TransientTime'])