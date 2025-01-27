# use the python interface for COPASI
# https://readthedocs.io/en/latest/API/html#module-task_timecourse
from basico import *

# load the model
model = load_model(model_name)

# transform the maximum runtime to maximum solver step number
if MaximumWallClock == float('inf'):
    MaximumWallClock=2**31 - 1;
else:
    # calculate solver speed
    import time

    steps = 100  # Number of steps for the sample simulation

    # Save the current stderr (we will exceed the number of steps)
    old_stderr = sys.stderr
    # Redirect stderr to os.devnull
    sys.stderr = open(os.devnull, 'w')

    # short test simulation
    try:
        start_time = time.time()
        tmp = run_time_course_with_output(output_selection=StatesTolog,a_tol=abstol,r_tol=reltol,start_time=StartTime,duration=StopTime,max_steps=steps,automatic=True,method=solver)
        end_time = time.time()
    finally:
        # Restore the original stderr
        sys.stderr = old_stderr

    # calculate speed
    total_time = end_time - start_time
    time_per_step = total_time / steps

    # maximum number of steps
    MaximumWallClock = MaximumWallClock/time_per_step;


# run simulation
if not times:
    simdata = run_time_course_with_output(output_selection=StatesTolog,a_tol=abstol,r_tol=reltol,start_time=StartTime,duration=StopTime,max_steps=MaximumWallClock,automatic=True,method=solver)
else:
    simdata = run_time_course_with_output(output_selection=StatesTolog,a_tol=abstol,r_tol=reltol,start_time=StartTime,duration=StopTime,max_steps=MaximumWallClock,values=times,automatic=False,method=solver)
