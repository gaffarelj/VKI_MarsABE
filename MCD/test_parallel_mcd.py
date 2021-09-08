from parallel_mcd import parallel_mcd as PMCD
import time
import numpy as np

# Test the module without loading the MCD in parallel for each Month
mcd_simple = PMCD(default_inputs=True, load_parallel=False)
mcd_simple.call(27, 11, -15, 35, 250e3, print_results=True)

input("Press ENTER to load modules in parallel next.")
mcd = PMCD(True, True)

# Test one call of the MCD and print the results
mcd.call(27, 11, -15, 35, 250e3, print_results=True)

# Test thousands of calls, changing the inputs every time
N = int(1e6)
t0 = time.time()
for i in range(N):
    date = np.random.uniform(0, 330)
    localtime = np.random.uniform(0, 24)
    lat = np.random.uniform(-90, 90)
    lon = np.random.uniform(-180, 180)
    h = np.random.uniform(200e3, 300e3)
    mcd.call(date, localtime, lat, lon, h)
call_time = time.time() - t0
print("Calling the MCD with different times took %.5f seconds." % (call_time))
print("This is an average of %.5f milliseconds per call" % (call_time*1e3/N)) # takes 0.05 ms per call using low res, and 0.114 ms per call using high res