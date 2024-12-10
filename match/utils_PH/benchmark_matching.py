import psutil
import datetime

def benchmark_out(context_msg, lapsed_time):
    print(context_msg)
    print('Total elapsed time to point: ' + str(lapsed_time))
    print('Absolute memory usage: ' + str(psutil.virtual_memory()[3]/(10**9)) + ' GB')
    print('')

    return None
