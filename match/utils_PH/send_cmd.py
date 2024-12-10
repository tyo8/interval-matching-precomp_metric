def send_cmd_windows(cmd) : # cmd is a string
    import subprocess
    'dont forget stdout = good output,  stderr = error like outputs e.g. help output for windows'
    out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('UTF-8').split('\n')
    # choosing first good output
    return out

def send_cmd_linux(cmd, spro_flag=False) :
    if spro_flag:
        import subprocess
        out = subprocess.Popen(
                cmd.split(' '), stdin=subprocess.PIPE, stdout=subprocess.PIPE
                ).stdout.read().decode().split('\n')
    else:
        import os
        out = os.popen(cmd).read().split('\n')

    return out

