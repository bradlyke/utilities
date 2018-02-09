#This just prints a time code line. Call it with an input string so that blocks
#of code have a start marker if you have to get up and leave. In the following
#example istring would be 'Making Hash Array'
#Example output:
#'Making Hash Array: 2017-11-29 | 14:20:02'
#Typically progressBar will be called in the loop that runs after this so the
#next line will be a filling bar.

#Call with tmark.tm(input_string)
#If you want to also use a timer to benchmark blocks of code use "timer=True"
#then call the function as <var> = tmark.tm(istring,timer=True). <var> then will
#be the start time and will be returned to the original calling program.

import time

def tm(istring,timer=False):
    t0 = time.time()
    t_start = time.strftime('%Y-%m-%d | %H:%M:%S')
    print('\n{}: {}'.format(istring,t_start))

    if timer == True:
        return t0
