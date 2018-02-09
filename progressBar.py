"""Creates a progress bar in the terminal that fills over one line only.

Parameters
----------
val : :class:'int'
	The current record number that the calling for loop is on.
endv : :class:'int'
	The total number of records that is being processed.
bar_length : :class:'int'
	The length, in characters, that the bar will set aside to fill.
        Originally designed for a terminal 96 characters wide.

Returns
----------
:class:'str'
	A string showing the current progress.

"""

import sys

#BENCHMARK NOTES: Running this program took 18.98 seconds for 1,000,000 iterations.
#It took 1.8977e-5 seconds per iteration.
#Loops will add ~2 seconds per 100,000 calls of this function.

def pbar(val,endv,bar_length=40):
    pad = len(str(endv))
    num_comma = int(pad/3)
    pad2 = pad + num_comma
    prct = float(val) / float(endv)
    arrow = '=' * int(round(prct * bar_length) - 1) + '>'
    spcs = ' ' * int(bar_length - len(arrow))

    vals = val + 1
    st_prct = int(round(prct * 100))
    ars = arrow + spcs
    st = '\rRecord: {:{},d}/{:{},d} | Percent: [{:{}}] {:{}}%'.format(vals,pad2,endv,pad2,ars,bar_length,st_prct,3)
    sys.stdout.write(st)
    sys.stdout.flush()
