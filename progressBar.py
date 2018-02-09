import sys

#val is the current "record number" the program is on.
#endv is the total number of records that is being processed.
#pad is the number of character spaces that should be saved for val and endv

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
