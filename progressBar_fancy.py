"""Creates a FANCY progress bar in the terminal that fills over one line only.

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

	The arrow type shows the "speed health" of the iterations, where _ means
	that each iteration is taking more than 10 seconds, - means each iteration
	is on the order of 1 second, and = means that each iteration is less than
	1 second.

"""

import sys
from colorama import Fore
from colorama import Back
from colorama import Style


def pbar(val,endv,bar_length=40,speed_tracker=0,estimator=False):
	pad = len(str(endv))
	num_comma = int(pad/3)
	pad2 = pad + num_comma
	prct = float(val + 1) / float(endv)
	if speed_tracker < 1.0:
		arrow = '=' * int(round(prct * bar_length) - 1) + '>'
	elif ((speed_tracker >= 1.0)&(speed_tracker < 10.0)):
		arrow = '-' * int(round(prct * bar_length) - 1) + '>'
	else:
		arrow = '_' * int(round(prct * bar_length) - 1) + '>'
	spcs = ' ' * int(bar_length - len(arrow))
	vals = val + 1

	#Time Estimation Block
	iter_left = endv - vals
	time_left = speed_tracker * iter_left
	hours_left = int(time_left / 3600.0)
	minutes_left = int((time_left % 3600.0)/60.0)
	seconds_left = int((time_left % 3600.0)%60.0)

	st_prct = int(round(prct * 100))
	ars = arrow + spcs
	if ((estimator==False)&(prct <= 0.33)):
		st = '\rRecord: {:{},d}/{:{},d} | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(vals,
			pad2,endv,pad2,Fore.RED,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)
	elif ((estimator==False)&(prct > 0.33)&(prct <= 0.67)):
		st = '\rRecord: {:{},d}/{:{},d} | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(vals,
			pad2,endv,pad2,Fore.YELLOW,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)
	elif ((estimator==False)&(prct > 0.67)):
		st = '\rRecord: {:{},d}/{:{},d} | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(vals,
			pad2,endv,pad2,Fore.GREEN,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)
	elif ((estimator==True)&(prct <= 0.33)):
		st = '\rEst. Time Left:{:2d}h {:2d}m {:2d}s | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(hours_left,
			minutes_left,seconds_left,Fore.RED,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)
	elif ((estimator==True)&(prct > 0.33)&(prct<=0.67)):
		st = '\rEst. Time Left:{:2d}h {:2d}m {:2d}s | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(hours_left,
			minutes_left,seconds_left,Fore.YELLOW,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)
	else:
		st = '\rEst. Time Left:{:2d}h {:2d}m {:2d}s | Percent: {}[{:{}}]{} {}{:{}}%{}'.format(hours_left,
			minutes_left,seconds_left,Fore.GREEN,ars,bar_length,Style.RESET_ALL,Style.BRIGHT,
			st_prct,3,Style.RESET_ALL)

	sys.stdout.write(st)
	sys.stdout.flush()
