def get_time_elapse(start_time):
	import time

	time_elapse = time.time() - start_time
	hours = time_elapse // 3600
	rest_t = time_elapse % 3600
	minutes = rest_t // 60
	seconds = rest_t % 60

	return hours, minutes, seconds

def report_time_elapse(hours, minutes, seconds):
	print("Time elapse: %d : %d : %.2f" % (hours, minutes, seconds), flush = True)
