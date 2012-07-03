import colorsys

def rainbow_colors(n = 8):
	return [colorsys.hsv_to_rgb(float(h)/n, 1.0, 1.0) for h in range(n)]

##########
# tests #
##########

def test_raibow_colors():
	r8 = rainbow_colors(8)
	assert(len(r8) == 8)
	assert(all([len(c) == 3 for c in r8]))
