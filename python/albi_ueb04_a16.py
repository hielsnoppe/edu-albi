data = {
	"example": {
			"states": ('Healthy', 'Fever'),
			"observations": ('normal', 'cold', 'dizzy'),
			"start_probability": {'Healthy': 0.6, 'Fever': 0.4},
			"transition_probability": {
				'Healthy' : {'Healthy': 0.7, 'Fever': 0.3},
				'Fever' : {'Healthy': 0.4, 'Fever': 0.6},
				},
			"emission_probability": {
				'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
				'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
				}
		},
	"a15": {
			"states": ("1", "2", "3", "4", "E"),
			"observations": ("C", "T", "G", "END"),
			"start_probability": {"1": 0.5, "2": 0, "3": 0.5, "4": 0, "E": 0},
			"transition_probability": {
				"1": {"Begin": 0, "1": 0, "2": 1.0, "3": 0, "4": 0, "E": 0},
				"2": {"Begin": 0, "1": 0, "2": 0.5, "3": 0.3, "4": 0, "E": 0.2},
				"3": {"Begin": 0, "1": 0, "2": 0, "3": 0, "4": 1.0, "E": 0},
				"4": {"Begin": 0, "1": 0.3, "2": 0, "3": 0, "4": 0.5, "E": 0.2},
				"E": {"Begin": 0, "1": 0, "2": 0, "3": 0, "4": 0, "E": 0}
				},
			"emission_probability": {
				"1": {"A": 0.2, "C": 0.1, "T": 0.2, "G": 0.5, "END": 0},
				"2": {"A": 0.5, "C": 0.2, "T": 0.1, "G": 0.2, "END": 0},
				"3": {"A": 0.2, "C": 0.5, "T": 0.2, "G": 0.1, "END": 0},
				"4": {"A": 0.1, "C": 0.2, "T": 0.5, "G": 0.2, "END": 0},
				"E": {"A": 0, "C": 0, "T": 0, "G": 0, "END": 1.0}
				}
		},
	"a16": {
			"states": ("F", "L"),
			"observations": tuple("315116246446644245311321631164152133625144543631656626566666\
651166453132651245636664631636663162326455236266666625151631\
222555441666566563564324364131513465146353411126414626253356\
366163666466232534413661661163252562462255265252266435353336\
233121625364414432335163243633665562466662632666612355245242"),
			"start_probability": {"F": 0.5, "L": 0.5},
			"transition_probability": {
				"F": {"F": 0.9, "L": 0.1},
				"L": {"F": 0.05, "L": 0.95}
				},
			"emission_probability": {
				"F": {"1": 1.0/6.0, "2": 1.0/6.0, "3": 1.0/6.0, "4": 1.0/6.0, "5": 1.0/6.0, "6": 1.0/6.0},
				"L": {"1": 0.1, "2": 0.1, "3": 0.1, "4": 0.1, "5": 0.1, "6": 0.5}
				}
		}
	}

# Helps visualize the steps of Viterbi.
def print_dptable(V):
    print "    ",
    for i in range(len(V)): print "%7d" % i,
    print
 
    for y in V[0].keys():
        print "%.5s: " % y,
        for t in range(len(V)):
            print "%.7s" % ("%f" % V[t][y]),
        print
 
def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]
 
    # Run Viterbi for t > 0
    for t in range(1,len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = max([(V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states])
            V[t][y] = prob
            newpath[y] = path[state] + [y]
 
        # Don't need to remember the old paths
        path = newpath
 
    print_dptable(V)
    (prob, state) = max([(V[len(obs) - 1][y], y) for y in states])
    return (prob, path[state])

def a15():
	return viterbi(data["a15"]["observations"],
                   data["a15"]["states"],
                   data["a15"]["start_probability"],
                   data["a15"]["transition_probability"],
                   data["a15"]["emission_probability"])
def a16(sub = 1):
	observations = data["a16"]["observations"]
	states = data["a16"]["states"]
	if sub == 2:
		observations = data["a16"]["observations"][::-1]
	if sub == 3:
		states = reversed(states)
	return viterbi(observations,
                   states,
                   data["a16"]["start_probability"],
                   data["a16"]["transition_probability"],
                   data["a16"]["emission_probability"])
def example():
    return viterbi(data["example"]["observations"],
                   data["example"]["states"],
                   data["example"]["start_probability"],
                   data["example"]["transition_probability"],
                   data["example"]["emission_probability"])

print("".join(a16(3)[1]))
