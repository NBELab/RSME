NEURON {
	POINT_PROCESS BcSyn
	RANGE tau, e, i
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	e = 0	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
}

STATE {
	g (uS)
}

INITIAL {
	g=0
}

BREAKPOINT {
	i = g*(v - e)
}


