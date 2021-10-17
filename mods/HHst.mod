
TITLE Stochastic Hodgkin and Huxley model & M-type potassium & T-and L-type Calcium channels incorporating channel noise .

COMMENT

Based on - Accurate and fast simulation of channel noise in conductance-based model neurons. Linaro, D., Storace, M., and Giugliano, M. 
Added: 
	Km T L channels
	fixed minor bugs and grouped variables into arrays 

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
} : end UNITS


NEURON {
    SUFFIX HHst
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
	USEION ca READ eca WRITE ica
	NONSPECIFIC_CURRENT ileak
	RANGE eleak, gleak,NFleak
    RANGE gnabar, gkbar
    RANGE gna, gk
	RANGE lm, lh, gl, glbar
	RANGE tm, th, gt, gtbar
	RANGE nm,  gkmbar
	RANGE m_exp, h_exp, n_exp, km_exp,lm_exp,tm_exp,lh_exp,th_exp
	
    GLOBAL gamma_na, gamma_k,gamma_km,gamma_l,gamma_t
	RANGE km_inf, tau_km,lm_inf,tm_inf,tau_lm,tau_tm,lh_inf,th_inf,tau_lh,tau_th	
    RANGE m_inf, h_inf, n_inf
    RANGE tau_m, tau_h, tau_n
    RANGE tau_zn,tau_zm,tau_zkm,tau_zt,tau_zl
    RANGE var_zn,var_zm,var_zkm,var_zt,var_zl
    RANGE noise_zn,noise_zm,noise_zkm,noise_zl,noise_zt
    RANGE Nna, Nk,Nkm,Nt,Nl
    GLOBAL seed    ,hslow,hfast
    RANGE mu_zn,mu_zm,mu_zkm,mu_zt,mu_zl
	:RANGE zl[4]
	GLOBAL vshift

	GLOBAL taukm,NF
    THREADSAFE
} : end NEURON


PARAMETER {
    gnabar  = 0.12   (S/cm2)
    gkbar   = 0.036  (S/cm2)
    glbar      = 0.0003 (S/cm2)  
    gtbar      = 0.0003 (S/cm2)  
	gkmbar = .002   	(S/cm2)
    gleak=0.00001		(S/cm2)
	eleak=-60 			(mV)
	NFleak=1
    gamma_na = 10  (pS)		: single channel sodium conductance
    gamma_k  = 10  (pS)		: single channel potassium conductance
    gamma_km  = 10  (pS)		: single channel potassium conductance
    gamma_t  = 10  (pS)		: single channel calcium conductance
    gamma_l  = 10  (pS)		: single channel calcium conductance
    seed = 1              : always use the same seed
	vshift=0				:Voltage shift of the recorded memebrane potential (to offset for liquid junction potential
	taukm=1					:speedup of Km channels
	NF=1					:Noise Factor (set to 0 to zero the noise part)
	hslow=100
	hfast=0.3
} : end PARAMETER


STATE {
    m h n nm tm th lm lh
    zn[3]
    zm[6]
	zkm
	zt[5]
	zl[5]
	zleak
} : end STATE


ASSIGNED {
    ina        (mA/cm2)
    il         (mA/cm2)
    it         (mA/cm2)
    ica         (mA/cm2)
	ikdr	(mA/cm2)
	ikm	 (mA/cm2)
	ik		(mA/cm2)
	ileak	(mA/cm2)
    gna        (S/cm2)
    gk         (S/cm2)
	gkm         (S/cm2)
	gl         (S/cm2)
	gt         (S/cm2)
    ena        (mV)
    ek         (mV)
    eca         (mV)
   
    dt         (ms)
    area       (um2)
    celsius    (degC)
    v          (mV)
        
    Nna  (1) : number of Na channels
    Nk   (1) : number of K channels
    Nkm   (1) : number of Km channels
    Nl   (1) : number of L channels
    Nt   (1) : number of T channels
    m_exp h_exp n_exp km_exp tm_exp lm_exp th_exp lh_exp
    m_inf h_inf n_inf km_inf tm_inf lm_inf th_inf lh_inf
    noise_zn[3]
    noise_zm[6]
	noise_zkm
	noise_zt[5]
	noise_zl[5]	
    var_zn[3]
    var_zm[6]
	var_zkm
	var_zt[5]
	var_zl[5]
    tau_m (ms) tau_h (ms) tau_n (ms) tau_km (ms) tau_lm (ms) tau_tm (ms) tau_th (ms) tau_lh (ms)
    tau_zn[3]
    tau_zm[6]
	tau_zkm
	tau_zt[5]
	tau_zl[5]
    mu_zn[3]
    mu_zm[6]
	mu_zkm
	mu_zt[5]
	mu_zl[5]
} : end ASSIGNED

INITIAL {
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   : area in um2 -> 1e-8*area in cm2; gnabar in S/cm2; gamma_na in pS -> 1e-12*gamma_na in S
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    Nkm = ceil(((1e-8)*area)*(gkmbar)/((1e-12)*gamma_km))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    Nt = ceil(((1e-8)*area)*(gtbar)/((1e-12)*gamma_t))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    Nl = ceil(((1e-8)*area)*(glbar)/((1e-12)*gamma_l))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    
    rates(v)
    m = m_inf
    h = h_inf
    n = n_inf
	nm=km_inf
	tm=tm_inf
	th=th_inf
	lm=lm_inf
	lh=lh_inf
	
	zn[0] = 0
	zn[1] = 0
	zn[2] = 0
	zn[3] = 0

    zm[0] = 0.
    zm[1] = 0.
    zm[2] = 0.
    zm[3] = 0.
    zm[4] = 0.
    zm[5] = 0.
    zm[6] = 0.
	
	zkm=0
	
	zt[0]=0
	zt[1]=0
	zt[2]=0
	zt[3]=0
	zt[4]=0

	zl[0]=0
	zl[1]=0
	zl[2]=0
	zl[3]=0
	zl[4]=0
	
	zleak=0
    set_seed(seed)
} : end INITIAL


BREAKPOINT {
    SOLVE states
    gna = gnabar * (m*m*m*h + (zm[0]+zm[1]+zm[2]+zm[3]+zm[4]+zm[5]+zm[6]))
:    if (gna < 0) {
		:gna = 0
:    }
    ina = gna * (v - ena)
    gk = gkbar * (n*n*n*n + (zn[0]+zn[1]+zn[2]+zn[3]))
:   if (gk < 0) {
:		gk = 0
:    }
    ikdr  = gk * (v - ek)
	gkm=gkmbar*(nm+zkm)
:	if (gkm < 0) {
:		gkm= 0
:    }
	ikm = gkm*(v-ek)
	ik=ikm+ikdr
	
    gt = gtbar * (tm*tm*th + (zt[0]+zt[1]+zt[2]+zt[3]+zt[4]))
:    if (gt < 0) {
:		gt = 0
:    }
    it  = gt * (v - eca)
    gl = glbar * (lm*lm*lh + (zl[0]+zl[1]+zl[2]+zl[3]+zl[4]))
:    if (gl < 0) {
:		gl = 0
:    }	
    il  = gl * (v - eca)
	ica=il+it
	ileak  = gleak*(1+zleak) * (v - eleak)
:	ileak  = (gleak) * (v - eleak)
} : end BREAKPOINT


PROCEDURE states() {
    rates(v+vshift)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
    n = n + n_exp * (n_inf - n)
	nm = nm + km_exp*(km_inf-nm)
    tm = tm + tm_exp * (tm_inf - tm)
    th = th + th_exp * (th_inf - th)
    lm = lm + lm_exp * (lm_inf - lm)
    lh = lh + lh_exp * (lh_inf - lh)

    VERBATIM
    return 0;
    ENDVERBATIM
} : end PROCEDURE states()


PROCEDURE rates(vm (mV)) { 
    LOCAL a,b,m3_inf,n4_inf,sum,one_minus_m,one_minus_h,one_minus_n,i
    
    UNITSOFF
    
    
 :NA m
	a =-.6 * vtrap((vm+30),-10)	
	b = 20 * (exp((-1*(vm+55))/18))
	tau_m = 1 / (a + b)
	m_inf = a * tau_m

    one_minus_m = 1. - m_inf
    m3_inf = m_inf*m_inf*m_inf

:NA h
	a = 0.4 * (exp((-1*(vm+50))/20))
	b = 6 / ( 1 + exp(-0.1 *(vm+20)))
	:b = 6 / ( 1 + exp(-(vm-50)/5))
	:tau_h = 1 / (a + b)
	:h_inf = a * tau_h
	:tau_h=1/(a+6 / ( 1 + exp(-(vm-50)/5)))
	:tau_h=hslow/(1+exp((vm+30)/4))+hfast
	tau_h=hslow/((1+exp((vm+30)/4))+(exp(-(vm+50)/2)))+hfast
	h_inf= 1/(1 + exp((vm + 44)/4))
	:h_inf = a / (a+b)

	one_minus_h = 1. - h_inf
    
    tau_zm[0] = tau_h
    tau_zm[1] = tau_m
    tau_zm[2] = tau_m/2
    tau_zm[3] = tau_m/3
    tau_zm[4] = tau_m*tau_h/(tau_m+tau_h)
    tau_zm[5] = tau_m*tau_h/(tau_m+2*tau_h)
    tau_zm[6] = tau_m*tau_h/(tau_m+3*tau_h)
    var_zm[0] = 1.0 / numchan(Nna) * m3_inf*m3_inf*h_inf * one_minus_h
    var_zm[1] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m
    var_zm[2] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m
    var_zm[3] = 1.0 / numchan(Nna) * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m
    var_zm[4] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h
    var_zm[5] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h
    var_zm[6] = 1.0 / numchan(Nna) * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h
	
	FROM i=0 TO 6 {
		mu_zm[i] = exp(-dt/tau_zm[i])
        noise_zm[i] = sqrt(var_zm[i] * (1-mu_zm[i]*mu_zm[i])) * normrand(0,NF)
		zm[i] = zm[i]*mu_zm[i] + noise_zm[i]
    }
  
   
:K n (non-inactivating, delayed rectifier)
	a = -0.02 * vtrap((vm+40),-10)
	b = 0.4 * (exp((-1*(vm + 50))/80))
	tau_n = 1 / (a + b)
	n_inf = a * tau_n
    one_minus_n = 1. - n_inf
    n4_inf = n_inf * n_inf * n_inf * n_inf
    tau_zn[0] = tau_n
    tau_zn[1] = tau_n/2
    tau_zn[2] = tau_n/3
    tau_zn[3] = tau_n/4
    var_zn[0] = 4.0/numchan(Nk) * n4_inf*n_inf*n_inf*n_inf * one_minus_n
    var_zn[1] = 6.0/numchan(Nk) * n4_inf*n_inf*n_inf * one_minus_n*one_minus_n
    var_zn[2] = 4.0/numchan(Nk) * n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n
    var_zn[3] = 1.0/numchan(Nk) * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n

	FROM i=0 TO 3 {
		mu_zn[i] = exp(-dt/tau_zn[i])
        noise_zn[i] = sqrt(var_zn[i] * (1-mu_zn[i]*mu_zn[i])) * normrand(0,NF)
		zn[i] = zn[i]*mu_zn[i] + noise_zn[i]
    }
:Km
    a = -.001/taukm * vtrap((vm+30),-9)
    b =.001/taukm * vtrap((vm+30),9) 
    tau_km = 1/(a+b)
	km_inf = a*tau_km
	tau_zkm = tau_km
	var_zkm=km_inf*(1-km_inf)/numchan(Nkm)
	mu_zkm = exp(-dt/tau_zkm)
    noise_zkm = sqrt(var_zkm * (1-mu_zkm*mu_zkm)) * normrand(0,NF)
	zkm = zkm*mu_zkm + noise_zkm
	
:L Ca (m)
	a = 0.055*vtrap(-(vm+27),3.8): (-27 - vm)/(exp((-27-vm)/3.8) - 1)
	b = 0.94*exp((-75-vm)/17)
    tau_lm = 1/(a+b)
	lm_inf = a*tau_lm
:L Ca (h)	
	a = 0.000457*exp((-13-vm)/50)
	b = 0.0065/(exp((-vm-15)/28) + 1)	
    tau_lh = 1/(a+b)
	lh_inf = a*tau_lh

    tau_zl[0] = tau_lm*tau_lh/(tau_lm+2*tau_lh)
    tau_zl[1] = tau_lm*tau_lh/(tau_lm+tau_lh)
    tau_zl[2] = tau_lh
    tau_zl[3] = tau_lm/2
    tau_zl[4] = tau_lm

    var_zl[0] = 1/numchan(Nl) * lm_inf^2*lh_inf*(1-lm_inf)^2*(1-lh_inf)
    var_zl[1] = 2/numchan(Nl) * lm_inf^3*lh_inf*(1-lm_inf)*(1-lh_inf)
    var_zl[2] = 1/numchan(Nl) * lm_inf^4*lh_inf*(1-lh_inf)
    var_zl[3] = 1/numchan(Nl) * lm_inf^2*lh_inf^2*(1-lm_inf)^2
    var_zl[4] = 2/numchan(Nl) * lm_inf^3*lh_inf^2*(1-lm_inf)

	FROM i=0 TO 4 {
		mu_zl[i] = exp(-dt/tau_zl[i])
        noise_zl[i] = sqrt(var_zl[i] * (1-mu_zl[i]*mu_zl[i])) * normrand(0,NF)
		zl[i] = zl[i]*mu_zl[i] + noise_zl[i]
    }	

:T Ca (m)
    tau_tm =( 4 / ( exp((vm+25)/20) + exp(-(vm+100)/15) ) ) 
	tm_inf = 1 / ( 1 + exp(-(vm+50)/7.4) )
:T Ca (h)	

    tau_th = ( 86/ ( exp((vm+46)/4) + exp(-(vm+405)/50) ) ) 
	th_inf = 1.0 / ( 1 + exp((vm+78)/5) )

    tau_zt[0] = tau_tm*tau_th/(tau_tm+2*tau_th)
    tau_zt[1] = tau_tm*tau_th/(tau_tm+tau_th)
    tau_zt[2] = tau_th
    tau_zt[3] = tau_tm/2
    tau_zt[4] = tau_tm

    var_zt[0] = 1/numchan(Nt) * tm_inf^2*th_inf*(1-tm_inf)^2*(1-th_inf)
    var_zt[1] = 2/numchan(Nt) * tm_inf^3*th_inf*(1-tm_inf)*(1-th_inf)
    var_zt[2] = 1/numchan(Nt) * tm_inf^4*th_inf*(1-th_inf)
    var_zt[3] = 1/numchan(Nt) * tm_inf^2*th_inf^2*(1-tm_inf)^2
    var_zt[4] = 2/numchan(Nt) * tm_inf^3*th_inf^2*(1-tm_inf)

	FROM i=0 TO 4 {
		mu_zt[i] = exp(-dt/tau_zt[i])
        noise_zt[i] = sqrt(var_zt[i] * (1-mu_zt[i]*mu_zt[i])) * normrand(0,NF)
		zt[i] = zt[i]*mu_zt[i] + noise_zt[i]
    }	
	
	zleak=normrand(0,NFleak)
	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	n_exp = 1 - exp(-dt/tau_n)
	km_exp= 1 - exp(-dt/tau_km)	
	lm_exp= 1 - exp(-dt/tau_lm)	
	lh_exp= 1 - exp(-dt/tau_lh)	
	tm_exp= 1 - exp(-dt/tau_tm)	
	th_exp= 1 - exp(-dt/tau_th)	
	
   UNITSON
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(exp(x/y) - 1) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}
FUNCTION mulnoise(mean,sd,power){
	LOCAL i,avg
	avg=1
	FROM i=1 TO power {
		avg=avg*normrand(mean,sd)
	}
	mulnoise=avg
}

FUNCTION numchan(Nchannels){
	if (Nchannels>0){
		numchan=(Nchannels)
	}else{
		numchan=1
	}
}