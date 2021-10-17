/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__caq
#define _nrn_initial _nrn_initial__caq
#define nrn_cur _nrn_cur__caq
#define _nrn_current _nrn_current__caq
#define nrn_jacob _nrn_jacob__caq
#define nrn_state _nrn_state__caq
#define _net_receive _net_receive__caq 
#define _f_settables _f_settables__caq 
#define settables settables__caq 
#define states states__caq 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define pcaqbar _p[0]
#define ica _p[1]
#define m _p[2]
#define eca _p[3]
#define minf _p[4]
#define cai _p[5]
#define cao _p[6]
#define Dm _p[7]
#define _g _p[8]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
#define mu	*_ppvar[4]._pval
#define _p_mu	_ppvar[4]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  4;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_settables(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_caq", _hoc_setdata,
 "efun_caq", _hoc_efun,
 "ghk_caq", _hoc_ghk,
 "settables_caq", _hoc_settables,
 0, 0
};
#define efun efun_caq
#define ghk ghk_caq
 extern double efun( double );
 extern double ghk( double , double , double );
 /* declare global and static user variables */
#define mshift mshift_caq
 double mshift = 0;
#define mtau mtau_caq
 double mtau = 1.13;
#define mslope mslope_caq
 double mslope = -6.6;
#define mvhalf mvhalf_caq
 double mvhalf = -9;
#define qfact qfact_caq
 double qfact = 3;
#define usetable usetable_caq
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_caq", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mvhalf_caq", "mV",
 "mslope_caq", "mV",
 "mtau_caq", "ms",
 "mshift_caq", "mV",
 "pcaqbar_caq", "cm/s",
 "ica_caq", "mA/cm2",
 "mu_caq", "1",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mvhalf_caq", &mvhalf_caq,
 "mslope_caq", &mslope_caq,
 "mtau_caq", &mtau_caq,
 "mshift_caq", &mshift_caq,
 "qfact_caq", &qfact_caq,
 "usetable_caq", &usetable_caq,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"caq",
 "pcaqbar_caq",
 0,
 "ica_caq",
 0,
 "m_caq",
 0,
 "mu_caq",
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	pcaqbar = 6e-06;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _caq_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "pointer");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 caq /Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/x86_64/caq.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
 static double *_t_minf;
static int _reset;
static char *modelname = "Q-type (not P) calcium channel for nucleus accumbens neuron ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(double);
static int settables(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(double);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   settables ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / ( mtau / qfact ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 settables ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mtau / qfact ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   settables ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mtau / qfact ))))*(- ( ( ( minf ) ) / ( mtau / qfact ) ) / ( ( ( ( - 1.0 ) ) ) / ( mtau / qfact ) ) - m) ;
   }
  return 0;
}
 static double _mfac_settables, _tmin_settables;
 static void _check_settables();
 static void _check_settables() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mshift;
  if (!usetable) {return;}
  if (_sav_mshift != mshift) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/201.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 202; _x += _dx, _i++) {
    _f_settables(_x);
    _t_minf[_i] = minf;
   }
   _sav_mshift = mshift;
  }
 }

 static int settables(double _lv){ _check_settables();
 _n_settables(_lv);
 return 0;
 }

 static void _n_settables(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
 if (isnan(_xi)) {
  minf = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 return; }
 if (_xi >= 201.) {
 minf = _t_minf[201];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 }

 
static int  _f_settables (  double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( ( _lv - mvhalf - mshift ) / mslope ) ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
    _r = 1.;
 settables (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk (  double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lz , _leci , _leco ;
 _lz = ( 1e-3 ) * 2.0 * FARADAY * _lv / ( R * ( celsius + 273.15 ) ) ;
   _leco = _lco * efun ( _threadargscomma_ _lz ) ;
   _leci = _lci * efun ( _threadargscomma_ - _lz ) ;
   _lghk = ( .001 ) * 2.0 * FARADAY * ( _leci - _leco ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   _r =  ghk (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  m = m0;
 {
   settables ( _threadargscomma_ v ) ;
   m = minf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ica = ghk ( _threadargscomma_ v , cai , cao ) * pcaqbar * m * m * ( 1.0 - ( mu - 1.0 ) * 0.5 ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
 { error =  states();
 if(error){fprintf(stderr,"at line 50 in file caq.mod:\n    SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
   _t_minf = makevector(202*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/mods/caq.mod";
static const char* nmodl_file_text = 
  "TITLE Q-type (not P) calcium channel for nucleus accumbens neuron \n"
  ": see comments at end of file\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(S) = (siemens)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	FARADAY = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX caq\n"
  "	USEION ca READ cai, cao WRITE ica\n"
  "	RANGE pcaqbar, ica\n"
  "		POINTER mu\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	pcaqbar = 6.0e-6(cm/s)		: vh = -100 mV, 120 ms pulse to 0 mV\n"
  "\n"
  "	mvhalf = -9.0	(mV)		: Churchill 1998, fig 5\n"
  "	mslope = -6.6	(mV)		: Churchill 1998, fig 5\n"
  "	mtau = 1.13	(ms)			: Randall 1995, fig 13\n"
  "	mshift = 0	(mV)\n"
  "	\n"
  "	qfact = 3					: m recorded at 22 C\n"
  "}\n"
  "\n"
  "ASSIGNED { \n"
  "    v 		(mV)\n"
  "    eca		(mV)\n"
  "    ica 	(mA/cm2)\n"
  "    minf\n"
  "\n"
  "    celsius	(degC)\n"
  "    cai		(mM)\n"
  "    cao		(mM)\n"
  "\n"
  "		mu (1)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    m\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE states METHOD cnexp\n"
  "    ica  = ghk(v,cai,cao) * pcaqbar * m * m *(1- (mu-1)*0.5)   : Kasai 92, Brown 93\n"
  "}						\n"
  "\n"
  "INITIAL {\n"
  "    settables(v)\n"
  "    m = minf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  "    settables(v)\n"
  "    m' = (minf - m) / (mtau/qfact)\n"
  "}\n"
  "\n"
  "PROCEDURE settables( v (mV) ) {\n"
  "	TABLE minf DEPEND mshift\n"
  "        FROM -100 TO 100 WITH 201\n"
  "\n"
  "		minf = 1  /  ( 1 + exp( (v-mvhalf-mshift) / mslope) )\n"
  "}\n"
  "\n"
  "\n"
  "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"
  "\n"
  ": ghk() borrowed from cachan.mod share file in Neuron\n"
  "FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {\n"
  "	LOCAL z, eci, eco\n"
  "	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))\n"
  "	eco = co*efun(z)\n"
  "	eci = ci*efun(-z)\n"
  "	:high cao charge moves inward\n"
  "	:negative potential charge moves inward\n"
  "	ghk = (.001)*2*FARADAY*(eci - eco)\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "	if (fabs(z) < 1e-4) {\n"
  "		efun = 1 - z/2\n"
  "	}else{\n"
  "		efun = z/(exp(z) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "Brown AM, Schwindt PC, Crill WE (1993) Voltage dependence and activation\n"
  "kinetics of pharmacologically defined components of the high-threshold\n"
  "calcium current in rat neocortical neurons. J Neurophysiol 70:1530-1543.\n"
  "\n"
  "Churchill D, Macvicar BA (1998) Biophysical and pharmacological\n"
  "characterization of voltage-dependent Ca2+ channels in neurons isolated\n"
  "from rat nucleus accumbens. J Neurophysiol 79:635-647.\n"
  "\n"
  "Kasai H, Neher E (1992) Dihydropyridine-sensitive and\n"
  "omega-conotoxin-sensitive calcium channels in a mammalian\n"
  "neuroblastoma-glioma cell line. J Physiol 448:161-188.\n"
  "\n"
  "Mermelstein PG, Foehring RC, Tkatch T, Song WJ, Baranauskas G, Surmeier\n"
  "DJ (1999) Properties of Q-type calcium channels in neostriatal and\n"
  "cortical neurons are correlated with beta subunit expression. J Neurosci\n"
  "19:7268-7277.\n"
  "\n"
  "Randall A, Tsien RW (1995) Pharmacological dissection of multiple types\n"
  "of Ca2+ channel currents in rat cerebellar granule neurons. J Neurosci\n"
  "15:2995-3012.\n"
  "\n"
  "Koch, C., and Segev, I., eds. (1998). Methods in Neuronal Modeling: From\n"
  "Ions to Networks, 2 edn (Cambridge, MA, MIT Press).\n"
  "\n"
  "Hille, B. (1992). Ionic Channels of Excitable Membranes, 2 edn\n"
  "(Sunderland, MA, Sinauer Associates Inc.).\n"
  "\n"
  "\n"
  "\n"
  "This is the w-agatoxin IVA (low conc (~20 nM) blocks P, high conc (~200 nM)\n"
  "blocks Q & P) sensitive current in fig 5 from Churchill - no P current.\n"
  "\n"
  "This current does not inactivate.  Mermelstein 99 - very slow\n"
  "inactivation (2 s). Churchill 98 - small-to-no inactivating component -\n"
  "so ignore inactivation.\n"
  "\n"
  "\n"
  "\n"
  "The standard HH model uses a linear approximation to the driving force\n"
  "for an ion: (Vm - ez).  This is ok for na and k, but not ca - calcium\n"
  "rectifies at high potentials because \n"
  "	1. internal and external concentrations of ca are so different,\n"
  "	making outward current flow much more difficult than inward \n"
  "	2. calcium is divalent so rectification is more sudden than for na\n"
  "	and k. (Hille 1992, pg 107)\n"
  "\n"
  "Accordingly, we need to replace the HH formulation with the GHK model,\n"
  "which accounts for this phenomenon.  The GHK equation is eq 6.6 in Koch\n"
  "1998, pg 217 - it expresses Ica in terms of Ca channel permeability\n"
  "(Perm,ca) times a mess. The mess can be circumvented using the ghk\n"
  "function below, which is included in the Neuron share files.  Perm,ca\n"
  "can be expressed in an HH-like fashion as \n"
  "	Perm,ca = pcabar * mca * mca 	(or however many m's and h's)\n"
  "where pcabar has dimensions of permeability but can be thought of as max\n"
  "conductance (Koch says it should be about 10^7 times smaller than the HH\n"
  "gbar - dont know) and mca is analagous to m (check out Koch 1998 pg 144)\n"
  "\n"
  "Calcium current can then be modeled as \n"
  "	ica = pcabar * mca * mca * ghk()\n"
  "\n"
  "Jason Moyer 2004 - jtmoyer@seas.upenn.edu\n"
  "ENDCOMMENT\n"
  "\n"
  ;
#endif
