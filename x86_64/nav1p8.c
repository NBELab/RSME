/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__nav1p8
#define _nrn_initial _nrn_initial__nav1p8
#define nrn_cur _nrn_cur__nav1p8
#define _nrn_current _nrn_current__nav1p8
#define nrn_jacob _nrn_jacob__nav1p8
#define nrn_state _nrn_state__nav1p8
#define _net_receive _net_receive__nav1p8 
#define states states__nav1p8 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gbar _p[0]
#define ena _p[1]
#define i _p[2]
#define m _p[3]
#define h _p[4]
#define g _p[5]
#define tau_h _p[6]
#define tau_m _p[7]
#define minf _p[8]
#define hinf _p[9]
#define Dm _p[10]
#define Dh _p[11]
#define v _p[12]
#define _g _p[13]
 
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
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_alphah(void);
 static void _hoc_alpham(void);
 static void _hoc_betah(void);
 static void _hoc_betam(void);
 static void _hoc_rates(void);
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_nav1p8", _hoc_setdata,
 "alphah_nav1p8", _hoc_alphah,
 "alpham_nav1p8", _hoc_alpham,
 "betah_nav1p8", _hoc_betah,
 "betam_nav1p8", _hoc_betam,
 "rates_nav1p8", _hoc_rates,
 0, 0
};
#define alphah alphah_nav1p8
#define alpham alpham_nav1p8
#define betah betah_nav1p8
#define betam betam_nav1p8
#define rates rates_nav1p8
 extern double alphah( _threadargsprotocomma_ double );
 extern double alpham( _threadargsprotocomma_ double );
 extern double betah( _threadargsprotocomma_ double );
 extern double betam( _threadargsprotocomma_ double );
 extern double rates( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define A_bh8 A_bh8_nav1p8
 double A_bh8 = 0.61714;
#define A_bm8 A_bm8_nav1p8
 double A_bm8 = 6.894;
#define A_ah8 A_ah8_nav1p8
 double A_ah8 = 0.013536;
#define A_am8 A_am8_nav1p8
 double A_am8 = 3.83;
#define B_bh8 B_bh8_nav1p8
 double B_bh8 = -21.8;
#define B_bm8 B_bm8_nav1p8
 double B_bm8 = 61.2;
#define B_ah8 B_ah8_nav1p8
 double B_ah8 = 105;
#define B_am8 B_am8_nav1p8
 double B_am8 = 2.58;
#define C_bh8 C_bh8_nav1p8
 double C_bh8 = -11.998;
#define C_bm8 C_bm8_nav1p8
 double C_bm8 = 19.8;
#define C_ah8 C_ah8_nav1p8
 double C_ah8 = 46.33;
#define C_am8 C_am8_nav1p8
 double C_am8 = -11.47;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "A_am8_nav1p8", "/ms",
 "B_am8_nav1p8", "mV",
 "C_am8_nav1p8", "mV",
 "A_ah8_nav1p8", "/ms",
 "B_ah8_nav1p8", "mV",
 "C_ah8_nav1p8", "mV",
 "A_bm8_nav1p8", "/ms",
 "B_bm8_nav1p8", "mV",
 "C_bm8_nav1p8", "mV",
 "A_bh8_nav1p8", "/ms",
 "B_bh8_nav1p8", "mV",
 "C_bh8_nav1p8", "mV",
 "ena_nav1p8", "mV",
 "i_nav1p8", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "A_am8_nav1p8", &A_am8_nav1p8,
 "B_am8_nav1p8", &B_am8_nav1p8,
 "C_am8_nav1p8", &C_am8_nav1p8,
 "A_ah8_nav1p8", &A_ah8_nav1p8,
 "B_ah8_nav1p8", &B_ah8_nav1p8,
 "C_ah8_nav1p8", &C_ah8_nav1p8,
 "A_bm8_nav1p8", &A_bm8_nav1p8,
 "B_bm8_nav1p8", &B_bm8_nav1p8,
 "C_bm8_nav1p8", &C_bm8_nav1p8,
 "A_bh8_nav1p8", &A_bh8_nav1p8,
 "B_bh8_nav1p8", &B_bh8_nav1p8,
 "C_bh8_nav1p8", &C_bh8_nav1p8,
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
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"nav1p8",
 "gbar_nav1p8",
 "ena_nav1p8",
 0,
 "i_nav1p8",
 0,
 "m_nav1p8",
 "h_nav1p8",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gbar = 2.2e-05;
 	ena = 79.6;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _nav1p8_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 nav1p8 /Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/x86_64/nav1p8.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / tau_m ;
   Dh = ( hinf - h ) / tau_h ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_h )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( minf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_h)))*(- ( ( ( hinf ) ) / tau_h ) / ( ( ( ( - 1.0 ) ) ) / tau_h ) - h) ;
   }
  return 0;
}
 
double alpham ( _threadargsprotocomma_ double _lVm ) {
   double _lalpham;
 _lalpham = A_am8 / ( 1.0 + exp ( ( _lVm + B_am8 ) / C_am8 ) ) ;
   
return _lalpham;
 }
 
static void _hoc_alpham(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alpham ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphah ( _threadargsprotocomma_ double _lVm ) {
   double _lalphah;
 _lalphah = A_ah8 * exp ( - ( _lVm + B_ah8 ) / C_ah8 ) ;
   
return _lalphah;
 }
 
static void _hoc_alphah(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphah ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betam ( _threadargsprotocomma_ double _lVm ) {
   double _lbetam;
 _lbetam = A_bm8 / ( 1.0 + exp ( ( _lVm + B_bm8 ) / C_bm8 ) ) ;
   
return _lbetam;
 }
 
static void _hoc_betam(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betam ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betah ( _threadargsprotocomma_ double _lVm ) {
   double _lbetah;
 _lbetah = A_bh8 / ( 1.0 + exp ( ( _lVm + B_bh8 ) / C_bh8 ) ) ;
   
return _lbetah;
 }
 
static void _hoc_betah(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betah ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double rates ( _threadargsprotocomma_ double _lVm ) {
   double _lrates;
 tau_m = 1.0 / ( alpham ( _threadargscomma_ _lVm ) + betam ( _threadargscomma_ _lVm ) ) ;
   minf = alpham ( _threadargscomma_ _lVm ) * tau_m ;
   tau_h = 1.0 / ( alphah ( _threadargscomma_ _lVm ) + betah ( _threadargscomma_ _lVm ) ) ;
   hinf = alphah ( _threadargscomma_ _lVm ) * tau_h ;
   
return _lrates;
 }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   m = alpham ( _threadargscomma_ v ) / ( alpham ( _threadargscomma_ v ) + betam ( _threadargscomma_ v ) ) ;
   h = alphah ( _threadargscomma_ v ) / ( alphah ( _threadargscomma_ v ) + betah ( _threadargscomma_ v ) ) ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gbar * pow( m , 3.0 ) * h ;
   i = g * ( v - ena ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/mods/nav1p8.mod";
static const char* nmodl_file_text = 
  ": nav1p8.mod is the NaV1.8 Na+ current from\n"
  ": Baker 2005, parameter assignments and formula's from page 854\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX nav1p8\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE gbar, ena\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(S) = (siemens)\n"
  "	(mV) = (millivolts)\n"
  "	(mA) = (milliamp)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar = 22e-6 : =220e-9/(100e-12*1e8) (S/cm2) : 220(nS)/100(um)^2\n"
  "	ena=79.6 (mV)\n"
  "\n"
  "	A_am8 = 3.83 (/ms) : A for alpha m(8 etc ...)\n"
  "	B_am8 = 2.58 (mV)\n"
  "	C_am8 = -11.47 (mV)\n"
  "\n"
  "	A_ah8 = 0.013536 (/ms) : A for alpha h\n"
  "	B_ah8 = 105 (mV)\n"
  "	C_ah8 = 46.33 (mV)\n"
  "\n"
  "	A_bm8 = 6.894 (/ms) : A for beta m\n"
  "	B_bm8 = 61.2 (mV)\n"
  "	C_bm8 = 19.8 (mV)\n"
  "\n"
  "	A_bh8 = 0.61714 (/ms)   : A for beta h\n"
  "	B_bh8 = -21.8 (mV)\n"
  "	C_bh8 = -11.998 (mV)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV) : NEURON provides this\n"
  "	i	(mA/cm2)\n"
  "	g	(S/cm2)\n"
  "	tau_h	(ms)\n"
  "	tau_m	(ms)\n"
  "	minf\n"
  "	hinf\n"
  "}\n"
  "\n"
  "STATE { m h }\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	g = gbar * m^3 * h\n"
  "	i = g * (v-ena)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	: assume that equilibrium has been reached\n"
  "	m = alpham(v)/(alpham(v)+betam(v))\n"
  "	h = alphah(v)/(alphah(v)+betah(v))\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	m' = (minf - m)/tau_m\n"
  "	h' = (hinf - h)/tau_h\n"
  "}\n"
  "\n"
  "FUNCTION alpham(Vm (mV)) (/ms) {\n"
  "	alpham=A_am8/(1+exp((Vm+B_am8)/C_am8))\n"
  "}\n"
  "\n"
  "FUNCTION alphah(Vm (mV)) (/ms) {\n"
  "	alphah=A_ah8*exp(-(Vm+B_ah8)/C_ah8)\n"
  "}\n"
  "\n"
  "FUNCTION betam(Vm (mV)) (/ms) {\n"
  "	betam=A_bm8/(1+exp((Vm+B_bm8)/C_bm8))\n"
  "}\n"
  "\n"
  "FUNCTION betah(Vm (mV)) (/ms) {\n"
  "	betah=A_bh8/(1+exp((Vm+B_bh8)/C_bh8))\n"
  "}\n"
  "\n"
  "FUNCTION rates(Vm (mV)) (/ms) {\n"
  "	tau_m = 1.0 / (alpham(Vm) + betam(Vm))\n"
  "	minf = alpham(Vm) * tau_m\n"
  "\n"
  "	tau_h = 1.0 / (alphah(Vm) + betah(Vm))\n"
  "	hinf = alphah(Vm) * tau_h\n"
  "}\n"
  ;
#endif
