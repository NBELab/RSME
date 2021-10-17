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
 
#define nrn_init _nrn_init__HHst
#define _nrn_initial _nrn_initial__HHst
#define nrn_cur _nrn_cur__HHst
#define _nrn_current _nrn_current__HHst
#define nrn_jacob _nrn_jacob__HHst
#define nrn_state _nrn_state__HHst
#define _net_receive _net_receive__HHst 
#define rates rates__HHst 
#define states states__HHst 
 
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
#define gnabar _p[0]
#define gkbar _p[1]
#define glbar _p[2]
#define gtbar _p[3]
#define gkmbar _p[4]
#define gleak _p[5]
#define eleak _p[6]
#define NFleak _p[7]
#define ileak _p[8]
#define gna _p[9]
#define gk _p[10]
#define gl _p[11]
#define gt _p[12]
#define Nna _p[13]
#define Nk _p[14]
#define Nkm _p[15]
#define Nl _p[16]
#define Nt _p[17]
#define m_exp _p[18]
#define h_exp _p[19]
#define n_exp _p[20]
#define km_exp _p[21]
#define tm_exp _p[22]
#define lm_exp _p[23]
#define th_exp _p[24]
#define lh_exp _p[25]
#define m_inf _p[26]
#define h_inf _p[27]
#define n_inf _p[28]
#define km_inf _p[29]
#define tm_inf _p[30]
#define lm_inf _p[31]
#define th_inf _p[32]
#define lh_inf _p[33]
#define noise_zn (_p + 34)
#define noise_zm (_p + 37)
#define noise_zkm _p[43]
#define noise_zt (_p + 44)
#define noise_zl (_p + 49)
#define var_zn (_p + 54)
#define var_zm (_p + 57)
#define var_zkm _p[63]
#define var_zt (_p + 64)
#define var_zl (_p + 69)
#define tau_m _p[74]
#define tau_h _p[75]
#define tau_n _p[76]
#define tau_km _p[77]
#define tau_lm _p[78]
#define tau_tm _p[79]
#define tau_th _p[80]
#define tau_lh _p[81]
#define tau_zn (_p + 82)
#define tau_zm (_p + 85)
#define tau_zkm _p[91]
#define tau_zt (_p + 92)
#define tau_zl (_p + 97)
#define mu_zn (_p + 102)
#define mu_zm (_p + 105)
#define mu_zkm _p[111]
#define mu_zt (_p + 112)
#define mu_zl (_p + 117)
#define m _p[122]
#define h _p[123]
#define n _p[124]
#define nm _p[125]
#define tm _p[126]
#define th _p[127]
#define lm _p[128]
#define lh _p[129]
#define zn (_p + 130)
#define zm (_p + 133)
#define zkm _p[139]
#define zt (_p + 140)
#define zl (_p + 145)
#define zleak _p[150]
#define Dm _p[151]
#define Dh _p[152]
#define Dn _p[153]
#define Dnm _p[154]
#define Dtm _p[155]
#define Dth _p[156]
#define Dlm _p[157]
#define Dlh _p[158]
#define Dzn (_p + 159)
#define Dzm (_p + 162)
#define Dzkm _p[168]
#define Dzt (_p + 169)
#define Dzl (_p + 174)
#define Dzleak _p[179]
#define ina _p[180]
#define il _p[181]
#define it _p[182]
#define ica _p[183]
#define ikdr _p[184]
#define ikm _p[185]
#define ik _p[186]
#define gkm _p[187]
#define ena _p[188]
#define ek _p[189]
#define eca _p[190]
#define v _p[191]
#define _g _p[192]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
#define _ion_eca	*_ppvar[6]._pval
#define _ion_ica	*_ppvar[7]._pval
#define _ion_dicadv	*_ppvar[8]._pval
#define area	*_ppvar[9]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_mulnoise(void);
 static void _hoc_numchan(void);
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_vtrap(void);
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
 "setdata_HHst", _hoc_setdata,
 "mulnoise_HHst", _hoc_mulnoise,
 "numchan_HHst", _hoc_numchan,
 "rates_HHst", _hoc_rates,
 "states_HHst", _hoc_states,
 "vtrap_HHst", _hoc_vtrap,
 0, 0
};
#define mulnoise mulnoise_HHst
#define numchan numchan_HHst
#define vtrap vtrap_HHst
 extern double mulnoise( _threadargsprotocomma_ double , double , double );
 extern double numchan( _threadargsprotocomma_ double );
 extern double vtrap( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define NF NF_HHst
 double NF = 1;
#define gamma_t gamma_t_HHst
 double gamma_t = 10;
#define gamma_l gamma_l_HHst
 double gamma_l = 10;
#define gamma_km gamma_km_HHst
 double gamma_km = 10;
#define gamma_k gamma_k_HHst
 double gamma_k = 10;
#define gamma_na gamma_na_HHst
 double gamma_na = 10;
#define hfast hfast_HHst
 double hfast = 0.3;
#define hslow hslow_HHst
 double hslow = 100;
#define seed seed_HHst
 double seed = 1;
#define taukm taukm_HHst
 double taukm = 1;
#define vshift vshift_HHst
 double vshift = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gamma_na_HHst", "pS",
 "gamma_k_HHst", "pS",
 "gamma_km_HHst", "pS",
 "gamma_t_HHst", "pS",
 "gamma_l_HHst", "pS",
 "gnabar_HHst", "S/cm2",
 "gkbar_HHst", "S/cm2",
 "glbar_HHst", "S/cm2",
 "gtbar_HHst", "S/cm2",
 "gkmbar_HHst", "S/cm2",
 "gleak_HHst", "S/cm2",
 "eleak_HHst", "mV",
 "ileak_HHst", "mA/cm2",
 "gna_HHst", "S/cm2",
 "gk_HHst", "S/cm2",
 "gl_HHst", "S/cm2",
 "gt_HHst", "S/cm2",
 "Nna_HHst", "1",
 "Nk_HHst", "1",
 "Nkm_HHst", "1",
 "Nl_HHst", "1",
 "Nt_HHst", "1",
 "tau_m_HHst", "ms",
 "tau_h_HHst", "ms",
 "tau_n_HHst", "ms",
 "tau_km_HHst", "ms",
 "tau_lm_HHst", "ms",
 "tau_tm_HHst", "ms",
 "tau_th_HHst", "ms",
 "tau_lh_HHst", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double lh0 = 0;
 static double lm0 = 0;
 static double m0 = 0;
 static double nm0 = 0;
 static double n0 = 0;
 static double th0 = 0;
 static double tm0 = 0;
 static double zleak0 = 0;
 static double zl0 = 0;
 static double zt0 = 0;
 static double zkm0 = 0;
 static double zm0 = 0;
 static double zn0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "gamma_na_HHst", &gamma_na_HHst,
 "gamma_k_HHst", &gamma_k_HHst,
 "gamma_km_HHst", &gamma_km_HHst,
 "gamma_t_HHst", &gamma_t_HHst,
 "gamma_l_HHst", &gamma_l_HHst,
 "seed_HHst", &seed_HHst,
 "vshift_HHst", &vshift_HHst,
 "taukm_HHst", &taukm_HHst,
 "NF_HHst", &NF_HHst,
 "hslow_HHst", &hslow_HHst,
 "hfast_HHst", &hfast_HHst,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"HHst",
 "gnabar_HHst",
 "gkbar_HHst",
 "glbar_HHst",
 "gtbar_HHst",
 "gkmbar_HHst",
 "gleak_HHst",
 "eleak_HHst",
 "NFleak_HHst",
 0,
 "ileak_HHst",
 "gna_HHst",
 "gk_HHst",
 "gl_HHst",
 "gt_HHst",
 "Nna_HHst",
 "Nk_HHst",
 "Nkm_HHst",
 "Nl_HHst",
 "Nt_HHst",
 "m_exp_HHst",
 "h_exp_HHst",
 "n_exp_HHst",
 "km_exp_HHst",
 "tm_exp_HHst",
 "lm_exp_HHst",
 "th_exp_HHst",
 "lh_exp_HHst",
 "m_inf_HHst",
 "h_inf_HHst",
 "n_inf_HHst",
 "km_inf_HHst",
 "tm_inf_HHst",
 "lm_inf_HHst",
 "th_inf_HHst",
 "lh_inf_HHst",
 "noise_zn_HHst[3]",
 "noise_zm_HHst[6]",
 "noise_zkm_HHst",
 "noise_zt_HHst[5]",
 "noise_zl_HHst[5]",
 "var_zn_HHst[3]",
 "var_zm_HHst[6]",
 "var_zkm_HHst",
 "var_zt_HHst[5]",
 "var_zl_HHst[5]",
 "tau_m_HHst",
 "tau_h_HHst",
 "tau_n_HHst",
 "tau_km_HHst",
 "tau_lm_HHst",
 "tau_tm_HHst",
 "tau_th_HHst",
 "tau_lh_HHst",
 "tau_zn_HHst[3]",
 "tau_zm_HHst[6]",
 "tau_zkm_HHst",
 "tau_zt_HHst[5]",
 "tau_zl_HHst[5]",
 "mu_zn_HHst[3]",
 "mu_zm_HHst[6]",
 "mu_zkm_HHst",
 "mu_zt_HHst[5]",
 "mu_zl_HHst[5]",
 0,
 "m_HHst",
 "h_HHst",
 "n_HHst",
 "nm_HHst",
 "tm_HHst",
 "th_HHst",
 "lm_HHst",
 "lh_HHst",
 "zn_HHst[3]",
 "zm_HHst[6]",
 "zkm_HHst",
 "zt_HHst[5]",
 "zl_HHst[5]",
 "zleak_HHst",
 0,
 0};
 extern Node* nrn_alloc_node_;
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 193, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.12;
 	gkbar = 0.036;
 	glbar = 0.0003;
 	gtbar = 0.0003;
 	gkmbar = 0.002;
 	gleak = 1e-05;
 	eleak = -60;
 	NFleak = 1;
 	_prop->param = _p;
 	_prop->param_size = 193;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 10, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 	_ppvar[9]._pval = &nrn_alloc_node_->_area; /* diam */
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[6]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[7]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[8]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _HHst_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 193, 10);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "area");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 HHst /Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/x86_64/HHst.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Stochastic Hodgkin and Huxley model & M-type potassium & T-and L-type Calcium channels incorporating channel noise .";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
static int states(_threadargsproto_);
 
static int  states ( _threadargsproto_ ) {
   rates ( _threadargscomma_ v + vshift ) ;
   m = m + m_exp * ( m_inf - m ) ;
   h = h + h_exp * ( h_inf - h ) ;
   n = n + n_exp * ( n_inf - n ) ;
   nm = nm + km_exp * ( km_inf - nm ) ;
   tm = tm + tm_exp * ( tm_inf - tm ) ;
   th = th + th_exp * ( th_inf - th ) ;
   lm = lm + lm_exp * ( lm_inf - lm ) ;
   lh = lh + lh_exp * ( lh_inf - lh ) ;
   
/*VERBATIM*/
    return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 states ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lvm ) {
   double _la , _lb , _lm3_inf , _ln4_inf , _lsum , _lone_minus_m , _lone_minus_h , _lone_minus_n , _li ;
  _la = - .6 * vtrap ( _threadargscomma_ ( _lvm + 30.0 ) , - 10.0 ) ;
   _lb = 20.0 * ( exp ( ( - 1.0 * ( _lvm + 55.0 ) ) / 18.0 ) ) ;
   tau_m = 1.0 / ( _la + _lb ) ;
   m_inf = _la * tau_m ;
   _lone_minus_m = 1. - m_inf ;
   _lm3_inf = m_inf * m_inf * m_inf ;
   _la = 0.4 * ( exp ( ( - 1.0 * ( _lvm + 50.0 ) ) / 20.0 ) ) ;
   _lb = 6.0 / ( 1.0 + exp ( - 0.1 * ( _lvm + 20.0 ) ) ) ;
   tau_h = hslow / ( ( 1.0 + exp ( ( _lvm + 30.0 ) / 4.0 ) ) + ( exp ( - ( _lvm + 50.0 ) / 2.0 ) ) ) + hfast ;
   h_inf = 1.0 / ( 1.0 + exp ( ( _lvm + 44.0 ) / 4.0 ) ) ;
   _lone_minus_h = 1. - h_inf ;
   tau_zm [ 0 ] = tau_h ;
   tau_zm [ 1 ] = tau_m ;
   tau_zm [ 2 ] = tau_m / 2.0 ;
   tau_zm [ 3 ] = tau_m / 3.0 ;
   tau_zm [ 4 ] = tau_m * tau_h / ( tau_m + tau_h ) ;
   tau_zm [ 5 ] = tau_m * tau_h / ( tau_m + 2.0 * tau_h ) ;
   tau_zm [ 6 ] = tau_m * tau_h / ( tau_m + 3.0 * tau_h ) ;
   var_zm [ 0 ] = 1.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * _lm3_inf * h_inf * _lone_minus_h ;
   var_zm [ 1 ] = 3.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * m_inf * m_inf * h_inf * h_inf * _lone_minus_m ;
   var_zm [ 2 ] = 3.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * m_inf * h_inf * h_inf * _lone_minus_m * _lone_minus_m ;
   var_zm [ 3 ] = 1.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * h_inf * h_inf * _lone_minus_m * _lone_minus_m * _lone_minus_m ;
   var_zm [ 4 ] = 3.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * m_inf * m_inf * h_inf * _lone_minus_m * _lone_minus_h ;
   var_zm [ 5 ] = 3.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * m_inf * h_inf * _lone_minus_m * _lone_minus_m * _lone_minus_h ;
   var_zm [ 6 ] = 1.0 / numchan ( _threadargscomma_ Nna ) * _lm3_inf * h_inf * _lone_minus_m * _lone_minus_m * _lone_minus_m * _lone_minus_h ;
   {int  _li ;for ( _li = 0 ; _li <= 6 ; _li ++ ) {
     mu_zm [ _li ] = exp ( - dt / tau_zm [ _li ] ) ;
     noise_zm [ _li ] = sqrt ( var_zm [ _li ] * ( 1.0 - mu_zm [ _li ] * mu_zm [ _li ] ) ) * normrand ( 0.0 , NF ) ;
     zm [ _li ] = zm [ _li ] * mu_zm [ _li ] + noise_zm [ _li ] ;
     } }
   _la = - 0.02 * vtrap ( _threadargscomma_ ( _lvm + 40.0 ) , - 10.0 ) ;
   _lb = 0.4 * ( exp ( ( - 1.0 * ( _lvm + 50.0 ) ) / 80.0 ) ) ;
   tau_n = 1.0 / ( _la + _lb ) ;
   n_inf = _la * tau_n ;
   _lone_minus_n = 1. - n_inf ;
   _ln4_inf = n_inf * n_inf * n_inf * n_inf ;
   tau_zn [ 0 ] = tau_n ;
   tau_zn [ 1 ] = tau_n / 2.0 ;
   tau_zn [ 2 ] = tau_n / 3.0 ;
   tau_zn [ 3 ] = tau_n / 4.0 ;
   var_zn [ 0 ] = 4.0 / numchan ( _threadargscomma_ Nk ) * _ln4_inf * n_inf * n_inf * n_inf * _lone_minus_n ;
   var_zn [ 1 ] = 6.0 / numchan ( _threadargscomma_ Nk ) * _ln4_inf * n_inf * n_inf * _lone_minus_n * _lone_minus_n ;
   var_zn [ 2 ] = 4.0 / numchan ( _threadargscomma_ Nk ) * _ln4_inf * n_inf * _lone_minus_n * _lone_minus_n * _lone_minus_n ;
   var_zn [ 3 ] = 1.0 / numchan ( _threadargscomma_ Nk ) * _ln4_inf * _lone_minus_n * _lone_minus_n * _lone_minus_n * _lone_minus_n ;
   {int  _li ;for ( _li = 0 ; _li <= 3 ; _li ++ ) {
     mu_zn [ _li ] = exp ( - dt / tau_zn [ _li ] ) ;
     noise_zn [ _li ] = sqrt ( var_zn [ _li ] * ( 1.0 - mu_zn [ _li ] * mu_zn [ _li ] ) ) * normrand ( 0.0 , NF ) ;
     zn [ _li ] = zn [ _li ] * mu_zn [ _li ] + noise_zn [ _li ] ;
     } }
   _la = - .001 / taukm * vtrap ( _threadargscomma_ ( _lvm + 30.0 ) , - 9.0 ) ;
   _lb = .001 / taukm * vtrap ( _threadargscomma_ ( _lvm + 30.0 ) , 9.0 ) ;
   tau_km = 1.0 / ( _la + _lb ) ;
   km_inf = _la * tau_km ;
   tau_zkm = tau_km ;
   var_zkm = km_inf * ( 1.0 - km_inf ) / numchan ( _threadargscomma_ Nkm ) ;
   mu_zkm = exp ( - dt / tau_zkm ) ;
   noise_zkm = sqrt ( var_zkm * ( 1.0 - mu_zkm * mu_zkm ) ) * normrand ( 0.0 , NF ) ;
   zkm = zkm * mu_zkm + noise_zkm ;
   _la = 0.055 * vtrap ( _threadargscomma_ - ( _lvm + 27.0 ) , 3.8 ) ;
   _lb = 0.94 * exp ( ( - 75.0 - _lvm ) / 17.0 ) ;
   tau_lm = 1.0 / ( _la + _lb ) ;
   lm_inf = _la * tau_lm ;
   _la = 0.000457 * exp ( ( - 13.0 - _lvm ) / 50.0 ) ;
   _lb = 0.0065 / ( exp ( ( - _lvm - 15.0 ) / 28.0 ) + 1.0 ) ;
   tau_lh = 1.0 / ( _la + _lb ) ;
   lh_inf = _la * tau_lh ;
   tau_zl [ 0 ] = tau_lm * tau_lh / ( tau_lm + 2.0 * tau_lh ) ;
   tau_zl [ 1 ] = tau_lm * tau_lh / ( tau_lm + tau_lh ) ;
   tau_zl [ 2 ] = tau_lh ;
   tau_zl [ 3 ] = tau_lm / 2.0 ;
   tau_zl [ 4 ] = tau_lm ;
   var_zl [ 0 ] = 1.0 / numchan ( _threadargscomma_ Nl ) * pow( lm_inf , 2.0 ) * lh_inf * pow( ( 1.0 - lm_inf ) , 2.0 ) * ( 1.0 - lh_inf ) ;
   var_zl [ 1 ] = 2.0 / numchan ( _threadargscomma_ Nl ) * pow( lm_inf , 3.0 ) * lh_inf * ( 1.0 - lm_inf ) * ( 1.0 - lh_inf ) ;
   var_zl [ 2 ] = 1.0 / numchan ( _threadargscomma_ Nl ) * pow( lm_inf , 4.0 ) * lh_inf * ( 1.0 - lh_inf ) ;
   var_zl [ 3 ] = 1.0 / numchan ( _threadargscomma_ Nl ) * pow( lm_inf , 2.0 ) * pow( lh_inf , 2.0 ) * pow( ( 1.0 - lm_inf ) , 2.0 ) ;
   var_zl [ 4 ] = 2.0 / numchan ( _threadargscomma_ Nl ) * pow( lm_inf , 3.0 ) * pow( lh_inf , 2.0 ) * ( 1.0 - lm_inf ) ;
   {int  _li ;for ( _li = 0 ; _li <= 4 ; _li ++ ) {
     mu_zl [ _li ] = exp ( - dt / tau_zl [ _li ] ) ;
     noise_zl [ _li ] = sqrt ( var_zl [ _li ] * ( 1.0 - mu_zl [ _li ] * mu_zl [ _li ] ) ) * normrand ( 0.0 , NF ) ;
     zl [ _li ] = zl [ _li ] * mu_zl [ _li ] + noise_zl [ _li ] ;
     } }
   tau_tm = ( 4.0 / ( exp ( ( _lvm + 25.0 ) / 20.0 ) + exp ( - ( _lvm + 100.0 ) / 15.0 ) ) ) ;
   tm_inf = 1.0 / ( 1.0 + exp ( - ( _lvm + 50.0 ) / 7.4 ) ) ;
   tau_th = ( 86.0 / ( exp ( ( _lvm + 46.0 ) / 4.0 ) + exp ( - ( _lvm + 405.0 ) / 50.0 ) ) ) ;
   th_inf = 1.0 / ( 1.0 + exp ( ( _lvm + 78.0 ) / 5.0 ) ) ;
   tau_zt [ 0 ] = tau_tm * tau_th / ( tau_tm + 2.0 * tau_th ) ;
   tau_zt [ 1 ] = tau_tm * tau_th / ( tau_tm + tau_th ) ;
   tau_zt [ 2 ] = tau_th ;
   tau_zt [ 3 ] = tau_tm / 2.0 ;
   tau_zt [ 4 ] = tau_tm ;
   var_zt [ 0 ] = 1.0 / numchan ( _threadargscomma_ Nt ) * pow( tm_inf , 2.0 ) * th_inf * pow( ( 1.0 - tm_inf ) , 2.0 ) * ( 1.0 - th_inf ) ;
   var_zt [ 1 ] = 2.0 / numchan ( _threadargscomma_ Nt ) * pow( tm_inf , 3.0 ) * th_inf * ( 1.0 - tm_inf ) * ( 1.0 - th_inf ) ;
   var_zt [ 2 ] = 1.0 / numchan ( _threadargscomma_ Nt ) * pow( tm_inf , 4.0 ) * th_inf * ( 1.0 - th_inf ) ;
   var_zt [ 3 ] = 1.0 / numchan ( _threadargscomma_ Nt ) * pow( tm_inf , 2.0 ) * pow( th_inf , 2.0 ) * pow( ( 1.0 - tm_inf ) , 2.0 ) ;
   var_zt [ 4 ] = 2.0 / numchan ( _threadargscomma_ Nt ) * pow( tm_inf , 3.0 ) * pow( th_inf , 2.0 ) * ( 1.0 - tm_inf ) ;
   {int  _li ;for ( _li = 0 ; _li <= 4 ; _li ++ ) {
     mu_zt [ _li ] = exp ( - dt / tau_zt [ _li ] ) ;
     noise_zt [ _li ] = sqrt ( var_zt [ _li ] * ( 1.0 - mu_zt [ _li ] * mu_zt [ _li ] ) ) * normrand ( 0.0 , NF ) ;
     zt [ _li ] = zt [ _li ] * mu_zt [ _li ] + noise_zt [ _li ] ;
     } }
   zleak = normrand ( 0.0 , NFleak ) ;
   m_exp = 1.0 - exp ( - dt / tau_m ) ;
   h_exp = 1.0 - exp ( - dt / tau_h ) ;
   n_exp = 1.0 - exp ( - dt / tau_n ) ;
   km_exp = 1.0 - exp ( - dt / tau_km ) ;
   lm_exp = 1.0 - exp ( - dt / tau_lm ) ;
   lh_exp = 1.0 - exp ( - dt / tau_lh ) ;
   tm_exp = 1.0 - exp ( - dt / tau_tm ) ;
   th_exp = 1.0 - exp ( - dt / tau_th ) ;
     return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( exp ( _lx / _ly ) - 1.0 ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double mulnoise ( _threadargsprotocomma_ double _lmean , double _lsd , double _lpower ) {
   double _lmulnoise;
 double _li , _lavg ;
 _lavg = 1.0 ;
   {int  _li ;for ( _li = 1 ; _li <= ((int) _lpower ) ; _li ++ ) {
     _lavg = _lavg * normrand ( _lmean , _lsd ) ;
     } }
   _lmulnoise = _lavg ;
   
return _lmulnoise;
 }
 
static void _hoc_mulnoise(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  mulnoise ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double numchan ( _threadargsprotocomma_ double _lNchannels ) {
   double _lnumchan;
 if ( _lNchannels > 0.0 ) {
     _lnumchan = ( _lNchannels ) ;
     }
   else {
     _lnumchan = 1.0 ;
     }
   
return _lnumchan;
 }
 
static void _hoc_numchan(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  numchan ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("HHst", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 6, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 7, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 8, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  lh = lh0;
  lm = lm0;
  m = m0;
  n = n0;
  nm = nm0;
  th = th0;
  tm = tm0;
  zleak = zleak0;
 for (_i=0; _i<5; _i++) zl[_i] = zl0;
 for (_i=0; _i<5; _i++) zt[_i] = zt0;
  zkm = zkm0;
 for (_i=0; _i<6; _i++) zm[_i] = zm0;
 for (_i=0; _i<3; _i++) zn[_i] = zn0;
 {
   Nna = ceil ( ( ( 1e-8 ) * area ) * ( gnabar ) / ( ( 1e-12 ) * gamma_na ) ) ;
   Nk = ceil ( ( ( 1e-8 ) * area ) * ( gkbar ) / ( ( 1e-12 ) * gamma_k ) ) ;
   Nkm = ceil ( ( ( 1e-8 ) * area ) * ( gkmbar ) / ( ( 1e-12 ) * gamma_km ) ) ;
   Nt = ceil ( ( ( 1e-8 ) * area ) * ( gtbar ) / ( ( 1e-12 ) * gamma_t ) ) ;
   Nl = ceil ( ( ( 1e-8 ) * area ) * ( glbar ) / ( ( 1e-12 ) * gamma_l ) ) ;
   rates ( _threadargscomma_ v ) ;
   m = m_inf ;
   h = h_inf ;
   n = n_inf ;
   nm = km_inf ;
   tm = tm_inf ;
   th = th_inf ;
   lm = lm_inf ;
   lh = lh_inf ;
   zn [ 0 ] = 0.0 ;
   zn [ 1 ] = 0.0 ;
   zn [ 2 ] = 0.0 ;
   zn [ 3 ] = 0.0 ;
   zm [ 0 ] = 0. ;
   zm [ 1 ] = 0. ;
   zm [ 2 ] = 0. ;
   zm [ 3 ] = 0. ;
   zm [ 4 ] = 0. ;
   zm [ 5 ] = 0. ;
   zm [ 6 ] = 0. ;
   zkm = 0.0 ;
   zt [ 0 ] = 0.0 ;
   zt [ 1 ] = 0.0 ;
   zt [ 2 ] = 0.0 ;
   zt [ 3 ] = 0.0 ;
   zt [ 4 ] = 0.0 ;
   zl [ 0 ] = 0.0 ;
   zl [ 1 ] = 0.0 ;
   zl [ 2 ] = 0.0 ;
   zl [ 3 ] = 0.0 ;
   zl [ 4 ] = 0.0 ;
   zleak = 0.0 ;
   set_seed ( seed ) ;
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
  ena = _ion_ena;
  ek = _ion_ek;
  eca = _ion_eca;
 initmodel(_p, _ppvar, _thread, _nt);
   }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnabar * ( m * m * m * h + ( zm [ 0 ] + zm [ 1 ] + zm [ 2 ] + zm [ 3 ] + zm [ 4 ] + zm [ 5 ] + zm [ 6 ] ) ) ;
   ina = gna * ( v - ena ) ;
   gk = gkbar * ( n * n * n * n + ( zn [ 0 ] + zn [ 1 ] + zn [ 2 ] + zn [ 3 ] ) ) ;
   ikdr = gk * ( v - ek ) ;
   gkm = gkmbar * ( nm + zkm ) ;
   ikm = gkm * ( v - ek ) ;
   ik = ikm + ikdr ;
   gt = gtbar * ( tm * tm * th + ( zt [ 0 ] + zt [ 1 ] + zt [ 2 ] + zt [ 3 ] + zt [ 4 ] ) ) ;
   it = gt * ( v - eca ) ;
   gl = glbar * ( lm * lm * lh + ( zl [ 0 ] + zl [ 1 ] + zl [ 2 ] + zl [ 3 ] + zl [ 4 ] ) ) ;
   il = gl * ( v - eca ) ;
   ica = il + it ;
   ileak = gleak * ( 1.0 + zleak ) * ( v - eleak ) ;
   }
 _current += ina;
 _current += ik;
 _current += ica;
 _current += ileak;

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
  ena = _ion_ena;
  ek = _ion_ek;
  eca = _ion_eca;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
 double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
  _ion_ica += ica ;
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
  ena = _ion_ena;
  ek = _ion_ek;
  eca = _ion_eca;
 {  { states(_p, _ppvar, _thread, _nt); }
  }   }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/elishai/Dropbox/1AFiles/NBEL G2/RSME/RSME_dis/mods/HHst.mod";
static const char* nmodl_file_text = 
  "\n"
  "TITLE Stochastic Hodgkin and Huxley model & M-type potassium & T-and L-type Calcium channels incorporating channel noise .\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Based on - Accurate and fast simulation of channel noise in conductance-based model neurons. Linaro, D., Storace, M., and Giugliano, M. \n"
  "Added: \n"
  "	Km T L channels\n"
  "	fixed minor bugs and grouped variables into arrays \n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "    (S)  = (siemens)\n"
  "    (pS) = (picosiemens)\n"
  "    (um) = (micrometer)\n"
  "} : end UNITS\n"
  "\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX HHst\n"
  "    USEION na READ ena WRITE ina\n"
  "    USEION k READ ek WRITE ik\n"
  "	USEION ca READ eca WRITE ica\n"
  "	NONSPECIFIC_CURRENT ileak\n"
  "	RANGE eleak, gleak,NFleak\n"
  "    RANGE gnabar, gkbar\n"
  "    RANGE gna, gk\n"
  "	RANGE lm, lh, gl, glbar\n"
  "	RANGE tm, th, gt, gtbar\n"
  "	RANGE nm,  gkmbar\n"
  "	RANGE m_exp, h_exp, n_exp, km_exp,lm_exp,tm_exp,lh_exp,th_exp\n"
  "	\n"
  "    GLOBAL gamma_na, gamma_k,gamma_km,gamma_l,gamma_t\n"
  "	RANGE km_inf, tau_km,lm_inf,tm_inf,tau_lm,tau_tm,lh_inf,th_inf,tau_lh,tau_th	\n"
  "    RANGE m_inf, h_inf, n_inf\n"
  "    RANGE tau_m, tau_h, tau_n\n"
  "    RANGE tau_zn,tau_zm,tau_zkm,tau_zt,tau_zl\n"
  "    RANGE var_zn,var_zm,var_zkm,var_zt,var_zl\n"
  "    RANGE noise_zn,noise_zm,noise_zkm,noise_zl,noise_zt\n"
  "    RANGE Nna, Nk,Nkm,Nt,Nl\n"
  "    GLOBAL seed    ,hslow,hfast\n"
  "    RANGE mu_zn,mu_zm,mu_zkm,mu_zt,mu_zl\n"
  "	:RANGE zl[4]\n"
  "	GLOBAL vshift\n"
  "\n"
  "	GLOBAL taukm,NF\n"
  "    THREADSAFE\n"
  "} : end NEURON\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "    gnabar  = 0.12   (S/cm2)\n"
  "    gkbar   = 0.036  (S/cm2)\n"
  "    glbar      = 0.0003 (S/cm2)  \n"
  "    gtbar      = 0.0003 (S/cm2)  \n"
  "	gkmbar = .002   	(S/cm2)\n"
  "    gleak=0.00001		(S/cm2)\n"
  "	eleak=-60 			(mV)\n"
  "	NFleak=1\n"
  "    gamma_na = 10  (pS)		: single channel sodium conductance\n"
  "    gamma_k  = 10  (pS)		: single channel potassium conductance\n"
  "    gamma_km  = 10  (pS)		: single channel potassium conductance\n"
  "    gamma_t  = 10  (pS)		: single channel calcium conductance\n"
  "    gamma_l  = 10  (pS)		: single channel calcium conductance\n"
  "    seed = 1              : always use the same seed\n"
  "	vshift=0				:Voltage shift of the recorded memebrane potential (to offset for liquid junction potential\n"
  "	taukm=1					:speedup of Km channels\n"
  "	NF=1					:Noise Factor (set to 0 to zero the noise part)\n"
  "	hslow=100\n"
  "	hfast=0.3\n"
  "} : end PARAMETER\n"
  "\n"
  "\n"
  "STATE {\n"
  "    m h n nm tm th lm lh\n"
  "    zn[3]\n"
  "    zm[6]\n"
  "	zkm\n"
  "	zt[5]\n"
  "	zl[5]\n"
  "	zleak\n"
  "} : end STATE\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "    ina        (mA/cm2)\n"
  "    il         (mA/cm2)\n"
  "    it         (mA/cm2)\n"
  "    ica         (mA/cm2)\n"
  "	ikdr	(mA/cm2)\n"
  "	ikm	 (mA/cm2)\n"
  "	ik		(mA/cm2)\n"
  "	ileak	(mA/cm2)\n"
  "    gna        (S/cm2)\n"
  "    gk         (S/cm2)\n"
  "	gkm         (S/cm2)\n"
  "	gl         (S/cm2)\n"
  "	gt         (S/cm2)\n"
  "    ena        (mV)\n"
  "    ek         (mV)\n"
  "    eca         (mV)\n"
  "   \n"
  "    dt         (ms)\n"
  "    area       (um2)\n"
  "    celsius    (degC)\n"
  "    v          (mV)\n"
  "        \n"
  "    Nna  (1) : number of Na channels\n"
  "    Nk   (1) : number of K channels\n"
  "    Nkm   (1) : number of Km channels\n"
  "    Nl   (1) : number of L channels\n"
  "    Nt   (1) : number of T channels\n"
  "    m_exp h_exp n_exp km_exp tm_exp lm_exp th_exp lh_exp\n"
  "    m_inf h_inf n_inf km_inf tm_inf lm_inf th_inf lh_inf\n"
  "    noise_zn[3]\n"
  "    noise_zm[6]\n"
  "	noise_zkm\n"
  "	noise_zt[5]\n"
  "	noise_zl[5]	\n"
  "    var_zn[3]\n"
  "    var_zm[6]\n"
  "	var_zkm\n"
  "	var_zt[5]\n"
  "	var_zl[5]\n"
  "    tau_m (ms) tau_h (ms) tau_n (ms) tau_km (ms) tau_lm (ms) tau_tm (ms) tau_th (ms) tau_lh (ms)\n"
  "    tau_zn[3]\n"
  "    tau_zm[6]\n"
  "	tau_zkm\n"
  "	tau_zt[5]\n"
  "	tau_zl[5]\n"
  "    mu_zn[3]\n"
  "    mu_zm[6]\n"
  "	mu_zkm\n"
  "	mu_zt[5]\n"
  "	mu_zl[5]\n"
  "} : end ASSIGNED\n"
  "\n"
  "INITIAL {\n"
  "    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   : area in um2 -> 1e-8*area in cm2; gnabar in S/cm2; gamma_na in pS -> 1e-12*gamma_na in S\n"
  "    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S\n"
  "    Nkm = ceil(((1e-8)*area)*(gkmbar)/((1e-12)*gamma_km))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S\n"
  "    Nt = ceil(((1e-8)*area)*(gtbar)/((1e-12)*gamma_t))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S\n"
  "    Nl = ceil(((1e-8)*area)*(glbar)/((1e-12)*gamma_l))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S\n"
  "    \n"
  "    rates(v)\n"
  "    m = m_inf\n"
  "    h = h_inf\n"
  "    n = n_inf\n"
  "	nm=km_inf\n"
  "	tm=tm_inf\n"
  "	th=th_inf\n"
  "	lm=lm_inf\n"
  "	lh=lh_inf\n"
  "	\n"
  "	zn[0] = 0\n"
  "	zn[1] = 0\n"
  "	zn[2] = 0\n"
  "	zn[3] = 0\n"
  "\n"
  "    zm[0] = 0.\n"
  "    zm[1] = 0.\n"
  "    zm[2] = 0.\n"
  "    zm[3] = 0.\n"
  "    zm[4] = 0.\n"
  "    zm[5] = 0.\n"
  "    zm[6] = 0.\n"
  "	\n"
  "	zkm=0\n"
  "	\n"
  "	zt[0]=0\n"
  "	zt[1]=0\n"
  "	zt[2]=0\n"
  "	zt[3]=0\n"
  "	zt[4]=0\n"
  "\n"
  "	zl[0]=0\n"
  "	zl[1]=0\n"
  "	zl[2]=0\n"
  "	zl[3]=0\n"
  "	zl[4]=0\n"
  "	\n"
  "	zleak=0\n"
  "    set_seed(seed)\n"
  "} : end INITIAL\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE states\n"
  "    gna = gnabar * (m*m*m*h + (zm[0]+zm[1]+zm[2]+zm[3]+zm[4]+zm[5]+zm[6]))\n"
  ":    if (gna < 0) {\n"
  "		:gna = 0\n"
  ":    }\n"
  "    ina = gna * (v - ena)\n"
  "    gk = gkbar * (n*n*n*n + (zn[0]+zn[1]+zn[2]+zn[3]))\n"
  ":   if (gk < 0) {\n"
  ":		gk = 0\n"
  ":    }\n"
  "    ikdr  = gk * (v - ek)\n"
  "	gkm=gkmbar*(nm+zkm)\n"
  ":	if (gkm < 0) {\n"
  ":		gkm= 0\n"
  ":    }\n"
  "	ikm = gkm*(v-ek)\n"
  "	ik=ikm+ikdr\n"
  "	\n"
  "    gt = gtbar * (tm*tm*th + (zt[0]+zt[1]+zt[2]+zt[3]+zt[4]))\n"
  ":    if (gt < 0) {\n"
  ":		gt = 0\n"
  ":    }\n"
  "    it  = gt * (v - eca)\n"
  "    gl = glbar * (lm*lm*lh + (zl[0]+zl[1]+zl[2]+zl[3]+zl[4]))\n"
  ":    if (gl < 0) {\n"
  ":		gl = 0\n"
  ":    }	\n"
  "    il  = gl * (v - eca)\n"
  "	ica=il+it\n"
  "	ileak  = gleak*(1+zleak) * (v - eleak)\n"
  ":	ileak  = (gleak) * (v - eleak)\n"
  "} : end BREAKPOINT\n"
  "\n"
  "\n"
  "PROCEDURE states() {\n"
  "    rates(v+vshift)\n"
  "	m = m + m_exp * (m_inf - m)\n"
  "	h = h + h_exp * (h_inf - h)\n"
  "    n = n + n_exp * (n_inf - n)\n"
  "	nm = nm + km_exp*(km_inf-nm)\n"
  "    tm = tm + tm_exp * (tm_inf - tm)\n"
  "    th = th + th_exp * (th_inf - th)\n"
  "    lm = lm + lm_exp * (lm_inf - lm)\n"
  "    lh = lh + lh_exp * (lh_inf - lh)\n"
  "\n"
  "    VERBATIM\n"
  "    return 0;\n"
  "    ENDVERBATIM\n"
  "} : end PROCEDURE states()\n"
  "\n"
  "\n"
  "PROCEDURE rates(vm (mV)) { \n"
  "    LOCAL a,b,m3_inf,n4_inf,sum,one_minus_m,one_minus_h,one_minus_n,i\n"
  "    \n"
  "    UNITSOFF\n"
  "    \n"
  "    \n"
  " :NA m\n"
  "	a =-.6 * vtrap((vm+30),-10)	\n"
  "	b = 20 * (exp((-1*(vm+55))/18))\n"
  "	tau_m = 1 / (a + b)\n"
  "	m_inf = a * tau_m\n"
  "\n"
  "    one_minus_m = 1. - m_inf\n"
  "    m3_inf = m_inf*m_inf*m_inf\n"
  "\n"
  ":NA h\n"
  "	a = 0.4 * (exp((-1*(vm+50))/20))\n"
  "	b = 6 / ( 1 + exp(-0.1 *(vm+20)))\n"
  "	:b = 6 / ( 1 + exp(-(vm-50)/5))\n"
  "	:tau_h = 1 / (a + b)\n"
  "	:h_inf = a * tau_h\n"
  "	:tau_h=1/(a+6 / ( 1 + exp(-(vm-50)/5)))\n"
  "	:tau_h=hslow/(1+exp((vm+30)/4))+hfast\n"
  "	tau_h=hslow/((1+exp((vm+30)/4))+(exp(-(vm+50)/2)))+hfast\n"
  "	h_inf= 1/(1 + exp((vm + 44)/4))\n"
  "	:h_inf = a / (a+b)\n"
  "\n"
  "	one_minus_h = 1. - h_inf\n"
  "    \n"
  "    tau_zm[0] = tau_h\n"
  "    tau_zm[1] = tau_m\n"
  "    tau_zm[2] = tau_m/2\n"
  "    tau_zm[3] = tau_m/3\n"
  "    tau_zm[4] = tau_m*tau_h/(tau_m+tau_h)\n"
  "    tau_zm[5] = tau_m*tau_h/(tau_m+2*tau_h)\n"
  "    tau_zm[6] = tau_m*tau_h/(tau_m+3*tau_h)\n"
  "    var_zm[0] = 1.0 / numchan(Nna) * m3_inf*m3_inf*h_inf * one_minus_h\n"
  "    var_zm[1] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m\n"
  "    var_zm[2] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m\n"
  "    var_zm[3] = 1.0 / numchan(Nna) * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m\n"
  "    var_zm[4] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h\n"
  "    var_zm[5] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h\n"
  "    var_zm[6] = 1.0 / numchan(Nna) * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h\n"
  "	\n"
  "	FROM i=0 TO 6 {\n"
  "		mu_zm[i] = exp(-dt/tau_zm[i])\n"
  "        noise_zm[i] = sqrt(var_zm[i] * (1-mu_zm[i]*mu_zm[i])) * normrand(0,NF)\n"
  "		zm[i] = zm[i]*mu_zm[i] + noise_zm[i]\n"
  "    }\n"
  "  \n"
  "   \n"
  ":K n (non-inactivating, delayed rectifier)\n"
  "	a = -0.02 * vtrap((vm+40),-10)\n"
  "	b = 0.4 * (exp((-1*(vm + 50))/80))\n"
  "	tau_n = 1 / (a + b)\n"
  "	n_inf = a * tau_n\n"
  "    one_minus_n = 1. - n_inf\n"
  "    n4_inf = n_inf * n_inf * n_inf * n_inf\n"
  "    tau_zn[0] = tau_n\n"
  "    tau_zn[1] = tau_n/2\n"
  "    tau_zn[2] = tau_n/3\n"
  "    tau_zn[3] = tau_n/4\n"
  "    var_zn[0] = 4.0/numchan(Nk) * n4_inf*n_inf*n_inf*n_inf * one_minus_n\n"
  "    var_zn[1] = 6.0/numchan(Nk) * n4_inf*n_inf*n_inf * one_minus_n*one_minus_n\n"
  "    var_zn[2] = 4.0/numchan(Nk) * n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n\n"
  "    var_zn[3] = 1.0/numchan(Nk) * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n\n"
  "\n"
  "	FROM i=0 TO 3 {\n"
  "		mu_zn[i] = exp(-dt/tau_zn[i])\n"
  "        noise_zn[i] = sqrt(var_zn[i] * (1-mu_zn[i]*mu_zn[i])) * normrand(0,NF)\n"
  "		zn[i] = zn[i]*mu_zn[i] + noise_zn[i]\n"
  "    }\n"
  ":Km\n"
  "    a = -.001/taukm * vtrap((vm+30),-9)\n"
  "    b =.001/taukm * vtrap((vm+30),9) \n"
  "    tau_km = 1/(a+b)\n"
  "	km_inf = a*tau_km\n"
  "	tau_zkm = tau_km\n"
  "	var_zkm=km_inf*(1-km_inf)/numchan(Nkm)\n"
  "	mu_zkm = exp(-dt/tau_zkm)\n"
  "    noise_zkm = sqrt(var_zkm * (1-mu_zkm*mu_zkm)) * normrand(0,NF)\n"
  "	zkm = zkm*mu_zkm + noise_zkm\n"
  "	\n"
  ":L Ca (m)\n"
  "	a = 0.055*vtrap(-(vm+27),3.8): (-27 - vm)/(exp((-27-vm)/3.8) - 1)\n"
  "	b = 0.94*exp((-75-vm)/17)\n"
  "    tau_lm = 1/(a+b)\n"
  "	lm_inf = a*tau_lm\n"
  ":L Ca (h)	\n"
  "	a = 0.000457*exp((-13-vm)/50)\n"
  "	b = 0.0065/(exp((-vm-15)/28) + 1)	\n"
  "    tau_lh = 1/(a+b)\n"
  "	lh_inf = a*tau_lh\n"
  "\n"
  "    tau_zl[0] = tau_lm*tau_lh/(tau_lm+2*tau_lh)\n"
  "    tau_zl[1] = tau_lm*tau_lh/(tau_lm+tau_lh)\n"
  "    tau_zl[2] = tau_lh\n"
  "    tau_zl[3] = tau_lm/2\n"
  "    tau_zl[4] = tau_lm\n"
  "\n"
  "    var_zl[0] = 1/numchan(Nl) * lm_inf^2*lh_inf*(1-lm_inf)^2*(1-lh_inf)\n"
  "    var_zl[1] = 2/numchan(Nl) * lm_inf^3*lh_inf*(1-lm_inf)*(1-lh_inf)\n"
  "    var_zl[2] = 1/numchan(Nl) * lm_inf^4*lh_inf*(1-lh_inf)\n"
  "    var_zl[3] = 1/numchan(Nl) * lm_inf^2*lh_inf^2*(1-lm_inf)^2\n"
  "    var_zl[4] = 2/numchan(Nl) * lm_inf^3*lh_inf^2*(1-lm_inf)\n"
  "\n"
  "	FROM i=0 TO 4 {\n"
  "		mu_zl[i] = exp(-dt/tau_zl[i])\n"
  "        noise_zl[i] = sqrt(var_zl[i] * (1-mu_zl[i]*mu_zl[i])) * normrand(0,NF)\n"
  "		zl[i] = zl[i]*mu_zl[i] + noise_zl[i]\n"
  "    }	\n"
  "\n"
  ":T Ca (m)\n"
  "    tau_tm =( 4 / ( exp((vm+25)/20) + exp(-(vm+100)/15) ) ) \n"
  "	tm_inf = 1 / ( 1 + exp(-(vm+50)/7.4) )\n"
  ":T Ca (h)	\n"
  "\n"
  "    tau_th = ( 86/ ( exp((vm+46)/4) + exp(-(vm+405)/50) ) ) \n"
  "	th_inf = 1.0 / ( 1 + exp((vm+78)/5) )\n"
  "\n"
  "    tau_zt[0] = tau_tm*tau_th/(tau_tm+2*tau_th)\n"
  "    tau_zt[1] = tau_tm*tau_th/(tau_tm+tau_th)\n"
  "    tau_zt[2] = tau_th\n"
  "    tau_zt[3] = tau_tm/2\n"
  "    tau_zt[4] = tau_tm\n"
  "\n"
  "    var_zt[0] = 1/numchan(Nt) * tm_inf^2*th_inf*(1-tm_inf)^2*(1-th_inf)\n"
  "    var_zt[1] = 2/numchan(Nt) * tm_inf^3*th_inf*(1-tm_inf)*(1-th_inf)\n"
  "    var_zt[2] = 1/numchan(Nt) * tm_inf^4*th_inf*(1-th_inf)\n"
  "    var_zt[3] = 1/numchan(Nt) * tm_inf^2*th_inf^2*(1-tm_inf)^2\n"
  "    var_zt[4] = 2/numchan(Nt) * tm_inf^3*th_inf^2*(1-tm_inf)\n"
  "\n"
  "	FROM i=0 TO 4 {\n"
  "		mu_zt[i] = exp(-dt/tau_zt[i])\n"
  "        noise_zt[i] = sqrt(var_zt[i] * (1-mu_zt[i]*mu_zt[i])) * normrand(0,NF)\n"
  "		zt[i] = zt[i]*mu_zt[i] + noise_zt[i]\n"
  "    }	\n"
  "	\n"
  "	zleak=normrand(0,NFleak)\n"
  "	m_exp = 1 - exp(-dt/tau_m)\n"
  "	h_exp = 1 - exp(-dt/tau_h)\n"
  "	n_exp = 1 - exp(-dt/tau_n)\n"
  "	km_exp= 1 - exp(-dt/tau_km)	\n"
  "	lm_exp= 1 - exp(-dt/tau_lm)	\n"
  "	lh_exp= 1 - exp(-dt/tau_lh)	\n"
  "	tm_exp= 1 - exp(-dt/tau_tm)	\n"
  "	th_exp= 1 - exp(-dt/tau_th)	\n"
  "	\n"
  "   UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "    if (fabs(exp(x/y) - 1) < 1e-6) {\n"
  "        vtrap = y*(1 - x/y/2)\n"
  "    }else{\n"
  "        vtrap = x/(exp(x/y) - 1)\n"
  "    }\n"
  "}\n"
  "FUNCTION mulnoise(mean,sd,power){\n"
  "	LOCAL i,avg\n"
  "	avg=1\n"
  "	FROM i=1 TO power {\n"
  "		avg=avg*normrand(mean,sd)\n"
  "	}\n"
  "	mulnoise=avg\n"
  "}\n"
  "\n"
  "FUNCTION numchan(Nchannels){\n"
  "	if (Nchannels>0){\n"
  "		numchan=(Nchannels)\n"
  "	}else{\n"
  "		numchan=1\n"
  "	}\n"
  "}\n"
  ;
#endif
