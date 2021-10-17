#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _BcSyn_reg(void);
extern void _CaDynamics_reg(void);
extern void _HHst_reg(void);
extern void _Kv3_1_reg(void);
extern void _ca_ch_reg(void);
extern void _cad_reg(void);
extern void _cagk_reg(void);
extern void _cal2_reg(void);
extern void _can_reg(void);
extern void _cap_ch_reg(void);
extern void _caq_reg(void);
extern void _gradsyn_reg(void);
extern void _kdr_reg(void);
extern void _nav1p8_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," mods/BcSyn.mod");
    fprintf(stderr," mods/CaDynamics.mod");
    fprintf(stderr," mods/HHst.mod");
    fprintf(stderr," mods/Kv3_1.mod");
    fprintf(stderr," mods/ca_ch.mod");
    fprintf(stderr," mods/cad.mod");
    fprintf(stderr," mods/cagk.mod");
    fprintf(stderr," mods/cal2.mod");
    fprintf(stderr," mods/can.mod");
    fprintf(stderr," mods/cap_ch.mod");
    fprintf(stderr," mods/caq.mod");
    fprintf(stderr," mods/gradsyn.mod");
    fprintf(stderr," mods/kdr.mod");
    fprintf(stderr," mods/nav1p8.mod");
    fprintf(stderr," mods/vecevent.mod");
    fprintf(stderr, "\n");
  }
  _BcSyn_reg();
  _CaDynamics_reg();
  _HHst_reg();
  _Kv3_1_reg();
  _ca_ch_reg();
  _cad_reg();
  _cagk_reg();
  _cal2_reg();
  _can_reg();
  _cap_ch_reg();
  _caq_reg();
  _gradsyn_reg();
  _kdr_reg();
  _nav1p8_reg();
  _vecevent_reg();
}
