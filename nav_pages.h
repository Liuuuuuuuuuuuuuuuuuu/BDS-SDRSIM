#ifndef NAV_PAGES_H
#define NAV_PAGES_H
#include <stdint.h>
#include "icd_fields.h"
void assemble_d1_subframe(int sfid, const B1I_D1_Frame *frm,
                          int week, uint32_t sow, uint32_t words[10]);
void assemble_d1_subframe2(const B1I_D1_Frame *frm, uint32_t words[10]);
void assemble_d1_subframe3(const B1I_D1_Frame *frm, uint32_t words[10]);
#endif
