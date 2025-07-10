#include <stdint.h>
#include <string.h>
#include "icd_fields.h"
#include "nav_words.h"

/*
 * Placeholder implementation that assembles a D1 subframe using
 * the provided fields. Only the first two words are constructed
 * to demonstrate the new build_word() helper. Remaining words are
 * left zero so that existing functionality is preserved.
 */
void assemble_d1_subframe(int sfid, const B1I_D1_Frame *frm,
                          int week, uint32_t sow, uint32_t words[10])
{
    memset(words, 0, sizeof(uint32_t) * 10);
    if (!frm) return;

    if (sfid == 1) {
        uint32_t w1_payload = (0x712 << 3 | sfid) << 8 | ((sow >> 12) & 0xFF);
        words[0] = build_word(w1_payload, 26);
        uint32_t w2_payload = ((sow & 0xFFF) << 13) | (week & 0x1FFF);
        words[1] = build_word(w2_payload, 26);
    }
}
