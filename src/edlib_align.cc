#include "edlib_align.h"

int edlib_align_HW(char *query, int qlen, char *target, int tlen, int *start,
                   int *end, int k) {
  int ed = -1;
  EdlibAlignResult result =
      edlibAlign(query, qlen, target, tlen,
                 edlibNewAlignConfig(k, EDLIB_MODE_HW, EDLIB_TASK_LOC,
                                     additionalEqualities, 5));
  if (result.status == EDLIB_STATUS_OK) {
    ed = result.editDistance;
    if (ed >= 0) {
      *start = result.startLocations[0];
      *end = result.endLocations[0];
    }
  }
  edlibFreeAlignResult(result);
  return ed;
}