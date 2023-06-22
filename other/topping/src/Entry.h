#ifndef ENTRY_H
#define ENTRY_H
#include <stdint.h>
#include "jpsp_oracle.h"

struct xyLoc {
  int16_t x;
  int16_t y;
};

void PreprocessMap(std::vector<bool> &bits, int width, int height, const char *filename);
void *PrepareForSearch(std::vector<bool> &bits, int width, int height, const char *filename);
void GetPath(void *data, xyLoc s, xyLoc g, std::vector<xyLoc> &path, warthog::jpsp_oracle& oracle);//, int &callCPD);
const char *GetName();

#endif
