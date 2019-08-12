#include <util.h>

void fill_int(int *dst, int nbr, int val) {
  if (nbr % 1 != 0) {
    *dst++ = val;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = val;
    *dst++ = val;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = val;
    dst[1] = val;
    dst[2] = val;
    dst[3] = val;
    dst += 4;
  }
}

void fill_dbl(double *dst, int nbr, double val) {
  if (nbr % 2 != 0) {
    *dst++ = val;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = val;
    *dst++ = val;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = val;
    dst[1] = val;
    dst[2] = val;
    dst[3] = val;
    dst += 4;
  }
}

void fill_float(float *dst, int nbr, float val) {
  if (nbr % 2 != 0) {
    *dst++ = val;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = val;
    *dst++ = val;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = val;
    dst[1] = val;
    dst[2] = val;
    dst[3] = val;
    dst += 4;
  }
}

void fill_char(char *dst, int nbr, char val) {
  if (nbr % 2 != 0) {
    *dst++ = val;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = val;
    *dst++ = val;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = val;
    dst[1] = val;
    dst[2] = val;
    dst[3] = val;
    dst += 4;
  }
}

void acopy_int(const int *src, int *dst, int nbr) {
  if (nbr % 2 != 0) {
    *dst++ = *src++;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = *src++;
    *dst++ = *src++;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst += 4;
    src += 4;
  }
}

void acopy_dbl(const double *src, double *dst, int nbr) {
  if (nbr % 2 != 0) {
    *dst++ = *src++;
  }

  nbr >>= 1;

  if (nbr % 2 != 0) {
    *dst++ = *src++;
    *dst++ = *src++;
  }

  nbr >>= 1;

  while (nbr--) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst += 4;
    src += 4;
  }
}
