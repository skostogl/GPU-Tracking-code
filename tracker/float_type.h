#ifndef FLOAT_TYPE_H
#define FLOAT_TYPE_H

#define DOUBLE_PREC

#ifndef DOUBLE_PREC
  typedef float Tfloat;
#else
  typedef double Tfloat;
  //#undef DOUBLE_PREC
#endif

#endif
