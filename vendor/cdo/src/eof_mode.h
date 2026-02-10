#ifndef EOF_MODE_H
#define EOF_MODE_H

enum T_WEIGHT_MODE
{
  WEIGHT_OFF,
  WEIGHT_ON
};
enum T_EIGEN_MODE
{
  JACOBI,
  DANIELSON_LANCZOS
};

enum T_EIGEN_MODE get_eigenmode(void);
enum T_WEIGHT_MODE get_weightmode(void);

#endif
