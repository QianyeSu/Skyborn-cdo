#ifndef REFERENCE_TO_GRID_H
#define REFERENCE_TO_GRID_H

struct RefGrid
{
  int gridID = -1;
  bool exists = false;    // true: a reference exists.
  bool isValid = false;   // true: the reference could be dereferenced.
  bool notFound = false;  // true: the reference was not found.
};

// int referenceToGrid(int gridID);
RefGrid dereferenceGrid(int gridID);

#endif
