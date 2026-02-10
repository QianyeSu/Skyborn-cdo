#ifndef UTIL_FILEEXTENSIONS_H
#define UTIL_FILEEXTENSIONS_H

const char *filetypeext(int filetype);
void repl_filetypeext(char *file, const char *oldext, const char *newext);

#endif
