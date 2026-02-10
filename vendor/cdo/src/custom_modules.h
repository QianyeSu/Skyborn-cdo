#ifndef CUSTOM_MODULES_H
#define CUSTOM_MODULES_H

#ifdef CUSTOM_MODULES
void load_custom_module(std::string path);
void load_custom_modules(std::string folder_path);
void close_library_handles();
#endif

#endif
