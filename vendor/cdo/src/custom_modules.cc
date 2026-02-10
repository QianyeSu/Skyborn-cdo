#ifdef CUSTOM_MODULES
/***
  loads all custom modules in a specified folder
  @param folder_path custom modules folder
*/
void
load_custom_modules(std::string const &folder_path)
{
  DIR *dir = opendir(folder_path.c_str());
  std::string file_path;
  std::regex library_regex("(.*\\.so)");
  if (dir != nullptr)
    {
      struct dirent *ent = readdir(dir);
      while (ent != nullptr)
        {
          if (std::regex_match(ent->d_name, library_regex))
            {
              file_path = folder_path + "/" + ent->d_name;
              load_custom_module(file_path);
            }
          ent = readdir(dir);
        }
    }
  else { std::cerr << "Could not find " << folder_path << "for loading custom modules" << std::endl; }
}

/***
 * Loads a custom module from given path.
 * Modules must contain a (TODO: rename function) init_custom_module function
 * Program exits if a module could not be loaded.
 * @param module file path
 */
void
load_custom_module(std::string const &file_path)
{
  void (*init_custom_module)();
  void *lib_handle = dlopen(file_path.c_str(), RTLD_LAZY);
  custom_modules_lib_handles.push_back(lib_handle);
  if (!lib_handle)
    {
      std::cerr << "Cannot open library: " << dlerror() << std::endl;
      return;
    }

  dlerror();
  init_custom_module = (void (*)()) dlsym(lib_handle, "init_custom_module");
  const char *dlsym_error = dlerror();

  if (dlsym_error)
    {
      std::cerr << "Cannot load symbol 'init_custom_module': " << dlsym_error << std::endl;
      dlclose(lib_handle);
      return;
    }

  init_custom_module();
  std::cout << "loaded custom module from '" << file_path << "'" << std::endl;
}

/***
  closes the handles for the loaded custum modules
*/
void
close_library_handles()
{
  for (unsigned long i = 0, n = custom_modules_lib_handles.size(); i < n; ++i) { dlclose(custom_modules_lib_handles[i]); }
}
#endif
