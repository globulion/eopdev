Contributing to oep-dev
=======================

OepDev is a plugin to Psi4. Therefore it should follow the programming etiquette of Psi4. Also,
oep-dev has additional programming tips to make the code more versatile and easy in further development.
Here, I emphasise on most important aspects regarding the **programming rules**:

1. Main routine and libraries
-----------------------------

Oep-dev has only *one* source file in the plugin base directory, i.e., `main.cc`. This is the main
driver routine that handles the functionality of the whole OEP testing platform: specifies options for 
Psi4 input file and implemented test routines based on the options. Other sources are stored
in `oepdev/libNAME*` directories where `NAME` is the name of the library with sources and header files.

> Things to remember:
>
>   1. It is *not* permitted to place any new source files in the plugin base directory (i.e., where `main.cc` resides).
>   2. Any additional source has to be placed in `oepdev/libNAME*` directory (either existing one or a new one; in the 
>      latter case remember to add the new `*.cc` files to `CMakeLists.txt` in the plugin base directory.


2. Header files in libraries
----------------------------

Header files are handy in obtaining a quick glimpse of the functionality within certain library. Each library
directory should contain at least one header file in oep-dev. However, header files can be problematic if not managed properly. 

> Things to remember:
>
>   1. **Header preprocessor variable**. Define the preprocessor variable specyfying the existence of include 
>      of the particular header file. The format of such is
>      ```c++
>      #ifndef MODULE_LIBRARY_HEADER_h
>      #define MODULE_LIBRARY_HEADER_h
>      // other code
>      #endif // MODULE_LIBRARY_HEADER_h
>      ```
>      Last line is the **end** of the header file. The preprocessor variables represents
>      the directory tree `oepdev/MODULE/LIBRARY/HEADER.h' structure (where `oepdev` is the base plugin directory). 
>        * `MODULE` is the plugin module name (e.g. `oepdev`, the 
>           name of the module directory)
>        * `LIBRARY` is the name of the library (e.g. `libutil`, should be the same as library directory name)
>        * `HEADER` is the name of the header in library directory (e.g. `diis` for `diis.h` header file)
>   2. **Set local namespace**. To prevent naming clashes with other modules and with Psi4 remember to specify
>      *only one* the local namespace of the format:



*********
