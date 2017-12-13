Contributing to oep-dev
=======================

OepDev is a plugin to Psi4. Therefore it should follow the programming etiquette of Psi4. Also,
oep-dev has additional programming tips to make the code more versatile and easy in further development.
Here, I emphasise on most important aspects regarding the **programming rules**:

Main routine and libraries
--------------------------

Oep-dev has only *one* source file in the plugin base directory, i.e., `main.cc`. This is the main
driver routine that handles the functionality of the whole OEP testing platform: specifies options for 
Psi4 input file and implemented test routines based on the options. Other sources are stored
in `oepdev/libNAME*` directories where `NAME` is the name of the library with sources and header files.

Things to remember:

  1. **No other sources in base directory.** 
     It is not permitted to place any new source or other files in the plugin base directory 
     (i.e., where `main.cc` resides). 
  2. **Sources in library directories.** 
     Any additional source code has to be placed in `oepdev/libNAME*` directory (either existing one or a new one; in the 
     latter case remember to add the new `*.cc` files to `CMakeLists.txt` in the plugin base directory.
  3. **Miscellanea in special directories.** 
     If you want to add documentation, put it in the `doc` directory. 
     If you want to add graphics, put it in the `images` directory.


Header files in libraries
-------------------------

Header files are handy in obtaining a quick glimpse of the functionality within certain library. Each library
directory should contain at least one header file in oep-dev. However, header files can be problematic if not managed properly. 

Things to remember:

   1. **Header preprocessor variable**. Define the preprocessor variable specyfying the existence of include 
      of the particular header file. The format of such is
      ```c++
      #ifndef MODULE_LIBRARY_HEADER_h
      #define MODULE_LIBRARY_HEADER_h
      // rest of your code goes here
      #endif // MODULE_LIBRARY_HEADER_h
      ```
      Last line is the **end** of the header file. The preprocessor variables represents
      the directory tree `oepdev/MODULE/LIBRARY/HEADER.h` structure (where `oepdev` is the base plugin directory). 
        * `MODULE` is the plugin module name (e.g. `oepdev`, the 
           name of the module directory)
        * `LIBRARY` is the name of the library (e.g. `libutil`, should be the same as library directory name)
        * `HEADER` is the name of the header in library directory (e.g. `diis` for `diis.h` header file)
   2. **Set local namespace**. To prevent naming clashes with other modules and with Psi4 remember to specify
      *only one* local namespace per library of the following format:
      ```c++
      namespace MODULE_LIBRARY {
      // your code goes here
      } // EndNameSpace MODULE_LIBRARY
      ```
      Example of namespace is `oepdev_liboep` in module `oepdev` and library `liboep`. 

Documenting the code
--------------------

Code has to be documented (at best at a time it is being created). The place for documentation 
is always in header files. Additional documentation can be also placed in source files. Leaving a chunk of code
for a production run without documentation is unacceptable. 

Use Doxygen style for documentation all the time. Remember that it supports markdown which can make the documentation
even more clear and easy to understand.
Additionally you can create a nice `.rst` documentation file for Sphinx program.

Things to remember:

   1. **Descriptions of classes, structures, global functions, etc**. Each programming object should have a description.
   2. **Documentation for function arguments and return object**. 
      Usage of functions and class methods should be explained by providing the description of all arguments 
      (use `\param` and `\return` Doxygen keywords).
   3. **One-line description of class member variables**. Any class member variable should be preceded by 
      a one-liner documentation (starting from `///`).
   4. **Do not be afraid of long names in the code**. Self-documenting code is a bless!

Naming conventions
------------------

Naming is important because it helps to create more readable and clear self-documented code. 

Things to remember:

   1. **Do not be afraid of long names in the code, but avoid redundancy**. Examples of good and bad names:
      * good name: `get_density_matrix`; bad name: `get_matrix`. Unless there is only one type of matrix
        a particular objects can store, `matrix` is not a good name for a getter method. 
      * good name: `class Wavefunction`, bad name: `class WFN`
      * good name: `int numberOfErrorVectors`, bad name: `int nvec`, bad name: `the_number_of_error_vectors`
      * good name: `class EFPotential`, probably bad name: `class EffectiveFragmentPotential`.
        The latter might be understood by some people as a class that inherits from `EffectiveFragment` class. 
        If it is not the case, compromise between abbreviation and long description is OK.
   2. **Short names are OK in special situations**. In cases meaning of a particular variable is obvious and
      it is frequently used in the code locally, it can be named shortly. Examples are:
      * `i` when iterating
      * `no` number of occupied orbitals, `nv` number of virtual orbitals, etc.
   3. **Clumped names for variables and dashed names for functions**. Try to distinguish between variable name 
      like `sizeOfOEPTypeList` and a method name `get_matrix()` (neither `size_of_OEP_type_list`, nor `getMatrix()`).
      This is little bit cosmetics, but helps in managing the code when it grows.
   4. **Class names start from capital letter**. However, avoid only capital letters in class names, unless it is obvious.
      Avoid also dashes in class names (they are reserved for global functions and class methods). Examples: 
      * good name: `DIISManager`, bad name: `DIIS`.
      * good name: `EETCouplingSolver`, bad name: `EETSolver`, very bad: `EET`.
*********
