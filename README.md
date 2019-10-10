# GRD - Geometric Relation Distribution
#### Copyright (C) 2018 Dario Lodi Rizzini.

OVERVIEW
-------------------------------------------------

Library **grd** implements the Geometric Relation Distribution, a signature 
encoding the spatial distribution of a planar point sets that can be used 
to address loop closure in 2D LiDAR localization and mapping. 
The library has been kept to a minimal design. 

If you use this library, please cite the following paper: 

<!--
D. Lodi Rizzini, F. Galasso, S. Caselli. 
Geometric Relation Distribution for Place Recognition,
IEEE Robotics and Automation Letters (RA-L), Accepted paper, 2019. 
-->

```  
   @article{lodirizzini2019ral,
     author={Lodi Rizzini, D. and Galasso, F. and Caselli, S.},
     title={{Geometric Relation Distribution for Place Recognition}},
     journal={IEEE Robotics and Automation Letters (RA-L)},
     volume={4},
     issue={2},
     pages={523-529},
     year={2019},
     publisher={IEEE},
     issn = {2377-3766},
     doi = {10.1109/LRA.2019.2891432},
     note = {DOI 10.1109/LRA.2019.2891432},
   }  
```  


or the most relevant associated publications by visiting: 
http://rimlab.ce.unipr.it/


DEPENDENCIES
-------------------------------------------------

The software depends on the following external libraries

- Boost (submodule lexical_cast)
- Eigen 3.0 

Other dependencies are placed in directory thirdparty. 
Some examples require the external application "gnuplot" to display 
results. 


HOW TO COMPILE
-------------------------------------------------

Let ${grd_ROOT} be the install directory of your local copy 
of library grd. 
The following standard commands are required to compile it:

-  cd ${grd_ROOT}
-  mkdir build
-  cd build
-  cmake ..
-  make

You can also install the library into a system directory. 
To change the install directory you must set cmake environment
variable ${CMAKE_INSTALL_PREFIX} (e.g. using command "ccmake .."
before calling "cmake .."). 
Its default value on UNIX-like/Linux systems is "/usr/local".
After compiling library grd, run the command:

-  sudo make install

The command "sudo" is required only if ${CMAKE_INSTALL_PREFIX} 
is a system diretory managed by administrator user root.
Such command copies:
- header files of ${grd_ROOT}/include/grd to
   ${CMAKE_INSTALL_PREFIX}/include/grd/
- library files ${grd_ROOT}/lib/libgrd.a to
   ${CMAKE_INSTALL_PREFIX}/lib/
- cmake script ${grd_ROOT}/cmake_modules/grdConfig.cmake to
   ${CMAKE_INSTALL_PREFIX}/share/grd/


HOW TO USE LIBRARY grd IN YOUR PROJECT
-------------------------------------------------

If library grd has been installed in system directory "/usr/local",
then it is straighforward to use it in your projects.
You need to add the following lines to your project as in this example:


> CMAKE_MINIMUM_REQUIRED(VERSION 2.8)  
> PROJECT(foobar)  
> 
> find_package(grd REQUIRED)  
> message(STATUS "grd_FOUND ${grd_FOUND}")  
> message(STATUS "grd_INCLUDE_DIRS ${grd_INCLUDE_DIRS}")  
> message(STATUS "grd_LIBRARY_DIRS ${grd_LIBRARY_DIRS}")  
> message(STATUS "grd_LIBRARIES ${grd_LIBRARIES}")  
>  
> if(${grd_FOUND})   
>   include_directories(${grd_INCLUDE_DIRS})  
>   link_directories(${grd_LIBRARY_DIRS})  
> endif()  
> 
> add_executable(foobar foobar.cpp)  
> target_link_libraries(foobar ${grd_LIBRARIES})  

The above example uses the variables defined in grdConfig.cmake:

-  grd_FOUND - system has grd module
-  grd_INCLUDE_DIRS - the grd include directories
-  grd_LIBRARY_DIRS - the grd library directories
-  grd_LIBRARIES - link these to use grd



