#  *****************************************************************************
#   File : zzz.R
#         ************************************************************
#   Description : 
#       Initialization of the library
#   Version : 3.1
#   Date : 2014-12-07
#         ************************************************************
#   Author : Julien Damon <julien.damon@free.fr>
#   License : LGPL
#   URL: http://julien.damon.free.fr
#  *****************************************************************************

#  *****************************************************************************
#   Title : .First.lib
#         ************************************************************
#   Description : 
#       Initialization of the library
#   Version : 3.1
#   Date : 2014-12-07
#  *****************************************************************************
".onAttach" <- function(lib, pkg)
{
  packageStartupMessage("far library : Modelization for Functional AutoRegressive processes\n")
  packageStartupMessage("version 0.6-4 (2014-12-07)\n")
}