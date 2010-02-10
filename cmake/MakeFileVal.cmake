MACRO(MAKEFILE_VAL file name varname)

FILE(WRITE make.stub "default:")
FILE(APPEND make.stub "	echo ${name} >make.out")
FILE(APPEND make.stub "")
FILE(APPEND make.stub "include " ${file})

# Should use ${CMAKE_MAKE_PROGRAM}, or is that set to
# something other than nmake when using visual studio?
IF(WIN32)
  EXEC_PROGRAM(nmake ARGS /f make.stub)
ELSE(WIN32)
  EXEC_PROGRAM(make ARGS -f make.stub)
ENDIF(WIN32) 

FILE(READ make.out ${varname})
FILE(REMOVE make.stub make.out)

ENDMACRO(MAKEFILE_VAL)
