#!/usr/bin/env python

import sys, os

if len(sys.argv) != 2:
    print 'Usage: ./create_test.py classname'

classname = sys.argv[1]
print classname
testfilename = classname+'-test.cpp'
print testfilename
testfile = open(testfilename, 'w')

template = open('TEMPLATE')
template = template.read()
template = template.replace('CLASSNAME', classname)
template = template.replace('FILENAME', testfilename)

testfile.write(template)

# adding the new unit test to svn
os.system("svn add " + testfilename)
# setting svn properties for new unit test
os.system("svn propset svn:keywords Id " + testfilename)

cmakefile = open('CMakeLists.txt', 'r')
cmake = cmakefile.read()
cmakefile.close()

srcs = 'SET(SRCS_'+classname.upper()+' '+testfilename+')'
print srcs
cmake = cmake.replace('#### Sources', '#### Sources\n'+srcs)


test = 'ADD_MGFP_TEST("'+classname+'" test_'+classname.lower()+' ${SRCS_'+classname.upper()+'})'
print test
cmake = cmake.replace('#### Tests', '#### Tests\n'+test)

cmakefile = open('CMakeLists.txt', 'w')
cmakefile.write(cmake)
cmakefile.close()
