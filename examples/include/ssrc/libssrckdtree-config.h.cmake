#ifndef _LIBSSRCKDTREE_CONFIG_H
#define _LIBSSRCKDTREE_CONFIG_H

/* Enable debug code. */
/* #undef LIBSSRCKDTREE_DEBUG */

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_DLFCN_H 

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_INTTYPES_H  

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_MEMORY_H  

/* Define if the compiler implements namespaces. */
#cmakedefine LIBSSRCKDTREE_HAVE_NAMESPACES

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_STDINT_H  

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_STDLIB_H  

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_STRINGS_H  

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_STRING_H  

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_SYS_STAT_H  

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_SYS_TYPES_H  

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine LIBSSRCKDTREE_HAVE_UNISTD_H  

/* Name of package */
#cmakedefine LIBSSRCKDTREE_PACKAGE  "libssrckdtree" 

/* Define to the address where bug reports for this package should be sent. */
#cmakedefine LIBSSRCKDTREE_PACKAGE_BUGREPORT   

/* Define to the full name of this package. */
#cmakedefine LIBSSRCKDTREE_PACKAGE_NAME  "libssrckdtree" 

/* Define to the full name and version of this package. */
#cmakedefine LIBSSRCKDTREE_PACKAGE_STRING  "libssrckdtree 1.0.2" 

/* Define to the one symbol short name of this package. */
#cmakedefine LIBSSRCKDTREE_PACKAGE_TARNAME  "libssrckdtree" 

/* Define to the version of this package. */
#cmakedefine LIBSSRCKDTREE_PACKAGE_VERSION  "1.0.2" 

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine LIBSSRCKDTREE_STDC_HEADERS  

/* Version number of package */
#cmakedefine LIBSSRCKDTREE_VERSION  "1.0.2" 


#if !defined(NS_KD)
#  define NS_KD kd
#endif

#if !defined(KD_DEFINE_NAMESPACE)
#  define KD_DEFINE_NAMESPACE(name) NS_KD::name
#endif

#if !defined(NS_KD_DECL_PREFIX)
#  define NS_KD_DECL_PREFIX \
namespace kd {
#endif

#if !defined(NS_KD_DECL_SUFFIX)
#  define NS_KD_DECL_SUFFIX \
}
#endif

 
/* _LIBSSRCKDTREE_CONFIG_H */
#endif
