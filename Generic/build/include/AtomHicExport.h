
#ifndef ATOMHIC_EXPORT_H
#define ATOMHIC_EXPORT_H

#ifdef ATOMHIC_STATIC_DEFINE
#  define ATOMHIC_EXPORT
#  define ATOMHIC_NO_EXPORT
#else
#  ifndef ATOMHIC_EXPORT
#    ifdef AtomHic_EXPORTS
        /* We are building this library */
#      define ATOMHIC_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define ATOMHIC_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef ATOMHIC_NO_EXPORT
#    define ATOMHIC_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef ATOMHIC_DEPRECATED
#  define ATOMHIC_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef ATOMHIC_DEPRECATED_EXPORT
#  define ATOMHIC_DEPRECATED_EXPORT ATOMHIC_EXPORT ATOMHIC_DEPRECATED
#endif

#ifndef ATOMHIC_DEPRECATED_NO_EXPORT
#  define ATOMHIC_DEPRECATED_NO_EXPORT ATOMHIC_NO_EXPORT ATOMHIC_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef ATOMHIC_NO_DEPRECATED
#    define ATOMHIC_NO_DEPRECATED
#  endif
#endif

#endif /* ATOMHIC_EXPORT_H */
