#pragma warning(disable:4251)

#if defined (_WIN32) || (_WIN64)
#ifdef BASIS_SOURCE
#define BASISDLL __declspec(dllexport)
#else
#define BASISDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define BASISDLL
#define DLLEXPORT
#endif

