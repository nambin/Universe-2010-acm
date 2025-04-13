#pragma warning(disable:4251)

#ifdef _WIN32
#ifdef NATURE_SOURCE
#define NATUREDLL __declspec(dllexport)
#else
#define NATUREDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define NATUREDLL
#define DLLEXPORT
#endif