#pragma warning(disable:4251)

#ifdef _WIN32
#ifdef ILLUMINATION_SOURCE
#define ILLUMINATIONDLL __declspec(dllexport)
#else
#define ILLUMINATIONDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define ILLUMINATIONDLL
#define DLLEXPORT
#endif