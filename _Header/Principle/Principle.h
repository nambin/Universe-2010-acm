#pragma warning(disable:4251)

#ifdef _WIN32
#ifdef PRINCIPLE_SOURCE
#define PRINCIPLEDLL __declspec(dllexport)
#else
#define PRINCIPLEDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define PRINCIPLEDLL
#define DLLEXPORT
#endif
