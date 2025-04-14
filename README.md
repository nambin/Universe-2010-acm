# Liquid Simulation Libraries 

This repository contains C++ libraries used for the liquid simulation techniques presented in our ACM TOG 2010 publications [1] [2].

**Status:**
* This code was originally developed in a Windows Visual Studio 2010 environment.
* It has **not been actively maintained since 2010**.
* It was uploaded to GitHub in April 2025 primarily for archival and reference purposes. Compatibility with modern compilers and systems is unlikely without modification.

**Contents:**
Most implementations reside in the `_Header` and `_Source` directories. The libraries cover:
* **Basic Math & Graphics:** Vectors, BiTree, triangular mesh handling, color representation, XML/OBJ/PLY file handlers, Intel TBB (Threading Building Blocks).
* **Numerical Solvers:** Matrix solvers (including Conjugate Gradient Descent), Level Set methods, Computational Fluid Dynamics (CFD) solvers.
* **3D Rendering:** Photon map rendering implementation [3], Ray tracing, Realtime shader with GLSL

**References:**
* [1] Detail-Preserving Fully-Eulerian Interface Tracking Framework. ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2010), Vol. 29, Issue. 6, 176, 2010.
* [2] https://youtu.be/ToOJT6B7YMM
* [3] Realistic Image Synthesis Using Photon Mapping by Henrik Wann Jensen
