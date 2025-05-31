# FireFly

FireFly is a C project structured with `include` and `src` folders.  
This repository contains the source code and a PowerShell build script (`build.ps1`) to compile the project.
The results are saved as a series of .vtk files which may be visualized using ParaView.

The project is **heavily inspired by and references** *Low-Speed Aerodynamics* (2nd Edition) by **Katz and Plotkin**, specifically **Program No. 16: Unsteady Rectangular Lifting Surface (Vortex Lattice Method)** found in **Appendix D.3**.  
It serves as a modern interpretation and exploration of unsteady aerodynamic modeling based on the concepts in that text.

## Build Instructions for Windows Users

To compile this project on Windows, you must download and install the OpenBLAS Windows binary (ZIP), and add the path of the OpenBLAS bin directory to your system's PATH variables. The provided PowerShell build script assumes that OpenBLAS has been extracted to the C:\ drive. To compile, run from PowerShell:

```powershell
.\build.ps1
