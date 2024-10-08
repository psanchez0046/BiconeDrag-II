
----------------------------------------
----------------------------------------
Readme file for BiconeDrag
----------------------------------------
----------------------------------------

========================================
MATLAB (GNU Octave) package codes:
----------------------------------------
This package contains four files:
- "example_script_bicone.m": is an example of a script that calls postprocessing_Bicone_CPC
- "postprocessingBiconeCPC.m": is the main function and calculates the rheological properties from amplitude ratio and phase lag experimental data of a rotational rheometer with a bicone fixture.
- "solve_NS_bicono.m": is a function that solves the Navier-Stokes equations with no-slip and Boussinesq_Scriven boundary conditions with second order centered finite differences for the Bicone bob-cylindrical cup configuration.
- "GetFilenames.m": is a function that returns a cell array with the names of the files having the specific pattern contained in the specified filepath.

REQUIREMENTS:
-------------
Check that MATLAB (or GNU Octave) is installed on your computer.
You can download the latest released version of GNU Octave for your operating system from: https://www.gnu.org/software/octave/download or https://ftp.gnu.org/gnu/octave/.

EXECUTION:
----------
1. All subroutines and files must be in the Current Folder or in a folder contained in the Search Path.
2. Configure parameters in your_script_file by editing "example_script_bicone.m" and changing the input parameters.
3. Execute the script file as follows:
	- On MATLAB GUI (or Octave GUI) version:
		- type on Command Window: run your_script_file.
		- open your_script_file.m on the editor and click run (or press F5)
	- On command line open cmd (Windows) or Terminal (Linux and OS X) and type:
		- matlab -r your_script_file.
		- octave your_script_file.

		
========================================
Python3 package codes:
----------------------------------------
This package contains two files:
- "example_script_bicone_3.py" : is an example of a script that calls postprocessing_Bicone_CPC contained in bicfuns_3 module.
- "bicfuns_3.py": a module that contains the definition of solve_NS_bicono and postprocessingBiconeCPC functions.

REQUIREMENTS:
-------------
Check that Python3 with numpy and scipy modules is installed on your computer.
You can download all versions of Python for your operating system from: https://www.python.org/downloads/.

EXECUTION:
----------
1. All subroutines must be in the current folder or in a folder contained in the system path.
2. Configure parameters in your_script_file by editing example_script_bicone_3.py and changing the input parameters.
3. Execute the script file: open cmd (Windows) or Terminal (Linux and OS X) and type: python your_script_file.

TESTING THE PACKAGE:
--------------------
You can use "test_exp.txt" as an experiment file for testing the package by executing the original example_script_bicone.
This experiment file has been generated using a pure viscous interface \eta_s^* = \eta_s' with values |\eta_s^*| varying from 1e-6 to 1 in 12 logarithmically spaced steps.
