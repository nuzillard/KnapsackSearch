This file reports what I did to facilitate the use of the python scripts
I developped.

In the following text, *HOME* stands for the name of the directory 
directly assessible from my session on Windows 10.
For me, *HOME* is C:\Users\jmn .

In *HOME*\Documents I created a directory named PythonLib
and I appended *HOME*\Documents\PythonLib to the value of the system variable PYTHONPATH.
You will have to create this variable if it does not exist already.

I copied all my python modules in directory PythonLib.
In this way, my python modules can be run from any current directory, using `python -m module_name` instead of
`python module_name.py`

I created the directory *HOME*\NMRShiftDB with nmrshiftdb2\, Demo.class, predictorc.jar, and predictSdf.bat inside.
This version of predictSdf.bat works only if the value of the new system variable NMRShiftDB
is set to *HOME*\NMRShiftDB and if *HOME*\NMRShiftDB is added to the
value of the system variable PATH. 
The value of a new system variable MY_JAVA is set to the full path to java.exe,
presently for me, *HOME*\Documents\jdk1.8.0_202\bin\java.exe , but may be simply java for a regular installation of the Java
runtime environment.
In this way predictSdf.bat can be run from any current directory as its execution depends on the location
of the java interpreter and of the nmrshiftdb2-related files.
