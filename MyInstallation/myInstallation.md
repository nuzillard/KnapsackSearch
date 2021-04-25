This file reports what I did to facilitate the use of the python scripts
I developped.

In the following text, *HOME* stands for the name of the directory 
directly assessible from my session on Windows 10.
For me, *HOME* is C:\Users\jmn

In *HOME*\Documents I created a directory named PythonLib
and I appended *HOME*\Documents\PythonLib to the value of the system variable PYTHONPATH.
You will have to create this variable if it does not exist already.

I copied all my python modules in this directory.
In this way, my python modules can be run from any current directory.

I created the directory *HOME*\NMRShiftDB with nmrshiftdb2\, Demo.class, predictorc.jar
and a modified version of predictSdf.bat inside. A copy of the latter is located besides of this file.
This version of predictSdf.bat works only if the value of system variable NMRShiftDB
is set to *HOME*\NMRShiftDB and if *HOME*\NMRShiftDB is added to the
value of the system variable PATH. 
The value of a new system variable MY_JAVA is set to the path to java.exe,
presently for me, *HOME*\Documents\jdk1.8.0_202\bin\java.exe .
In this way, predictSdf.bat can be run from any current directory.
