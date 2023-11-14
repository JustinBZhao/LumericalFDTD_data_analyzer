This is a plugin library for processing Lumerical datasets with ease.

For more information on Lumerical dataset, see [here](https://optics.ansys.com/hc/en-us/articles/360034409554-Introduction-to-Lumerical-datasets).

This library uses ["redblueTecplot"](https://www.mathworks.com/matlabcentral/fileexchange/69559-diverging-redblue-colormap-from-tecplot) from MATLAB file exchange. Copyright (c) 2018, Fernando Zigunov.
All rights reserved.

## Handle class

This class is a handle class. It means that if you have a class object ``obj``, and you directly assign the object to another variable:
```MATLAB
another_obj = obj;
```
``another_obj`` will be referring to the same underlying object. If you modify the underlying data in ``obj`` (let's say it has a data member ``data``:
```MATLAB
obj.data = 5;
```
Then ``another_obj.data`` will also become ``5``.
If you need a deep copy of the underlying data into another variable, you should do:
```MATLAB
another_obj = obj.copy();
```
## Dataset
In Lumerical, dataset is a convenient object bundling attributes data and its corresponding parameter data together. There are two types of dataset:
- Matrix dataset, which is the regular version of the dataset.
- Rectilinear dataset, which always includes x,y,z parameters, which means that the attributes are spatially dependent.

This MATLAB library handles the two types of dataset in a similar manner. x,y,z parameters behave very similar to other parameters.

## Getting started
To use this library, you can clone this repository (assuming you know how to use Git) or simply download the latest version of the files inside this repository. __'src'__ folder contains all the source codes needed to use the library -- simply class definitions. __'tests'__ folder contains testing related files. They are mainly for me to test the code, but if you want, you can test them on your machine too.
When you use this library in your script, make sure that the __'src'__ folder is in your MATLAB path. This can be done in a few ways.
- You can release all source codes from the __'src'__ folder into your working directory (current MATLAB folder). I don't recommend doing this, because it makes your working folder messy.
- You can move the __'src'__ folder into your working directory or a sub-folder in your working directory. When running the code, make sure to add the path of the __'src'__ folder to the MATLAB path using ``addpath`` [function](https://www.mathworks.com/help/matlab/ref/addpath.html). In this case, you can use the relative path to your working directory, for example:
```MATLAB
addpath("\src");
```
if the __'src'__ folder is directly under your working directory.
- You can move the __'src'__ folder into a location on your machine. In this case, you need to include the full path when using ``addpath``, for example:
```MATLAB
addpath("D:\[your_path]\src");
```
## Load dataset into MATLAB class
To work with this library, you first need to convert the MATLAB dataset exported from Lumerical into the specific format defined by the LumericalDataset class. This is very simple. Suppose you loaded the .mat file (containing the data) into MATLAB workspace, you should have a ``struct`` variable (``dataset``) representing one of your datasets. Simply do:
```MATLAB
converted_dataset = LumericalDataset.createObject(dataset);
```
And ``converted_dataset`` is the converted result. Now you can perform various operations on this object.
### Error while loading dataset
Occasionally, you might see that the dataset was not converted, instead, an error was thrown. the ``createObject`` function performs extensive checking to ensure that the target being converted is indeed a dataset. Otherwise, it is carefully (hopefully) written to report specific errors regarding why the target cannot be recognized as a dataset. In almost all situations, the dataset directly exported from Lumerical should pass all the checkings. However, there are very rare situations where I don't think the dataset has any value to be further studied, even though it is legal in Lumerical (for example, a dataset with no attribute), so I reject it.
These errors are here to help you, and to make sure that various assumtions (invariants) about the data inside the class holds. Don't be scared away.
## Perform operations on the converted dataset
Once you have successfully converted the dataset, the good things begin: there are many operations (class methods) you can perform on the dataset, including making 1D and 2D plots, merging datasets, etc. For example, 1D and 2D visualizations inside Lumerical visualizer can also be done here in MATLAB, but with prettier looking figures and better interactivity. You can also manipulate data in the dataset in several useful ways, which you CAN'T do in Lumerical.

