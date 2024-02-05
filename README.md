## General information
LumericalDataset library v 1.0
This is a MATLAB library for processing Lumerical datasets with ease.
For more information on Lumerical dataset, see [here](https://optics.ansys.com/hc/en-us/articles/360034409554-Introduction-to-Lumerical-datasets).
This library uses ["redblueTecplot"](https://www.mathworks.com/matlabcentral/fileexchange/69559-diverging-redblue-colormap-from-tecplot) from MATLAB file exchange. Copyright (c) 2018, Fernando Zigunov.
## Lumerical dataset
In Lumerical, datasets are structured data objects that collect a set of related matrices into a single convenient object. Dataset contains information about one or more physical quantities against some orthogonal coordinates. There are two types of dataset:
- Matrix dataset, which is the regular version of the dataset.
- Rectilinear dataset, which always includes x,y,z as 3 of the coordinates (plus other possible parameters too), which means that the attributes data are mapped spatially in 3D.
This library attempts to handle both types of dataset in a similar manner, by treating x,y,z the same way as other normal parameters.
## Handle class
The library converts the Lumerical Dataset into a class object. This class is a handle class, and you need to be careful when copying the object: 
+ Copy the handle
When you directly assign a class object ``obj`` to another object variable ``another_obj``, you are creating another handle that refers to the same underlying object.
```MATLAB
another_obj = obj;
```
If you make changes to the underlying data in ``obj`` (let's say it has a member ``a``:
```MATLAB
obj.a = 5;
```
Then ``another_obj.a`` will also become ``5``.
+ True deep copy
If you need create a different class object with the same underlying data values (make a deep copy), you should do this:
```MATLAB
another_obj = obj.copy();
```
## Getting started
To use this library, you can clone this repository (assuming you know how to use Git) or simply download the latest version of the necessary files inside this repository. __'src'__ folder contains all the source codes. Those are all you need to use the library. __'tests'__ folder contains testing related files. They are mainly for me to test the code, but if you want, you can test them in your MATLAB environment too.
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
### Terminologies
+ Parameter
A parameter describes one of the coordinate directions in a dataset. A dataset can have unlimited number of parameters (at least 1, however). A parameter must be a 1D vector.
Each parameter can be expanded to an interdependent parameter set. An interdependent set controls only one coordinate dimension, but with many alias names each associated with unique coordinate 1D vector values. An example is "lambda" and "f" (wavelength and frequency).
+ Position vector
Position vectors are the three special parameters - x, y and z. Every rectilinear dataset must have at least these three parameters.
+ Attribute
An attribute is the actual data contained in the multi-dimensional coordinate space described by the parameters. A dataset can have unlimited number of attributes (at least 1, however). An attribute usually has many dimensions.
+ Scalar vs. vector
Each attribute can be either scalar or vector. A vector attribute corresponds to an Euclidean vector. As a result, it need to store the magnitude of the value along all 3 Cartesian coordinates. Therefore, an additional dimension of length 3 is required for these attributes. In contrast, a scalar attribute does not have this extra dimension. 
+ Parameter slice index
Slice index is an important concept when dealing with making plots out of the dataset. If we have a dataset with 3 parameters, making a 2D plot means picking 2 out of 3 parameters to plot against. In this case, we need to specify how we want to slice through the remaining parameter, which is the slice index.
+ Attribute component
We cannot display all information of a vector attribute when making plots out of it, except for vector plots. Therefore, we always need to decide to plot along one Cartesian coordinate ("x", "y", or "z"), or the overall magnitude of the vector ("magnitude").
+ Scalar operation
In additional to an attribute being scalar or vector, its data can also be complex. This is very common for quantities of oscillating nature. To plot these complex values, select one of the following scalar operations: take real part ("real"), take imaginary part ("imag"), take modulus ("abs") or calculate the phase angle ("angle"). For a complex number x+yi, the angle is arctan(y/x).
## Load dataset into MATLAB class
To work with this library, you first need to convert the MATLAB dataset exported from Lumerical into the specific format defined by the LumericalDataset class. This is very simple. Suppose you loaded the .mat file (containing the data) into MATLAB workspace, you should have a ``struct`` variable (``dataset``) representing one of your datasets. Simply do:
```MATLAB
converted_dataset = LumericalDataset.createObject(dataset);
```
And ``converted_dataset`` is the converted result. Now you can perform various operations on this object.
### Error while loading dataset
Occasionally, you might see that the dataset was not converted, instead, an error was thrown. the ``createObject`` function performs extensive checking to ensure that the target being converted is indeed a dataset. Otherwise, it is carefully (hopefully) written to report specific errors regarding why the target cannot be recognized as a dataset. In almost all situations, the dataset directly exported from Lumerical should pass all the checkings. However, there are very rare situations where I don't think the dataset has any value to be further studied, even though it is legal in Lumerical (for example, a dataset with no attribute), so I will reject it.

These errors are here to help you, and to make sure that various assumtions (invariants) about the data inside the class holds. Don't be scared away.
## Perform operations on the converted dataset
Once you have successfully converted the dataset, the good things begin: there are many operations (class methods) you can perform on the dataset, including making 1D and 2D plots, merging datasets, etc. For example, 1D and 2D visualizations inside Lumerical visualizer can also be done here in MATLAB, but with prettier looking figures and better interactivity. You can also manipulate data in the dataset in several useful ways, which you CAN'T do in Lumerical.
## List of functions supported in this library
Below is the list of functions that have been implemented in this library, organized by category. However, I might update and add some new functionalities in the future versions.
### Data loading
+ [static] ``obj = LumericalDataset.createObject(lum_dataset)``
This is a static function. ``lum_dataset`` is the struct that is contained inside Lumerical's exported ``.mat`` file. This function will correctly resolve the type of the dataset (matrix dataset or rectilinear dataset) automatically. Any conversion that does not result in an error gives you a valid object.
+ [static] ``converted_obj = createObjectFromMat(mat_name, varargin)``
This is a static function. Instead of passing the struct object as an input argument, you can directly specify the name of the ``.mat`` file in ``mat_name``. You can also specify variables names to be loaded in additional arguments.
Note that in this case, ``converted_object`` is a struct containing each converted variable inside the original ``.mat`` file. Variables that are not Lumerical datasets are not converted. This is especially useful when the ``.mat`` file contains a few datasets, as well as some other non-dataset information.
+ ``obj = MatrixDataset(lum_dataset)`` or ``obj = RectilinearDataset(lum_dataset)`` (not recommended)
If you know which type of dataset you are converting (matrix dataset or rectilinear dataset), you can choose the correct syntax to manually perform the conversion. You **DO NOT** have to do it this way!
### Display information
+ ``showInformation(obj)``
This function will print the information of the dataset onto the command window. The information includes a list of parameters, attributes, etc.
### Get simple data
+ ``result = getParameterData(obj, parameter_name)``
Get the 1D vector data corresponding to a specific parameter. Note that for rectilinear datasets, to retrieve x, y, or z, you should directly call ``x = obj.x`` or similar instead of using this function.
+ ``result = getAttributeData(obj, attribute_name)`` (not recommended)
Get the data corresponding to a specific attribute. Use of this function is not recommended because the returned ``result`` will be a multi-dimensional matrix that is exactly how it is stored inside the dataset. Without knowing the details of the implementation of this library, the ``result`` is not useful.
### Set operations
+ ``setParameterSliceIndex(obj, varargin)``
Set the slice index of a parameter, including x, y, or z in rectilinear datasets. You can specify an arbitrary number of parameters, in the form of "parameter name-index" pair.
The slice index for each parameter is 1 by default. The value persists throughout the lifetime of the dataset object.
+ ``setAttributeComponent(obj, attribute_name, component)``
Set the component to plot for a vector attribute. ``component`` must be "x", "y", "z" or "magnitude". Cannot set this property for any scalar attribute (will receive an error).
The attribute component for each vector attribute is set to "magnitude" by default. And the selection of this property persists throughout the lifetime of the dataset object.
### Get data for plot
+ ``[xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name)``
Get the data for 1D plot. ``xdata`` is the parameter vector, ``ydata`` is the attribute data against that parameter, after slicing off all other dimensions. 
+ ``[xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name)``
Get the data for 2D plot. ``xdata`` is the parameter 1 vector, ``ydata`` is the parameter 2 vector, and ``zdata`` is the attribute data, against these two parameters and after slicing off all other dimensions.
+ ``[x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name)``
Get the data for 3D plot. ``x``, ``y`` and ``z`` correspond to 3 parameters, and ``data`` is the attribute after slicing off all other dimensions. You can only get data for 3D plot, but not actually plotting in 3D.
### Make plots
+ ``hPlot = plotData1D(obj, parameter_name, attribute_name, scalar_operation, ax)``
Make a 1D plot. ``scalar_operation`` and ``ax`` are optional positional arguments. The handle to the line object is returned as ``hPlot``.
+ ``[hSurf, hClb] = plotData2D(obj, parameter1_name, parameter2_name, attribute_name, scalar_operation, ax)``
Make a 2D plot. ``scalar_operation`` and ``ax`` are optional positional arguments. The handles to the surface object and the colorbar object are returned as ``hSurf`` and ``hClb``.
This 2D plot is actually a surface plot in 2D view. To make it more realistic, the coordinates are slightly adjusted. ``XData``, ``YData`` and ``ZData`` in the surface object (by calling ``hSurf.XData``, for example) do not equal to ``xdata``, ``ydata`` or ``zdata`` from ``getPlot2DData()``.
### Dataset manipulation
+ ``new_obj = removeDimensions(obj, varargin);``
Remove a few dimensions (parameters) from the dataset. Each attribute data will be sliced along these dimensions based on slice index. Since x, y, and z must be present in a rectilinear dataset, removing them is equivalent to reducing the length of these position vectors to 1.
+ ``new_obj = mergeDataset(obj, other_obj, varargin);``
Merge 2 datasets together. The dataset must have everything the same, except one parameter data. There are a few options for this function.