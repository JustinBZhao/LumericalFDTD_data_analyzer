#LumericalDataset library v0.1.0
## General information
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
### Error while loading the dataset
Occasionally, you might see that the dataset was not converted, instead, an error was thrown. the ``createObject`` function performs extensive checking to ensure that the target being converted is indeed a dataset. Otherwise, it is carefully (hopefully) written to report specific errors regarding why the target cannot be recognized as a dataset. In almost all situations, the dataset directly exported from Lumerical should pass all the checkings. However, there are very rare situations where I don't think the dataset has any value to be further studied, even though it is legal in Lumerical (for example, a dataset with no attribute), so I will reject it.

These errors are here to help you, and to make sure that various assumtions (invariants) about the data inside the class holds. Don't be scared away.
## Perform operations on the converted dataset
Once you have successfully converted the dataset, the good things begin: there are many operations (class methods) you can perform on the dataset, including making 1D and 2D plots, merging datasets, etc. For example, 1D and 2D visualizations inside Lumerical visualizer can also be done here in MATLAB, but with prettier looking figures and better interactivity. You can also manipulate data in the dataset in several useful ways, which you CAN'T do in Lumerical.
## List of functions supported in this library
Below is the list of functions that have been implemented in this library, organized by category. However, I might update and add some new functionalities in the future versions.
### Data loading
+ [static] ``obj = LumericalDataset.createObject(lum_dataset)``
This is a static function. ``lum_dataset`` is the struct that is contained inside Lumerical's exported ``.mat`` file. This function will correctly resolve the type of the dataset (matrix dataset or rectilinear dataset) automatically. Any conversion that does not result in an error gives you a valid object.
+ [static] ``converted_obj = LumericalDataset.createObjectFromMat(mat_name, [*variable])``
This is a static function. Instead of passing the struct object as an input argument, you can directly specify the name of the ``.mat`` file in ``mat_name``. Additionally, you can also specify the list of variable names to be loaded in ``*variable``.
Note that in this case, ``converted_object`` is a struct containing each converted variable inside the original ``.mat`` file. Variables that are not Lumerical datasets are not converted. This is especially useful when the ``.mat`` file contains a few datasets, as well as some other non-dataset information.
+ ``obj = MatrixDataset(lum_dataset)`` or ``obj = RectilinearDataset(lum_dataset)`` (not recommended)
If you know which type of dataset you are converting (matrix dataset or rectilinear dataset), you can choose the correct syntax to manually perform the conversion. You **DO NOT** have to do it this way!
+ [static] ``new_obj = LumericalDataset.convertDatasetArray(obj_array, new_parameter_name, new_parameter_values)``
This is a static function. ``obj_array`` is an 1D array (or cell array) of unconverted or converted datasets. This function concatenates the datasets into one single dataset, ``new_obj``, if the parameters are exactly the same in each dataset and the attributes are the same except their values. After concatenation, a new parameter will be created, with ``new_parameter_name``, and ``new_parameter_values``. The resulting dataset will be ``new_obj``.
### Display information
+ ``obj.showInformation()``
This function will print out the basic information of the dataset, including the parameters and the attributes.
### Get raw data
+ ``result = obj.getParameterData(parameter_name)``
Retrieve the vector data corresponding to a specific parameter ``parameter_name``. For rectilinear datasets, you can also call ``x = obj.x``, for example, to retrieve x, y, or z data.
+ ``result = obj.getAttributeData(attribute_name)`` (not recommended)
Retrieve the data corresponding to a specific attribute ``attribute_name``. Use of this function is not recommended unless you understand the implementation detail of this library. Try ``getPlot1DData()``, ``getPlot2DData()`` or ``getPlot3DData()`` to retrieve attribute data corresponding to specific parameters.
### Set operations
+ ``obj.setParameterSliceIndex([**parameter_name=index])``
Set the slice index of one or more parameters, including x, y, or z in rectilinear datasets. Input should be in the form of "parameter\_name-index" pair(s).
The slice index for each parameter is 1 by default. The value persists throughout the lifetime of the dataset object.
+ ``obj.setParameterSliceValue([**parameter_name=value])``
Instead of the index, specify the value.
+ ``obj.setAttributeComponent(attribute_name, component)``
Set the component to plot for a vector attribute. ``component`` must be "x", "y", "z" or "magnitude". Cannot set this property for any scalar attribute (and will receive an error).
The attribute component for each vector attribute is set to "magnitude" by default. In addition, this property is stored in the dataset and will persists throughout the lifetime of the dataset object.
### Get data for making plots
+ ``[xdata, ydata] = obj.getPlot1DData(parameter_name, attribute_name)``
Retrieve the data for 1D plot. ``xdata`` is the parameter vector, ``ydata`` is the attribute data against that parameter, after slicing off all other dimensions. 
+ ``[xdata, ydata, zdata] = obj.getPlot2DData(parameter1_name, parameter2_name, attribute_name)``
Retrieve the data for 2D plot. ``xdata`` is the parameter\_1 vector, ``ydata`` is the parameter\_2 vector, and ``zdata`` is the attribute data, against these two parameters and after slicing off all other dimensions.
Note that for 2D plots, the MATLAB convention is that ``size(zdata) = [length(ydata), length(xdata)]``.
+ ``[x, y, z, data] = obj.getPlot3DData(parameter1_name, parameter2_name, parameter3_name, attribute_name)``
Retrieve the data for 3D plot. ``x``, ``y`` and ``z`` correspond to 3 parameters, and ``data`` is the attribute after slicing off all other dimensions. Currently, MATLAB does not have a way to image "3D" data. Therefore, this function does not have the corresponding plotting function.
Note that following MATLAB convention, ``size(data) = [length(y), length(x), length(z)]``.
+ ``interp_data = obj.getInterpolatedPlot2DData(parameter1_name, parameter2_name, attribute_name, Xq, Yq, method, extrapval, options)``
When retrieving the data for 2D plot, also perform an interpolation on to the new x (``Xq``) and new y (``Yq``). Optionally specify ``method``, ``extrapval`` and ``options``. Please see MATLAB document for ``interp2`` function for more details.
Since an interpolation is performed from ``Xq`` and ``Yq``, only the interpolated zdata ``interp_data`` is returned from the function.
### Make plots
+ ``hPlot = obj.plotData1D(parameter_name, attribute_name, scalar_operation, ax)``
Make a 1D plot with ``parameter_name`` against ``attribute_name``. Optionally, you can specify the scalar operation ("real" (default), "imag", "abs" or "angle") ``scalar_operation`` on the attribute data. You can also specify the axis ``ax`` to be plot on. The handle to the line object is returned as ``hPlot``.
+ ``[hSurf, hClb] = obj.plotData2D(parameter1_name, parameter2_name, attribute_name, scalar_operation, ax)``
Make a 2D plot with ``parameter1_name`` (x) and ``parameter2_name`` (y) against ``attribute_name``. Similarly, you can optionally provide ``scalar_operation`` and ``ax``. The handles to the surface object and the colorbar object are returned as ``hSurf`` and ``hClb``, respectively.
The 2D plot produced here is a surface plot in 2D view. However, ``hSurf.XData``, ``hSurf.YData`` and ``hSurf.ZData`` do not equal to ``xdata``, ``ydata`` or ``zdata`` from ``getPlot2DData()``[^1].
### Dataset manipulation
+ ``new_obj = obj.removeDimensions([*parameter_name]);``
Remove one or more parameters from the dataset specified by ``*parameter_names``. The dimensions in each attribute data corresponding to these parameters will be removed by slicing along these dimensions based on the slice index. x, y, and z positional vectors in a rectilinear dataset will be treated slightly differently, as they will be reduced to a singleton dimension instead of being completely removed [^2].
+ ``new_obj = obj.mergeDataset(other_obj, [**option=value]);``
Merge 2 datasets (``obj`` and ``other_obj``) together. Both datasets must have the same sets of parameters, except that one of the parameter data could be different. Can specify options in the form of option-value pairs.
"ParameterName": specify the name of the parameter to merge. Otherwise, it will be inferred from the content in the datasets.
"RemoveDuplicate": default is``false``. If ``true``, when merging the parameter data that is different between two datasets, Duplication values will be removed and the corresponding attributes data will be averaged to yield the new values. In this context, "duplication" is defined as values different within the tolerance level, specified by the "Tolerance" option.
"Sort": default is ``false``. If ``true``, the merge parameter data will also be sorted, along with the corresponding attributes data.
"Tolerance": default is [0, 0]. Specify as a two-element non-negative numeric vector in the form of [AbsoluteTolerance, RelativeTolerance]. The tolerance level is satisfied when either absolute or relative tolerance is satisfied.
**Note**: this function for rectilinear datasets is not yet implemented.

[^1]: This is because the surface plot sets the color of the cell edges, not the cell face area. Without additional treatment, the last row and last column of the attribute data will not be displayed in the surface plot.
[^2]: This is because rectilinear datasets require x, y and z positional vectors to be present. Even after they are removed, they will still show up in ``showInformation()`` outputs.
