function [parameters_info, attributes_info, attributes_component, xyz] = loadLumDataset(lum_dataset)
% Load Lumerical exported MATLAB dataset into custom class object

% Parameter and attribute data array cannot be empty.
% Missing required data will result in an error. However, more (useless)
% data will not trigger the error, and will simply be ignored.
% In a rectilinear dataset, non-positional vectors cannot be named x,y or
% z.
% If a parameter vector contains duplicate element or is not monotonic, no
% error/warning will be thrown. However, you cannot make 2D plot on that
% parameter.
% If a parameter or positional vectors (x,y,z) has complex data, it will be
% taken real part when converting.
% A parameter data does not need to be sorted.
% An empty matrix dataset is not accepted.
% A dataset without any attribute is not accepted.

