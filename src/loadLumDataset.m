function [parameters_info, attributes_info, attributes_component, xyz] = loadLumDataset(lum_dataset)
% Do not pre-assign enough memory. Arrays are small anyways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check contents in Lumerical_dataset.attributes (also if it includes
% 'geometry')
% Check if it is actually a dataset (check every single aspect of it)?
% Did not check if there are more things in the dataset
% Check duplication???
% Check x, y, z?????????
% Only supports scalar or vector attributes (tensor?)
% In matrixdataset, you can have a parameter called 'x', etc.
% In rectilineardataset, if you name a parameter 'x', it will be
% transformed to 'x_2'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First dimension: rectilinear dataset x*y*z
% Second dimension: if vector, 3
% Remaining dimensions: parameters

% Input dataset has to be a scalar
validateStructScalar(lum_dataset, "Input dataset must be a struct!");

% Check field: Lumerical_dataset
validateFieldInStruct(lum_dataset, 'Lumerical_dataset', "Input dataset does not have the field 'Lumerical_dataset'!");
validateStructScalar(lum_dataset.Lumerical_dataset, "Field 'Lumerical_dataset' is not a struct!");
% Check 'attribute' field
validateFieldInStruct(lum_dataset.Lumerical_dataset, 'attributes', "Field 'Lumerical_dataset' does not have the 'attributes' subfield!");
% check the contents in 'attributes' later
% Check 'parameters' field
validateFieldInStruct(lum_dataset.Lumerical_dataset, 'parameters', "Field 'Lumerical_dataset' does not have the 'parameters' subfield!")
% check the contents in 'parameters' later

% Decide between matrixdataset or rectilineardataset ('geometry' field)
if isfield(lum_dataset.Lumerical_dataset, 'geometry')
    dataset_type = 'rectilinear';
    disp("This data set is a rectilinear dataset!");
    if lum_dataset.Lumerical_dataset.geometry ~= "rectilinear"
        error("Wrong label in 'lum_dataset.geometry' for the rectilinear dataset!");
    end
else
    dataset_type = 'matrix';
    disp("This dataset is a matrix dataset!");
end

% If it is rectilinear, load x,y,z and organize them
xyz = struct;
if isequal(dataset_type, 'rectilinear')
    % Check if x, y and z data exist in the dataset
    validateFieldInStruct(lum_dataset, 'x', "No x field in the rectilinear dataset!");
    validateFieldInStruct(lum_dataset, 'y', "No y field in the rectilinear dataset!");
    validateFieldInStruct(lum_dataset, 'z', "No z field in the rectilinear dataset!");

    % Load x,y,z to xyz and remove those field from the dataset
    xyz.x = lum_dataset.x;
    xyz.y = lum_dataset.y;
    xyz.z = lum_dataset.z;
    xyz.size = [length(lum_dataset.x), length(lum_dataset.y), length(lum_dataset.z)];
    % Remove x,y,z field from the dataset. This ensures that if other
    % parameters have these names, an error will be issued when we try to
    % look for them in the dataset
    lum_dataset = rmfield(lum_dataset, 'x');
    lum_dataset = rmfield(lum_dataset, 'y');
    lum_dataset = rmfield(lum_dataset, 'z');
end

% Read parameters
% Load all parameters names and organize them
% Will group interdependent parameters together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if rectilinear, the names cant be x, y, z
% check parameter values

parameters = lum_dataset.Lumerical_dataset.parameters;
if ~iscell(parameters) || ~iscolumn(parameters)
    error("Field 'Lumerical_dataset.parameters' does not have the correct format!");
end
parameters_info = cell(length(parameters), 3);
for i = 1:length(parameters)
    parameter = parameters{i};
    % Check struct
    if ~isstruct(parameter) || ~isfield(parameter, 'variable') || ~isfield(parameter, 'name')
        error("The " + i + "th interdependent parameter set is not properly defined!");
    end

    % Initialize variables
    parameter_names = strings(1, length(parameter));
    parameter_values = [];
    parameter_length = nan(1, length(parameter));

    % Retrieve the interdependent parameter set
    for j = 1:length(parameter)
        validateTextScalar(parameter(j).variable, "The " + i + "th interdependent parameter set is not properly defined!");
        validateTextScalar(parameter(j).name, "The " + i + "th interdependent parameter set is not properly defined!");
        if ~isvarname(parameter(j).variable)
            error("The " + i + "th interdependent parameter set is not properly defined!");
        end
        % Illegal characters in the names are converted to '_'
        if ~isequal(parameter(j).variable, parameter(j).name)
            warning("Parameter name '" + parameter(j).name + ...
                "' contains illegal characters. Converted to '" + parameter(j).variable + "'.");
        end
        interdep_parameter_name = parameter(j).variable;
        parameter_names(j) = interdep_parameter_name;
        validateFieldInStruct(lum_dataset, interdep_parameter_name, ...
            "Parameter field '" + interdep_parameter_name + "' data not found!");
        value = lum_dataset.(interdep_parameter_name);
        if ~isnumeric(value) || ~isvector(value) || isempty(value)
            error("Parameter field '" + interdep_parameter_name + "' data is not a numeric vector!");
        end
        % Remove this field from the dataset, prevent duplicate names
        lum_dataset = rmfield(lum_dataset, interdep_parameter_name);
        % If parameter not monotonic increasing or decreasing, issue a
        % warning
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check parameter value data type
        if ~all(diff(value) >= 0) && ~all(diff(value) <= 0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plotting should not work if it is not monotonic
            warning("Data of parameter '" + interdep_parameter_name + "' is not monotonic!");
        end
        parameter_length(j) = length(value);
        parameter_values(:, j) = value;

    end
    if ~isscalar(unique(parameter_length))
        error("Interdependent parameters not the same length!");
    end
    parameters_info{i, 1} = parameter_names;
    parameters_info{i, 2} = parameter_values;
    parameters_info{i, 3} = parameter_length(1);
end

% Read attributes
% Attributes could be scalar attribute, or vector(x, y, z) attributes,
% depending on the type of the dataset
attributes = lum_dataset.Lumerical_dataset.attributes;
if ~isstruct(attributes) || ~isfield(attributes, 'variable') || ~isfield(attributes, 'name')
    error("The attributes are not properly defined!");
end

attributes_info = struct;
for i = 1:length(attributes)
    attribute = attributes(i);
    % attribute.variable
    % attribute.name
    % Same here, illegal characters converted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check these fields are both text
    validateTextScalar(attribute.variable, "The " + i + "th attribute is not properly defined!");
    validateTextScalar(attribute.name, "The " + i + "th attribute is not properly defined!");
    if ~isvarname(parameter(j).variable)
        error("The " + i + "th attribute is not properly defined!");
    end
    if ~isequal(attribute.variable, attribute.name)
                    warning("Attribute name '" + attribute.name + ...
                "' contains illegal characters. Converted to '" + attribute.variable + "'.");
    end
    validateFieldInStruct(lum_dataset, attribute.variable, "Attribute field '" + attribute.variable + "' not found!");
    attribute_value = lum_dataset.(attribute.variable);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check attribute value data type
    % Check first dimension: should correspond to dataset_type
    if isequal(dataset_type, 'rectilinear')
        if size(attribute_value, 1) ~= prod(xyz.size)
            error("Unexpected size of attribute at 1st dimension (rectilinear)!");
        end
    else % should assume dataset_type == 'matrix'
        if size(attribute_value, 1) ~= 1
            error("Unexpected size of attribute at 1st dimension (matrix)!");
        end
    end

    % Check second dimension: scalar or vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % May need to label it somewhere??????
    if size(attribute_value, 2) == 1
        disp("The attribute '" + attribute.variable + "' is scalar! ");
        attributes_component.(attribute.variable) = NaN;
    elseif size(attribute_value, 2) == 3
        disp("The attribute '" + attribute.variable + "' is vector! ");
        attributes_component.(attribute.variable) = 0; % default-magnitude
    else
        error("Unexpected size of attribute at 2nd dimension!");
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Should also check that there is no extra non-singular dimension
    % Check remaining dimensions, should correspond to each parameter
    for k = 1:size(parameters_info, 3)
        if size(attribute_value, k + 2) ~= parameters_info{k, 3}
            error("Unexpected size of attribute at dimension " + k + " !");
        end
    end

    attributes_info.(attribute.variable) = attribute_value;
end
end

%% Private functions that should be eventually merged with the class file
function validateTextScalar(input, errmsg)
% Validate input as text scalar and throw if not
if ~(ischar(input) && isrow(input)) && ~isStringScalar(input)
    throwAsCaller(MException('', errmsg));
end
end

function validateStructScalar(input, errmsg)
% Validate input as struct scalar and throw if not
if ~isstruct(input) || ~isscalar(input)
    throwAsCaller(MException('', errmsg));
end
end

function validateFieldInStruct(struct_in, field_in, errmsg)
if ~isfield(struct_in, field_in)
    throwAsCaller(MException('', errmsg));
end
end