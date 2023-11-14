function [parameters_info, attributes_info, attributes_component, xyz] = loadLumDataset(lum_dataset)
% Load Lumerical exported MATLAB dataset into custom class object

% Do not pre-assign enough memory. Arrays are small anyways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problems not solved %%%%%%%%%%%%%%%%
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
% Sort the parameters?
% Zero parameters, zero attributes???
% parameters complex?
% Cannot be empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not allow no attribute
% Invalid (parameter, attribute) names, like "1", still passes
% data can have NaN? Inf?
% May need to label attribute scalar or vector somewhere??????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% values can be complex (parameters, xyz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First dimension: rectilinear dataset x*y*z
% Second dimension: if vector, 3
% Remaining dimensions: parameters

% Input dataset has to be a struct scalar
validateStructScalar(lum_dataset, "Input dataset must be a struct scalar!");

% Check field: Lumerical_dataset
validateFieldInStruct(lum_dataset, 'Lumerical_dataset', "Input dataset does not have the field 'Lumerical_dataset'!");
validateStructScalar(lum_dataset.Lumerical_dataset, "Field 'Lumerical_dataset' is not a struct scalar!");
% Check 'attribute' field
validateFieldInStruct(lum_dataset.Lumerical_dataset, 'attributes', "Field 'Lumerical_dataset' does not have the 'attributes' subfield!");
% check the contents in 'attributes' later
% Check 'parameters' field
validateFieldInStruct(lum_dataset.Lumerical_dataset, 'parameters', "Field 'Lumerical_dataset' does not have the 'parameters' subfield!")
% check the contents in 'parameters' later

% Decide between matrixdataset or rectilineardataset ('geometry' field)
if isfield(lum_dataset.Lumerical_dataset, 'geometry')
    dataset_type = 'rectilinear';
    if ~isequal(lum_dataset.Lumerical_dataset.geometry, "rectilinear")
        error("Wrong label in 'lum_dataset.geometry' for the rectilinear dataset!");
    end
else
    dataset_type = 'matrix';
end

% If it is rectilinear, load x,y,z and organize them
xyz = struct;
if isequal(dataset_type, 'rectilinear')
    % Check if x, y and z data exist in the dataset
    validateFieldInStruct(lum_dataset, 'x', "No x data in the rectilinear dataset!");
    validateFieldInStruct(lum_dataset, 'y', "No y data in the rectilinear dataset!");
    validateFieldInStruct(lum_dataset, 'z', "No z data in the rectilinear dataset!");

    % Load x,y,z to xyz and remove those field from the dataset
    % Check x,y,z data.
    for axis = 'xyz'
        if ~isnumeric(lum_dataset.(axis)) || isempty(lum_dataset.(axis))
            error(axis + " data must be a numeric vector!");
        end
        if ~isvector(lum_dataset.(axis)) % 2+ dimensional matrix?
            warning("Parameter " + axis + " is multi-dimensional! Stretched to one dimension!");
        end
        if any(imag(lum_dataset.(axis))) % has imaginary part?
            warning("Parameter " + axis + " is complex! Takes the real part and proceed.");
        end
    end
    xyz.x = real(lum_dataset.x(:)); % (vectorize) convert to column vector
    xyz.y = real(lum_dataset.y(:));
    xyz.z = real(lum_dataset.z(:));
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
% if rectilinear, the names can't be x, y, z
% check parameter values

parameters = lum_dataset.Lumerical_dataset.parameters;
if ~(iscell(parameters) && iscolumn(parameters))
    error("Field 'Lumerical_dataset.parameters' should be a cell column vector!");
end
parameters_info = cell(length(parameters), 3);
for i = 1:length(parameters)
    parameter = parameters{i};
    % Check struct
    if ~(isstruct(parameter) && isfield(parameter, 'variable') && isfield(parameter, 'name'))
        error("The interdependent parameter set " + i + " is not properly defined!");
    end

    % Initialize variables
    parameter_names = strings(1, length(parameter));
    parameter_values = cell(1, length(parameter));
    parameter_length = nan(1, length(parameter));

    % Retrieve the interdependent parameter set
    for j = 1:length(parameter)
        validateTextScalar(parameter(j).variable, "The interdependent parameter set " + i + " names cannot be resolved!");
        interdep_parameter_name = parameter(j).variable;
        parameter_names(j) = interdep_parameter_name; % char array also works
        validateTextScalar(parameter(j).name, "The interdependent parameter set " + i + " names cannot be resolved!");
        if ~isvarname(parameter(j).variable) % must be legal variable name
            error("The interdependent parameter set " + i + " names cannot be resolved!");
        end
        % Illegal characters in the names are converted to '_'
        if ~isequal(parameter(j).variable, parameter(j).name)
            warning("Parameter name '" + parameter(j).name + ...
                "' contains illegal characters. Converted to '" + parameter(j).variable + "'.");
        end
        validateFieldInStruct(lum_dataset, interdep_parameter_name, ...
            "Parameter field '" + interdep_parameter_name + "' data not found!");
        value = lum_dataset.(interdep_parameter_name);
        % parameter will always be non-empty N-by-1 column vector. Check that.
        % In the check, can relax the condition to "non-empty vectors".
        validateNonEmptyNumericVector(value, ... %% instead check iscolumn?
            "Parameter field '" + interdep_parameter_name + "' data is not a numeric vector!");
        value = value(:); % convert to column vector, if applicable
        % Remove this field from the dataset, prevent duplicate names
        lum_dataset = rmfield(lum_dataset, interdep_parameter_name);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If parameter not monotonic increasing or decreasing, issue a
        % warning
        if ~all(diff(value) >= 0) && ~all(diff(value) <= 0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plotting should not work if it is not monotonic
            warning("Data of parameter '" + interdep_parameter_name + "' is not monotonic! Plots cannot be made correctly!");
        end
        parameter_length(j) = length(value);
        parameter_values{j} = value; % each value could be different length

    end
    if ~isscalar(unique(parameter_length)) % lengths all the same?
        error("Interdependent parameters data do not have the same length!");
    end
    parameters_info{i, 1} = parameter_names;
    parameters_info{i, 2} = cell2mat(parameter_values); % values same length, combine
    parameters_info{i, 3} = parameter_length(1);
end

% Read attributes
% Attributes could be scalar attribute, or vector(x, y, z) attributes,
% depending on the type of the dataset
attributes = lum_dataset.Lumerical_dataset.attributes;
if ~(isstruct(attributes) && isfield(attributes, 'variable') && isfield(attributes, 'name'))
    error("The attributes are not properly defined!");
end

attributes_info = struct;
for i = 1:length(attributes)
    attribute = attributes(i);
    % Same here, illegal characters in the names are converted to '_'
    % Check names are text and legal variable name
    validateTextScalar(attribute.variable, "One or more attribute names cannot be resolved!");
    validateTextScalar(attribute.name, "One or more attribute names cannot be resolved!");
    if ~isvarname(attribute.variable)
        error("One or more attribute names cannot be resolved!");
    end
    if ~isequal(attribute.variable, attribute.name)
        warning("Attribute name '" + attribute.name + ...
            "' contains illegal characters. Converted to '" + attribute.variable + "'.");
    end
    % Check attribute data exists
    validateFieldInStruct(lum_dataset, attribute.variable, "Attribute field '" + attribute.variable + "' data not found!");
    attribute_value = lum_dataset.(attribute.variable);
    % Verify attribute value non-empty numeric
    if ~isnumeric(attribute_value) || isempty(attribute_value)
        error("Attribute field '" + attribute.variable + "' data must be numeric!");
    end
    % Check first dimension: should equal to multiplied x,y,z lengths
    if isequal(dataset_type, 'rectilinear')
        if size(attribute_value, 1) ~= prod(xyz.size)
            error("Unexpected size for attribute field '" + attribute.variable + "' data at 1st dimension!");
        end
    else % should assume dataset_type == 'matrix'
        if size(attribute_value, 1) ~= 1
            error("Unexpected size for attribute field '" + attribute.variable + "' data at 1st dimension!");
        end
    end
    % Check second dimension: scalar or vector
    if size(attribute_value, 2) == 1
        attributes_component.(attribute.variable) = NaN;
    elseif size(attribute_value, 2) == 3
        attributes_component.(attribute.variable) = 0; % default-magnitude
    else
        error("Unexpected size for attribute field '" + attribute.variable + "' data at 2nd dimension!");
    end
    % Check remaining dimensions, should agree with each parameter length
    for k = 1:size(parameters_info, 3)
        if size(attribute_value, k + 2) ~= parameters_info{k, 3}
            error("Unexpected size for attribute field '" + attribute.variable + "' data at dimension " + (k+2) + " !");
        end
    end
    % If ndims-2 < number of parameters, that means there is no extra
    % dimension(s) in the attribute data
    if ndims(attribute_value) > size(parameters_info, 1) + 2
        error("Too many dimensions for attribute field '" + attribute.variable + "' data!");
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

function validateNonEmptyNumericVector(input, errmsg)
% Validate input as non-empty vector (N-by-1 or 1-by-N)
if ~(isnumeric(input) && isvector(input) && ~isempty(input))
    throwAsCaller(MException('', errmsg));
end
end