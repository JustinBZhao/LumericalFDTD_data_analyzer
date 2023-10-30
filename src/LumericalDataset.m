classdef (Abstract) LumericalDataset < matlab.mixin.Copyable
    % Lumerical dataset base class
    % When you want to make a copy of an object, call the copy() function
    % instead of using the assignment operator '='.

    properties (SetAccess = protected)
        attributes
        parameters
        attributes_component % select the current attribute component (x, y, z)
        % NaN-scalar 0-magnitude 1-x 2-y 3-z
        num_attributes
        parameters_indexes % selected the current parameter slice indexes
        num_parameters
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make only available in base class?
        function obj = createObject(lum_dataset)
            % Initialize the dataset object and choose the correct subclass
            % to invoke based on the type of the dataset

            %%%%%%%%%%%%%%%%%%%%%%%%
            % Do you check the input NOW???
            % This is problematic because it assumes the input is a struct
            if isfield(lum_dataset.Lumerical_dataset, 'geometry')
                obj = RectilinearDataset(lum_dataset);
            else
                obj = MatrixDataset(lum_dataset);
            end
        end
    end

    methods
        function obj = LumericalDataset(lum_dataset)
            % Constructor
            % Create a MATLAB data structure of a matrix dataset from the
            %   Lumerical exported MATLAB data

            %%%%%%%%%%%%%%%%%
            % Check the input is actually a lumerical dataset?
            % Duplicate parameter and attribute names?
            % Complex values?
            % Design assigning xyz

            % First load parameters
            [obj.parameters, obj.attributes, obj.attributes_component, ~] = loadLumDataset(lum_dataset);
            obj.num_attributes = length(fieldnames(obj.attributes));
            % If scalar, give NaN. If vector, assign 1.

            obj.num_parameters = size(obj.parameters, 1);
            % Initialize to all 1 (first value for each parameter)
            obj.parameters_indexes = ones(obj.num_parameters, 1);
            % obj.num_free_parameters = obj.num_parameters;
        end

        function showInformation(obj)
            % Display the information for an overview of this dataset
            % Print attributes information
            fprintf("%d attribute(s):\n", obj.num_attributes);
            attributes_fields = fieldnames(obj.attributes);
            max_str_length = max(cellfun(@length, attributes_fields)); % max name length
            for i = 1:obj.num_attributes
                field = attributes_fields{i};
                if isnan(obj.attributes_component.(field)) % scalar or vector
                    type = "scalar";
                else
                    type = "vector";
                end
                fprintf("%-*s | %s\n", max_str_length, field, type); % attribute names and types
            end
            fprintf("%d parameter(s):\n", obj.num_parameters);
            % Print parameters information
            params = obj.parameters(:, 1);
            max_digits = max(cellfun(@(x) length(num2str(x)), obj.parameters(:, 3))); % max number of points length
            for i = 1:obj.num_parameters
                fprintf("%*d point(s) | ", max_digits, obj.parameters{i, 3}); % number of points
                for param = params{i} % independent parameter names
                    fprintf("%s ", param);
                end
                fprintf('\n');
            end
        end

        function result = getParameterData(obj, parameter_name) % non-virtual
            % Get the data of a parameter
            para_loc = obj.iCheckAndFindParameter(parameter_name);
            result = obj.parameters{para_loc(1), 2}(:, para_loc(2));
        end

        function result = getAttributeData(obj, attribute_name) % non-virtual
            % Get the data for an attribute
            obj.iCheckAttributeExist(attribute_name);
            result = obj.attributes.(attribute_name);
        end

        function setAttributeComponent(obj, attribute_name, component) % non-virtual
            % set the component (x,y,z,magnitude) for an attribute
            obj.iCheckAttributeExist(attribute_name);
            if isnan(obj.attributes_component.(attribute_name))
                error("Scalar attribute cannot be modified!");
            end

            LumericalDataset.validateTextScalar(component, "Component name must be text!");
            switch component
                case "x"
                    obj.attributes_component.(attribute_name) = 1;
                case "y"
                    obj.attributes_component.(attribute_name) = 2;
                case "z"
                    obj.attributes_component.(attribute_name) = 3;
                case "magnitude"
                    obj.attributes_component.(attribute_name) = 0;
                otherwise
                    error("Invalid component name!");
            end
        end

        function [hAx, hPlot] = plotData1D(obj, parameter_name, attribute_name)  % non-virtual
            % This function makes a 1D plot based on the parameter name and
            % the attribute name.

            % Calls the respective polymorphic function for each object
            [xdata, ydata] = obj.getPlot1DData(parameter_name, attribute_name);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Problem (also for 2D):
            % complex numbers? real?
            % Unit conversion and scaling?
            figure;
            hPlot = plot(xdata, real(ydata));
            xlabel(parameter_name, 'Interpreter', 'none');
            ylabel(attribute_name, 'Interpreter', 'none');
            hAx = gca;
            set(hAx, 'FontSize', 16);
        end

        function [hAx, hSurf, hClbr] = plotData2D(obj, parameter1_name, parameter2_name, attribute_name)  % non-virtual
            % This function makes a 2D plot based on the parameter names
            % and the attribute name.

            % Calls the respective polymorphic function for each object
            [xdata, ydata, zdata] = obj.getPlot2DData(parameter1_name, parameter2_name, attribute_name);

            % Adjust data to make true "2D plot"
            % If length(xdata) == 1 or length(ydata) == 1, this program
            % will not report an error but will have nothing plotted
            xdata = ([xdata(1); xdata(:)] + [xdata(:); xdata(end)])/2;
            ydata = ([ydata(1); ydata(:)] + [ydata(:); ydata(end)])/2;
            zdata_new = zeros(length(ydata), length(xdata));
            zdata_new(1:end-1, 1:end-1) = zdata;
            zdata = zdata_new;

            figure;
            hSurf = surface(xdata*1e9, ydata*1e9, real(zdata));
            xlim([min(xdata*1e9), max(xdata*1e9)]);
            ylim([min(ydata*1e9), max(ydata*1e9)]);
            xlabel(parameter1_name + " (nm)", 'Interpreter', 'none');
            ylabel(parameter2_name + " (nm)", 'Interpreter', 'none');
            colormap jet;
            hClbr = colorbar;
            set(gca, 'FontSize', 16);
            hSurf.EdgeColor = 'none';
            hAx = gca;
            set(hAx, 'Layer', 'top');
            box on;
        end
    end

    methods (Abstract)
        obj = setParameterSliceIndex(obj, varargin);
        [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name);
        [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name);
        [x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name);
        new_obj = removeDimensions(obj, varargin);
    end

    methods (Access = protected)
        function p = iAddParametersToParser(obj, p)
            % Add all parameters and their range validations to the parser
            for i = 1:size(obj.parameters, 1)
                for para = obj.parameters{i, 1}
                    p.addParameter(para, NaN, @(x) LumericalDataset.validateIndex(x, obj.parameters{i, 3}));
                end
            end
        end

        function iAnalyzeAndSetParsedParameter(obj, p)
            % You can supply the same parameter multiple times and it will
            % keep the last occurance of it
            for i = 1:size(obj.parameters, 1)
                interdep_indexes = [];
                for para = obj.parameters{i, 1}
                    interdep_indexes(end+1) = p.Results.(para);  % retrieve
                end
                unique_index = unique(interdep_indexes(~isnan(interdep_indexes)));
                if numel(unique_index) > 1
                    throwAsCaller(MException('', "Interdependent parameters for '" + para + "' should select the same index!"));
                elseif numel(unique_index) == 1
                    obj.parameters_indexes(i) = unique_index;
                end
            end
        end

        function [para_slice_indexes, para_value_list, para_keep_indexes] = ...
                iGenerateParametersSliceIndexAndData(obj, para_value_list, varargin)
            % Generate index list that slices off passed in parameters
            % (para_slice_indexes), or index list that only keeps passed in
            % parameters (para_keep_indexes), and update the information of
            % the parameters passed in (para_value_list).
            % Certain elements in varargin can be empty, meaning that it
            % has already been parsed in the last step. In this case, the
            % corresponding element in para_value_list should already have
            % a value.
            % ---para_value_list---
            % First column: actual data for the parameter
            % Second column: the number of the parameter (1:num_parameters)
            try
                para_slice_indexes = num2cell(obj.parameters_indexes);
                for i = 1:length(varargin)
                    if isempty(varargin{i}) % already parsed
                        continue;
                    end
                    para_loc = obj.iCheckAndFindParameter(varargin{i});
                    para_slice_indexes{para_loc(1)} = ':';
                    % Extract found parameter data
                    para_value_list{i, 1} = obj.parameters{para_loc(1), 2}(:, para_loc(2));
                    para_value_list{i, 2} = para_loc(1);
                end

                para_keep_indexes = num2cell(obj.parameters_indexes);
                for i = 1:length(para_slice_indexes)
                    if isnumeric(para_slice_indexes{i})
                        para_keep_indexes{i} = ':';
                    end
                end
            catch ME
                throwAsCaller(ME)
            end
        end

        function location = iCheckAndFindParameter(obj, parameter_name)
            % Check if a parameter name exists and return location
            % Location: [row, column (within interdep)]
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                throwAsCaller(ME);
            end

            tf_found = false;
            location = zeros(1, 2);
            parameter_names = obj.parameters(:, 1);
            for i = 1:length(parameter_names)
                if ismember(parameter_name, parameter_names{i})
                    location(1) = i;
                    location(2) = find(parameter_names{i} == parameter_name);
                    tf_found = true;
                    break;
                end
            end
            if ~tf_found
                throwAsCaller(MException('', "Parameter name '" + parameter_name + "' not found!"));
            end
        end

        function iCheckParameterExist(obj, parameter_name)
            % Check if a parameter name exists
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                throwAsCaller(ME);
            end

            tf_found = false;
            parameter_names = obj.parameters(:, 1);
            for i = 1:length(parameter_names)
                if ismember(parameter_name, parameter_names{i})
                    tf_found = true;
                    break;
                end
            end
            if ~tf_found
                throwAsCaller(MException('', "Parameter name '" + parameter_name + "' not found!"));
            end
        end

        function iCheckAttributeExist(obj, attribute_name)
            % Check if an attribute name exists
            try
                LumericalDataset.validateTextScalar(attribute_name, "Attribute name must be text scalar!");
                LumericalDataset.validateFieldInStruct(obj.attributes, attribute_name, "This attribute '" + attribute_name + "' is not found!");
            catch ME
                throwAsCaller(ME);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deprecated, but this one takes in the value instead of the
        % indexes
        function parameters_value_list = parseParameterList(obj, parameters_value_list, arglist)
            % Check varargin (should be parameter name-value pair)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % In this version, only one interdependent parameter is allowed
            % to be declared. Further, pass in values rather than indexes.
            n = length(arglist);
            if mod(n, 2)
                error('Parameter name-value not in pairs!');
            end
            for i = 1:2:n
                parameter = arglist{i};
                value = arglist{i+1};
                para_loc = obj.iCheckAndFindParameter(parameter);
                if ~isnumeric(value) || ~isscalar(value)
                    error("The value for this parameter '" + parameter + "' is not a number!");
                end

                % Check if the parameter (or other interdep) has already been defined
                if ~isnan(parameters_value_list(para_loc(1), 1))
                    disp(parameter);
                    error('Repeat definition of this parameter!');
                end

                % Find the index of the closest value
                [~, value_index] = min(abs(obj.parameters{para_loc(1), 2}(:, para_loc(2)) - value));
                parameters_value_list(para_loc(1), 1) = para_loc(2);
                parameters_value_list(para_loc(1), 2) = value_index;
            end

            % Any not specified parameters are default to the first value
            parameters_value_list(isnan(parameters_value_list(:, 1)), :) = 1;
        end
    end

    methods (Static, Access = protected)
        % Helper functions to be shared with derived classes
        function validateTextScalar(input, errmsg)
            % Validate input as text scalar and throw if not
            if ~(ischar(input) && isrow(input)) && ~isStringScalar(input)
                throwAsCaller(MException('', errmsg));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Not used
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

        function validateIndex(idx, max_val)
            if ~isnumeric(idx) || floor(idx) ~= idx % test integer
                throwAsCaller(MException('', "Must be an integer!"));
            end
            if idx < 1 || idx > max_val
                throwAsCaller(MException('', "Must be between 1 and " + max_val + "!"));
            end
        end

        function result = sliceThroughAttributeVectorDim(attribute_data, component_index)
            % Slice off attribute data matrix based on the selected vector
            % component (NaN-scalar, 0-magnitude, 1-x, 2-y, 3-z)
            % on the second dimension
            if isnan(component_index) % NaN-scalar
                result = attribute_data;
            elseif component_index == 0 % 0-magnitude
                result = vecnorm(attribute_data, 2, 2);
            else % 1, 2, 3-x, y, z
                idx = repmat({':'}, 1, ndims(attribute_data)); % ugly
                idx{2} = component_index;
                result = attribute_data(idx{:});
            end
        end
    end
end