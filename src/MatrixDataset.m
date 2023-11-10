classdef MatrixDataset < LumericalDataset
    % Subclass for matrix dataset

    properties
    end

    methods
        function obj = MatrixDataset(lum_dataset)
            % Call superclass constructor
            obj = obj@LumericalDataset(lum_dataset);
        end

        function showInformation(obj)
            % Print one additional line
            fprintf("This dataset is a matrix dataset.\n");
            showInformation@LumericalDataset(obj);
        end

        function obj = setParameterSliceIndex(obj, varargin)
            % Set slice index of one or more parameters (specified as
            % name-value pair)

            % Initialize inputParser, add regular parameters
            p = inputParser();
            p.PartialMatching = false;
            p = obj.iAddParametersToParser(p);

            % Parse input argument
            try
                p.parse(varargin{:});
            catch ME
                throwAsCaller(ME);
            end

            % Analyze and set values for regular parameters
            % If multiple interdependent parameters are declared, their
            % indexes should be the same. Otherwise, report an error.
            obj.iAnalyzeAndSetParsedParameter(p);
        end

        function [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name)
            % Get x and y data for 1D plot
            para_value_list = cell(1, 2); % 1D, xdata
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter_name);
            xdata = para_value_list{1, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            ydata = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            ydata = squeeze(ydata(:, :, para_slice_indexes{:}));
        end

        function [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name)
            % Get x,y and z data for 2D plot
            if obj.num_parameters < 2
                error("Can't call this method on a dataset that has less than 2 parameters!");
            end
            para_value_list = cell(2, 2); % 2D, xdata & ydata
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter1_name, parameter2_name);
            xdata = para_value_list{1, 1};
            ydata = para_value_list{2, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            zdata = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            zdata = squeeze(zdata(:, :, para_slice_indexes{:}));

            % Check for interdependent parameters and rearrange zdata based
            % on the order of parameters
            if (para_value_list{1, 2} == para_value_list{2, 2})
                error("Can't plot against two interdependent parameters!");
            elseif (para_value_list{1, 2} > para_value_list{2 ,2}) % flipped
                zdata = zdata.'; % flip
            end

            zdata = zdata.'; % swap x,y dimensions
        end

        function [x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name)
            % Get x,y,z and data for 3D plot
            if obj.num_parameters < 3
                error("Can't call this method on a dataset that has less than 3 parameters!");
            end
            para_value_list = cell(3, 2); % 3D, x,y,z
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter1_name, parameter2_name, parameter3_name);
            x = para_value_list{1, 1};
            y = para_value_list{2, 1};
            z = para_value_list{3, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            data = squeeze(data(:, :, para_slice_indexes{:}));

            % Check for interdependent parameters and rearrange data based
            % on the order of parameters
            if (para_value_list{1, 2} == para_value_list{2, 2} || ...
                    para_value_list{2, 2} == para_value_list{3, 2} || ...
                    para_value_list{1, 2} == para_value_list{3, 2})
                error("Can't plot against two interdependent parameters!");
            end

            [~, order] = sort(cell2mat(para_value_list(:, 2)));
            [~, rank] = sort(order); % sort twice to get params ranking
            data = permute(data, rank); % rearrange based on params ranking
            data = permute(data, [2, 1, 3]); % swap x,y dimensions
        end

        function new_obj = removeDimensions(obj, varargin)
            % Remove some dimensions (parameters) of the dataset
            para_value_list = cell(nargin - 1, 2); % nargin includes obj
            [~, para_value_list, para_remove_indexes] = obj.iGenerateParametersSliceIndexAndData(para_value_list, varargin{:});

            new_obj = obj.copy();

            % Adjust attributes
            attributes_fields = fieldnames(obj.attributes);

            for i = 1:obj.num_attributes
                field = attributes_fields{i};
                attribute_data = new_obj.attributes.(field);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%% This is dumb implementation of removing dimension
                sz = size(attribute_data); % get the size before trimming
                attribute_data = attribute_data(:, :, para_remove_indexes{:});

                sz(cell2mat(para_value_list(:, 2)) + 2) = [];
                attribute_data = reshape(attribute_data, sz);
                new_obj.attributes.(field) = attribute_data;
            end

            % Remove parameters
            new_obj.parameters(cell2mat(para_value_list(:, 2)), :) = [];
            % Update parameters_indexes
            new_obj.parameters_indexes(cell2mat(para_value_list(:, 2))) = [];
            % Update num_parameters
            new_obj.num_parameters = new_obj.num_parameters - (nargin - 1);
        end

        function new_obj = mergeDataset(obj, other_obj, varargin)
            % Merge two datasets into one
            % Other than the specified parameter name, every other
            % parameter has to be exactly the same.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Same class? matrixdataset
            % Check everything else is the same.
            % What if the parameter is interdependent?
            % What if the parameter has overlap between objs?
            % What if the parameter are not same monotonic direction?
            % What will be the slicing (axis component) after merging?
            % Parameter orders: does not need to be the same???
            p = inputParser();
            p.PartialMatching = false;
            p.CaseSensitive = true;
            p.addParameter('RemoveDuplicate', false, ...
                @(x) validateInputText(x, [true, false], "Expected true or false."));
            p.addParameter('Tolerance', 0, @mustBeNonnegativeScalar); % default tolerance is 0
            p.addParameter('ParameterName', "", @mustBeTextScalar); % parameter name cannot be empty
            p.parse(varargin{:});

            function validateInputText(input, valid_options, errmsg)
                for option = valid_options
                    if isequal(input, option)
                        return;
                    end
                end
                error(errmsg);
            end

            function mustBeNonnegativeScalar(input)
                mustBeNonnegative(input);
                if ~isscalar(input)
                    error("Value must be scalar.");
                end
            end

            % Go through the dataset
            % Same attribute fields
            attribute_names = fieldnames(obj.attributes);
            if ~isequal(attribute_names, fieldnames(other_obj.attributes))
                error("Unable to merge: two datasets do not have the same attribute sets!");
            end
% Same type (scalar, vector) for each attribute
            % Shortcut: check attributes_component (both NaN or neither)
            for i = 1:length(attribute_names)
                name = attribute_names{i};
                if (isnan(obj.attributes_component.(name)) && ~isnan(other_obj.attributes_component.(name))) || ...
                    (~isnan(obj.attributes_component.(name)) && isnan(other_obj.attributes_component.(name)))
                    error("Unable to merge: at least one attribute in two datasets do not have the same type (scalar or vector)!");
                end
            end
            % Same parameter fields
            if ~isequal(obj.parameters(:, 1), other_obj.parameters(:, 1))
                error("Unable to merge: two datasets do not have the same parameter sets!");
            end

            % Parameter sizes: at most one can be different
            diff = cell2mat(obj.parameters(:, 3)) - cell2mat(other_obj.parameters(:, 3));
            if nnz(diff) == 0 % sizes all the same
                % Need to have user specified parameter name
                if p.Results.ParameterName == "" % not specified
                error("Unable to merge: unable to deduce the parameter to merge!");
                end
                parameter_name = p.Results.ParameterName;
            elseif nnz(diff) == 1 % one size different
                % if parameter name specified, check if it is in there
                % otherwise, it cannot be interdependent set
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp("aha");
            end
            disp(diff ~= 0);


            % ONly the size of each attribute data matters (size not need
            % to be the same)
            % If parameter name is not specified, find it
            if p.Results.ParameterName == "" % not specified
                % Go through the dataset
            end

            para_loc = obj.iCheckAndFindParameter(parameter_name);

            new_obj = obj.copy();

            %%%%%%%%%%%%%% Interdep complication not considered
            % Reverse monotonic decrease
            if ~all(diff(obj.parameters{para_loc(1), 2}(:, 1), 1, 1) > 0)
                obj_store = obj.parameters{para_loc(1), 2}(end:-1:1, :);
                obj_wrong = true;
            else
                obj_store = obj.parameters{para_loc(1), 2};
                obj_wrong = false;
            end
            if ~all(diff(other_obj.parameters{para_loc(1), 2}(:, 1), 1, 1) > 0)
                other_obj_store = other_obj.parameters{para_loc(1), 2}(end:-1:1, :);
                other_obj_wrong = true;
            else
                other_obj_store = other_obj.parameters{para_loc(1), 2};
                other_obj_wrong = false;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Temperary
            % WHAT tolerance?
            if abs(obj_store(end, 1) - other_obj_store(1, 1)) < 1e-18
                if obj_store(1, 1) < 400e-9 || obj_store(end, 1) > 1500e-9
                    cancel = -1;
                    new_obj.parameters{para_loc(1), 2} = [obj_store(1:end-1, :); other_obj_store];
                else
                    cancel = 1;
                    new_obj.parameters{para_loc(1), 2} = [obj_store; other_obj_store(2:end, :)];
                end
                new_obj.parameters{para_loc(1), 3} = obj.parameters{para_loc(1), 3} + other_obj.parameters{para_loc(1), 3} - 1;
            else
                cancel = 0;
                new_obj.parameters{para_loc(1), 2} = [obj_store; other_obj_store];
                new_obj.parameters{para_loc(1), 3} = obj.parameters{para_loc(1), 3} + other_obj.parameters{para_loc(1), 3};
            end

            % Adjust attributes
            attributes_fields = fieldnames(obj.attributes);

            % Also reverse monotonic decrease and THEN trim
            for i = 1:obj.num_attributes
                field = attributes_fields{i};
                if obj_wrong
                    obj_attri = flip(obj.attributes.(field), para_loc(1) + 2);
                else
                    obj_attri = obj.attributes.(field);
                end
                if other_obj_wrong
                    other_obj_attri = flip(other_obj.attributes.(field), para_loc(1) + 2);
                else
                    other_obj_attri = other_obj.attributes.(field);
                end

                if cancel == -1
                    % Remove last index along parameter dim
                    dim = para_loc(1) + 2;
                    idx = true(1, size(obj_attri, dim));
                    idx(end) = false;
                    subs = repmat({':'}, 1, ndims(obj_attri));
                    subs{dim} = idx;
                    obj_attri = obj_attri(subs{:});
                elseif cancel == 1
                    % Remove first index along parameter dim
                    dim = para_loc(1) + 2;
                    idx = true(1, size(other_obj_attri, dim));
                    idx(1) = false;
                    subs = repmat({':'}, 1, ndims(other_obj_attri));
                    subs{dim} = idx;
                    other_obj_attri = other_obj_attri(subs{:});
                end
                new_obj.attributes.(field) = cat(para_loc(1) + 2, obj_attri, other_obj_attri);
            end
        end
    end
end