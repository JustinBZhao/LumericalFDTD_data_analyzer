classdef MatrixDataset < LumericalDataset
    % Subclass for matrix dataset

    properties
    end

    methods
        function obj = MatrixDataset(lum_dataset)
            % Parse dataset structure, determine dataset type
            dataset_type = LumericalDataset.parseDatasetStructure(lum_dataset);
            if dataset_type ~= "matrix"
                error("This is not a matrix dataset!");
            end

            % Parse parameter and attribute names
            parameters_info = LumericalDataset.parseParametersName(lum_dataset);
            attributes_names = LumericalDataset.parseAttributesName(lum_dataset);
            % No duplicate names allowed
            allnames = [parameters_info{:, 1}, string(attributes_names).'];
            if length(unique(allnames)) ~= length(allnames)
                error("Parameters and attributes have duplicate names! Rejected!");
            end

            % Parse parameter data
            obj.parameters = LumericalDataset.parseParametersData(lum_dataset, parameters_info);
            obj.num_parameters = size(obj.parameters, 1);
            obj.parameters_indexes = ones(obj.num_parameters, 1); % initialize to all 1

            % Parse attribute data
            [obj.attributes, obj.attributes_component] = LumericalDataset.parseAttributesData(lum_dataset, attributes_names, obj.parameters, dataset_type);
            obj.num_attributes = length(fieldnames(obj.attributes));
        end

        function showInformation(obj)
            % Print matrix dataset header
            fprintf("This dataset is a matrix dataset.\n");
            fprintf('\n');
            showInformation@LumericalDataset(obj);
            fprintf('\n');
            % Print parameters information
            LumericalDataset.printParametersInfo(obj.parameters(:, 1), [obj.parameters{:, 3}], obj.parameters_indexes);
        end

        function result = getParameterData(obj, parameter_name)
            % Get the data of a parameter
            para_loc = obj.iCheckAndFindParameter(parameter_name);
            result = obj.parameters{para_loc(1), 2}(:, para_loc(2));
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

            if size(para_value_list, 1) >= size(obj.parameters, 1)
                error("Should have at least one parameter left after removal!");
            end

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
            % Merge two datasets into one along one parameter direction
            % All parameter names should be the same, and every other
            % parameter data should be the same. After merging, that
            % parameter slice index will be reset to 1.
            % Options:
            % 'ParameterName': specify the parameter to merge. If not
            % specified, it will be deduced when applicable.
            % 'RemoveDuplicate': default is false. Boolean specifying
            % whether or not to remove the duplicate parameter points when
            % merging. Duplication is defined as values within the
            % tolerance specified by 'Tolerence'. Both the parameter
            % (including other interdependent parameters) data and the
            % corresponding attribute data for duplicate points will be
            % averaged to yield the new data. The new data will be place at
            % the location of the first instance of duplicate values in the
            % resulting array, whether sorted or unsorted.
            % 'Sort': default is false. Boolean specifying whether to sort
            % the resulting parameter (along with the attribute data) or
            % not. Currently only supports monotonic increase sort. The
            % original order is preserved only for EXACTLY equal elements,
            % if you choose not to remove duplicates.
            % 'Tolerance': default is [0, 0]. Specified as two-element
            % non-negative numeric array corresponding to
            % [AbsoluteTolerance, RelativeTolerance]. The tolerance is
            % satisfied when either absolute or relative tolerance is
            % satisfied.

            % Parse name-value pair
            p = createParser();
            p.parse(varargin{:});

            abs_tol = p.Results.Tolerance(1);
            rel_tol = p.Results.Tolerance(2);

            % Go through the dataset
            if class(other_obj) ~= "MatrixDataset"
                error("Matrix dataset should merge with another matrix dataset!");
            end
            % Same attribute fields
            attributes_names = fieldnames(obj.attributes);
            if ~isequal(attributes_names, fieldnames(other_obj.attributes))
                error("Unable to merge: two datasets do not have the same attribute sets!");
            end
            % Same type (scalar, vector) for each attribute
            % Shortcut: check attributes_component (both NaN or neither)
            for i1 = 1:length(attributes_names)
                name = attributes_names{i1};
                if (isnan(obj.attributes_component.(name)) && ~isnan(other_obj.attributes_component.(name))) || ...
                        (~isnan(obj.attributes_component.(name)) && isnan(other_obj.attributes_component.(name)))
                    error("Unable to merge: at least one attribute in two datasets do not have the same type (scalar or vector)!");
                end
            end
            % Same parameter fields
            if ~isequal(obj.parameters(:, 1), other_obj.parameters(:, 1))
                error("Unable to merge: two datasets do not have the same parameter sets!");
            end

            % Parameter data: at most one can be different
            % Here, equal within tolerance with AbsTol=1e-12, RelTol=1e-10
            param_data_diff = false(obj.num_parameters, 1);
            for i2 = 1:obj.num_parameters
                param_data_diff(i2) = ~LumericalDataset.isequalWithinTol( ...
                    obj.parameters{i2, 2}, other_obj.parameters{i2, 2});
            end

            if nnz(param_data_diff) == 0 % sizes all the same
                % Need to have user specified parameter name
                if p.Results.ParameterName == "" % not specified
                    error("Unable to merge: unable to deduce the parameter to merge!");
                end
                parameter_name = p.Results.ParameterName;
            elseif nnz(param_data_diff) == 1 % one size different
                % if parameter name specified, check if it is in there
                % otherwise, it cannot be interdependent set
                names_data_diff = obj.parameters{param_data_diff ~= 0, 1}; % find the name of different size
                if p.Results.ParameterName == "" % name not specified
                    if length(names_data_diff) > 1
                        error("Unable to merge: unable to deduce the parameter to merge!");
                    elseif length(names_data_diff) == 1
                        parameter_name = names_data_diff;
                    end
                else % name specified
                    parameter_name = p.Results.ParameterName;
                    if ~ismember(parameter_name, names_data_diff) % specified name does not match this interdependent set
                        error("Cannot merge the specified parameter name!");
                    end
                end
            else % more than one parameter data different
                error("Unable to merge: more than one parameter have different data!");
            end

            % Now we know the parameter to merge. Start manipulating data
            para_loc = obj.iCheckAndFindParameter(parameter_name);

            % Combine into new_obj (new_obj ~ [obj other_obj])
            new_obj = createNewObj(obj, other_obj, attributes_names, para_loc);

            interdep_param_data = new_obj.parameters{para_loc(1), 2};

            if p.Results.RemoveDuplicate
                interdep_param_data = removeDuplicate(new_obj, interdep_param_data, attributes_names, para_loc);
            end

            % Sort parameter and attributes if requested
            if p.Results.Sort
                [~, sort_idx] = sort(interdep_param_data(:, para_loc(2)), 1);
                interdep_param_data = interdep_param_data(sort_idx, :);
                for i = 1:length(attributes_names)
                    name = attributes_names{i};
                    slicing_sort = repmat({':'}, 1, 2+obj.num_attributes);
                    slicing_sort{para_loc(1) + 2} = sort_idx;
                    new_obj.attributes.(name) = new_obj.attributes.(name)(slicing_sort{:});
                end
            end

            % Assign parameter back to dataset
            new_obj.parameters{para_loc(1), 2} = interdep_param_data;
            % Update parameter data length
            new_obj.parameters{para_loc(1), 3} = size(interdep_param_data, 1);
            % Update parameter slice index
            new_obj.parameters_indexes(para_loc(1)) = 1;

            function p = createParser()
                % Parse name-value pair
                p = inputParser();
                p.PartialMatching = false;
                p.CaseSensitive = true;
                p.addParameter('Sort', false, ...
                    @(x) validateInput(x, [true, false], "Expected true or false."));
                p.addParameter('RemoveDuplicate', false, ...
                    @(x) validateInput(x, [true, false], "Expected true or false."));
                p.addParameter('Tolerance', [0, 0], @mustBeNonnegativeDual);
                p.addParameter('ParameterName', "", @mustBeTextScalar); % parameter name cannot be empty
            end

            function validateInput(input, valid_options, errmsg)
                for option = valid_options
                    if isequal(input, option)
                        return;
                    end
                end
                error(errmsg);
            end

            function mustBeNonnegativeDual(input) % non-negative two-element numeric array
                mustBeNonnegative(input);
                if numel(input) ~= 2
                    error("Value must be a two-element array.");
                end
            end

            function new_obj = createNewObj(obj, other_obj, attributes_names, para_loc)
                % Combine into new_obj (new_obj ~ [obj other_obj])
                new_obj = obj.copy();
                new_obj.parameters{para_loc(1), 2} = cat(1, obj.parameters{para_loc(1), 2}, other_obj.parameters{para_loc(1), 2});
                for i3 = 1:length(attributes_names)
                    name = attributes_names{i3};
                    new_obj.attributes.(name) = cat(para_loc(1) + 2, obj.attributes.(name), other_obj.attributes.(name));
                end
            end

            function interdep_param_data = removeDuplicate(new_obj, interdep_param_data, attributes_names, para_loc)
                % Remove duplicate for the combined parameter and
                % corresponding attribute values

                % Sort the array for easier duplication checking
                [param_data_sorted, param_data_sorted_order] = sort(interdep_param_data(:, para_loc(2)), 1);

                % Traverse through sorted list, find all duplicate sets
                duplicates_set = cell(0, 0); % empty cell
                tf_continue_current_set = false;
                for i4 = 1:length(param_data_sorted)-1
                    % Here, I am assuming that for 3 adjacent elements
                    % [a,b,c], if a==b and b==c within the tolerance, then
                    % a,b,c are all in the same duplicate set. In reality,
                    % a==c within the tolerance might not hold.
                    if LumericalDataset.isequalWithinTol( ...
                            param_data_sorted(i4), param_data_sorted(i4+1), abs_tol, rel_tol)
                        if tf_continue_current_set % continue in the current set
                            % Store the index in the original unsorted array
                            duplicates_set{end}(end+1) = param_data_sorted_order(i4+1);
                        else % start of new set
                            duplicates_set{end+1}(1) = param_data_sorted_order(i4); % "push_back"
                            duplicates_set{end}(2) = param_data_sorted_order(i4+1);
                        end
                        tf_continue_current_set = true;
                        % Where the next group begin?
                    else
                        tf_continue_current_set = false;
                    end
                end
                % For each duplicate set:
                data_to_keep = true(size(interdep_param_data, 1), 1);
                for i5 = 1:length(duplicates_set)
                    % Sort each duplicate set
                    duplicates_set{i5} = sort(duplicates_set{i5});
                    % Calculate and assign average of parameter and attributes
                    % data to the first index location
                    interdep_param_data(duplicates_set{i5}(1), :) ...
                        = mean(interdep_param_data(duplicates_set{i5}, :), 1);
                    for j = 1:length(attributes_names)
                        name = attributes_names{j};
                        % Acquire the slice of attribute data and reassign
                        slicing_first = repmat({':'}, 1, 2+new_obj.num_attributes);
                        slicing_first{para_loc(1) + 2} = duplicates_set{i5}(1);
                        slicing_all = repmat({':'}, 1, 2+new_obj.num_attributes);
                        slicing_all{para_loc(1) + 2} = duplicates_set{i5};
                        new_obj.attributes.(name)(slicing_first{:}) = ...
                            mean(new_obj.attributes.(name)(slicing_all{:}), para_loc(1) + 2);
                    end
                    % Mark other instances as 'to be deleted'
                    data_to_keep(duplicates_set{i5}(2:end)) = false;
                end
                % Remove all deleted instances
                interdep_param_data = interdep_param_data(data_to_keep, :);
                for i6 = 1:length(attributes_names)
                    name = attributes_names{i6};
                    slicing_keep = repmat({':'}, 1, 2+new_obj.num_attributes);
                    slicing_keep{para_loc(1) + 2} = data_to_keep;
                    new_obj.attributes.(name) = new_obj.attributes.(name)(slicing_keep{:});
                end
            end
        end

        function new_obj = plus(obj, other_obj)
            % MatrixDataset + MatrixDataset
            if isa(obj, "MatrixDataset") && isa(other_obj, "MatrixDataset")
                new_obj = plus@LumericalDataset(obj, other_obj);
                return;
            end
            % Dataset with a number
            if isa(obj, "MatrixDataset") && isnumeric(other_obj) && isscalar(other_obj)
                % MatrixDataset + num
                new_obj = obj.copy();
                num = other_obj;
            elseif isa(other_obj, "MatrixDataset") && isnumeric(obj) && isscalar(obj)
                % num + MatrixDataset
                new_obj = other_obj.copy();
                num = obj;
            else % Wrong inputs
                error("You can only add a matrix dataset with another matrix dataset or a number!");
            end
            if isnan(num) || isinf(num) % Warning on NaN or Inf
                warning("The number is NaN or Inf!");
            end
            new_obj.operateOnAllAttributesData(@(dataset) dataset + num);
        end
    end

    methods (Access = protected)
        function iSetParameterSlice(obj, mode_flag, varargin)
            % Set slice position of one or more parameters. Name-value
            % pairs consist of parameter name and index or value.
            % If mode_flag is "index", than accept indexes.
            % If mode_flag is "value", than accept values.
            % This is a helper function and should not be used directly by
            % the users.

            % Initialize inputParser, add all parameters
            p = inputParser();
            p.PartialMatching = false;
            p = obj.iAddAllParametersToParser(p, mode_flag);

            % Parse input name-value pairs
            try
                p.parse(varargin{:});
            catch ME
                ME.throw();
            end

            % Analyze indexes or values provided for the parameters If
            % multiple interdependent parameters are declared, they must
            % resolve to the same index. Otherwise, report an error.
            parsed_index_list = obj.iAnalyzeParsedParameters(p, mode_flag);
            obj.iUpdateParametersSliceIndex(parsed_index_list);
        end

        function iCheckAllParametersEquivalent(obj, other_obj, absTol, relTol)
            % Check parameters (excluding x, y and z) in both MatrixDataset
            % objects equivalent, which includes parameter names and
            % parameter values (within certain tolerance level)
            if ~isequal(obj.parameters(:, 1), other_obj.parameters(:, 1))
                ME = MException('', "Parameter lists in both datasets are not equivalent!");
                ME.throwAsCaller();
            end
            for iP = 1:obj.num_parameters
                if nargin == 2
                    isequivalent = LumericalDataset.isequalWithinTol(obj.parameters{iP, 2}, other_obj.parameters{iP, 2});
                elseif nargin == 3
                    isequivalent = LumericalDataset.isequalWithinTol(obj.parameters{iP, 2}, other_obj.parameters{iP, 2}, absTol);
                elseif nargin == 4
                    isequivalent = LumericalDataset.isequalWithinTol(obj.parameters{iP, 2}, other_obj.parameters{iP, 2}, absTol, relTol);
                end
                if ~isequivalent
                    ME = MException('', "Parameter values in both datasets are not equal within the tolerance limit!");
                    ME.throwAsCaller();
                end
            end
        end
    end
end