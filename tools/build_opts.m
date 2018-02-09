function opts = build_opts(opts_def, opts_arg)
% input : opts_def must be a structure
%         opts_arg must be a structure
% ouput : opts structure whose field names are the ones of opts_def
%         value of each fieldname is either corresponding value in opts_arg
%         or default value of opts_def if opts_arg has no corresponding
%         field name

    opts = struct;
    opts_names = fieldnames(opts_def);
    for i = 1:length(opts_names)
        if isfield(opts_arg,opts_names{i})
            value = getfield(opts_arg,opts_names{i});
            opts = setfield(opts,opts_names{i},value);
        else
            value = getfield(opts_def,opts_names{i});
            opts = setfield(opts,opts_names{i},value);
        end
    end

end