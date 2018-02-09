% Generate all the comparators for a bitonic sorting network
% This is a n x 2 matrix such that the (i,j)-th coordinate represents the
% j-th input (first or second) input to the i-th comparator 
function [in_pairs] = bitonic_comparators(num_input)
    in_pairs = bitonic_sort(1:num_input);

% SORTING NETWORK WIRING RELATED CODE
% INPUTS:
%   input_range: the range of indices that we should be sorting
% OUTPUTS:
%   comparator_tuples = [in_idx(1), in_idx(2)]
function [comparator_tuples] = bitonic_sort(input_range)
    num_input = length(input_range);

    % base case
    if num_input == 1
        comparator_tuples = [];
        return
    end

    % otherwise, we split the thing into a power of 2 chunk and a rest chunk
    num_top_half = 2^floor(log2(num_input-1));
    top_half_idx = 1:num_top_half;
    bottom_half_idx = (num_top_half+1):num_input;

    top_tuples = bitonic_sort(input_range(top_half_idx));
    bottom_tuples = bitonic_sort(input_range(bottom_half_idx));             
    merged_tuples = bitonic_merge(input_range); 
    
    comparator_tuples = [top_tuples; bottom_tuples; merged_tuples];

function [comparator_tuples] = bitonic_merge(input_range)
    num_input = length(input_range);
    
    % this shouldn't happen, since we call bitonic_sort which catches 
    % num_input == 1 case
    if num_input == 1
        comparator_tuples = [];
        return
    end

    % we split the thing into a power of 2 chunk and a rest chunk
    num_top_half = 2^floor(log2(num_input-1));
    num_bottom_half = num_input - num_top_half;
    
    % initial crossover stage
    not_crossover_idx = input_range(1:num_top_half - num_bottom_half);
    crossover_input = input_range(num_top_half - num_bottom_half + 1: ...
                                  num_input);
    crossover_tuples = bitonic_crossover(crossover_input);
    
    if num_input == 2
        comparator_tuples = crossover_tuples;
        return
    end

    input_range = [not_crossover_idx crossover_input];
    top_half_idx = input_range(1:num_top_half);
    bottom_half_idx = input_range(num_top_half+1:num_input);

    % recursive merge stage
    top_half_tuples = recursive_merge(top_half_idx);
    bottom_half_tuples = recursive_merge(bottom_half_idx);

    comparator_tuples = [crossover_tuples; top_half_tuples; bottom_half_tuples];

function comparator_tuples = bitonic_crossover(input_range);
    num_input = length(input_range);
    if mod(num_input,2) ~= 0
        error('bitonic_crossover only handles even numbered input \n');
    end
    
    half = num_input/2;
    
    comparator_tuples = [];
    
    for i = 1:half
        comparator_tuples = [comparator_tuples;
                             input_range(i) input_range(num_input + 1 - i)];
    end
   
function comparator_tuples = recursive_merge(input_range);
    num_input = length(input_range);

    % base case
    if num_input == 1
        comparator_tuples = [];
        return
    end
   
    % we split the thing into a power of 2 chunk and a rest chunk
    num_top_half = 2^floor(log2(num_input-1));
    num_bottom_half = num_input - num_top_half;
    
    num_compares = num_bottom_half;
    
    compare_tuples = [];

    for i = 1:num_compares
        compare_tuples = [compare_tuples;
                          input_range(i) input_range(num_top_half + i)];
    end

    top_half_tuples = recursive_merge(input_range(1:num_top_half));
    bottom_half_tuples = recursive_merge(input_range(num_top_half+1:num_input));

    comparator_tuples = [compare_tuples; top_half_tuples; bottom_half_tuples];