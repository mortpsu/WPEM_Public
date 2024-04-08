function[out_mat]=ReshapeMatrix_DREM(long_mat)

    % define value columns
    val_col = size(long_mat,2);
    
    % define if long matrix array is 2 or 3 sets
    set_end = val_col - 1;
    
    % define two sets to search on 
    set1_col = 1;
    set2_col = set_end;
    
    % define set ids bus and economic region
    sets = unique(long_mat(:,set1_col:set2_col-1),'rows');
    
    % define unique list of sets of long matrix array
    list_set1 = unique(long_mat(:,set1_col)); % buses I of economic region R
    list_set2 = unique(long_mat(:,set2_col)); % hours D   
    
    % define length of sets of long matrix array
    len_set1 = length(list_set1);      % buses I of economic region R
    len_set2 = length(list_set2);      % hours D 

    % pre-allocate wide matrix array
    wide_mat = zeros(len_set1, len_set2);
    
    % loop to reshape long matrix arrary from gams 
    %       value col (i or r x d, 1) to (i or r x d) wide matrix arrary
    for i = 1:1:len_set1       % bus i
        for j = 1:1:len_set2   % hour d
            
            indices       = find((long_mat(:,set1_col) == list_set1(i)) & (long_mat(:,set2_col) == list_set2(j))); % find bus i and hour j 
            wide_mat(i,j) = long_mat(indices,val_col);
            
        end
    end
    
 out_mat = [sets, wide_mat];