function y=basis_to_index(vector)

    y=vector*(2.^(length(vector)-1:-1:0))';

end