function reord_c_c = F_reorder(c_c)
    %UNTITLED4 Summary of this function goes here
    %   Detailed explanation goes here
    u_cc = unique(c_c);
    c_c2 = c_c;
    c_c3 = c_c;
    classe = 0;
    for i=1:length(u_cc)
        classe = classe + 1;
        c = c_c2(1);
        ind = find(c_c3 == c);
        for  j= 1:length(ind)
            c_c(ind(j)) = classe;
        end
        c_c2(c_c2 == c) = [];
    end
    reord_c_c = c_c;
end

