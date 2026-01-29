skip = [1 2; 2 1; 3 1; 3 2];
for i = 1:3
    for j = 1:3
        if ismember([i j], skip, 'rows')
            continue;
        end
        disp([i,j]);
    end
end
%%
