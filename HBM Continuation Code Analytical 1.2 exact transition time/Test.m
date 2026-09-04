flag = [2,2,2,2,1,1,1,2,2,2,2,-1,-1,-1,2,2,0,0]';
diffs = [diff(flag); flag(1) - flag(end)]  
trans_idx = find(diffs ~= 0);
