function [cel,wel] = spot_etime(t2,t1)
    cel = t2{1} - t1{1};
    wel = etime(t2{2},t1{2});
end