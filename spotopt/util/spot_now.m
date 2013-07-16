function [time] = spot_now()
    time{1} = cputime;
    time{2} = clock;
end