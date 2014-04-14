% Test spot_mex_msspoly_check_canonical.
% TODO(mmt): Write tests for the above function.
% Test spot_mex_msspoly_make_canonical_combine_coeffs.
key = [ 1 1 1
        1 1 1
        2 1 1
        2 1 1
        3 3 3];
value = [ 1 2 -100 -200 1000 ].';

[k, kO, vO] = spot_mex_msspoly_make_canonical_combine_coeffs(key, value);

if k ~= 3, error('Wrong number of keys.'); end
if ~all(all(kO(1:3,:) == [1 1 1; 2 1 1; 3 3 3])), error('Wrong keys.'); end
if ~all(vO(1:3) == [3; -300; 1000]), error('Wrong values.'); end

key = [ 1 1 1
        2 1 1
        2 1 1
        3 3 3
        3 3 3];
value = [ 1000i 1-100i 2-200i 1i-200 2i-100 ].';
[k, kO, vO] = spot_mex_msspoly_make_canonical_combine_coeffs(key, value);

if k ~= 3, error('Wrong number of keys.'); end
if ~all(all(kO(1:3,:) == [1 1 1; 2 1 1; 3 3 3])), error('Wrong keys.'); end
if ~all(vO(1:3) == [1000i; 3-300i; 3i-300]), error('Wrong values.'); end

% Test spot_mex_msspoly_make_canonical_combine_powers.
powIn = [ 1 2 0
          4 5 0
          7 0 0];
varIn = [ 22 23 24
          22 22 0
          24 27 0 ];
[pow, var] = spot_mex_msspoly_make_canonical_combine_powers(powIn, ...
                                                  varIn);
if ~isequal(pow, [1 2 0; 9 0 0; 7 0 0]), error('Wrong powers.'); end
if ~isequal(var, [22 23 0; 22 0 0; 24 0 0]), error('Wrong variables.'); end