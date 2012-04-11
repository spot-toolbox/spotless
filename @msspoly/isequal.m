function y=isequal(a,b)

if ~isa(b,'msspoly'),
    y=1;
    return
end
y=isequal(size(a),size(b))&&isequal(gets(a),gets(b));