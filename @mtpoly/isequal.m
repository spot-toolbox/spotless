function y=isequal(a,b)

if ~isa(b,'mtpoly'),
    y=0;
    return
end
y=isequal(size(a),size(b))&&isequal(gets(a),gets(b));
