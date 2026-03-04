function out = make_python_bytes(in)

if length(in) <= 1
    out = py.bytes(uint8(in));
else
    s = sprintf('%d ',uint8(in));
    x = py.numpy.fromstring(s, py.numpy.int8, int8(-1), char(' '));
    out = py.bytes(x);
end

