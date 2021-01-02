# https://stackoverflow.com/a/27930367/10844976
function print_rgb(r, g, b, t)
    print("\e[1m\e[38;2;$r;$g;$b;249m",t)
end
