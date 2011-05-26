modules=%w{flexpart_helpers multi_rconc_w}

c_modules=modules.map {|m| m+".c"}
o_modules=modules.map {|m| m+".o"}
executable="flexpart2bin"

puts "Compiling programs"
`gcc -c #{c_modules.join(" ")}`
puts "Linking programs"
`gcc  #{o_modules.join(" ")} -o #{executable} `

puts `#{executable}`
