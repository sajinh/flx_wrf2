ROOT_DIR ='/home/saji'
FLX_ROOT=File.join(ROOT_DIR,"IceLand")
Dir.chdir(FLX_ROOT) do
  puts "Cleaning up output directory"
  system("rm -rf ./output/*")
  puts "Running FlexPart"
  system("./flexpart_wrf")

   Dir.chdir("./plot") do
     puts "Making input list"
     system("ruby make_inlist.rb")
     puts "Creating DDep data"
     system("ruby ddep.rb")
     puts "Creating WDep data"
     system("ruby wdep.rb")
     puts "Creating Conc data"
     system("ruby conc.rb")
   end
end

