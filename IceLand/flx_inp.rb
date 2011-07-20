require 'parse_times'
require 'fileutils'

# this script calls methods described in parse_times.rb
#  - it performs the following
#    + parses times from the wrfout_* files with the help of ncdump
#    + creates pathnames and AVAILABLE.* files
#    + creates the files options/COMMAND, options/RELEASES
#    + options/OUTGRID should be modified manually
#    + an NCL script print_latlon.ncl can help in determining
#    + the output grid dimensions
#  - see FLXPART_WRF documents for further instructions

# WRF initialization time (in local time!!! groan!!)
wrf_init=[2011,3,12]   # needed to located WRF output directory
                       # - this variable only used in the calling
                       # - script and helps construct name of
                       # - WRF output directory

# Output date for First file
first_wrfout_time = Time.utc(2011,5,22,0,0,0)
                  # this is also used within this script only
                  # - I write out WRF output every so many hours
                  # - this var is used to select the first file
                  # - i want to use for running FLXPART

# Starting date and time for FLXPART simulation (in UTC, please... so sorry)
flx_start = Time.utc(2011,5,22,6,0,0)
flx_end   = Time.utc(2011,5,22,12,0)
          # These specify start and end of release of particles
          # - currently it is set to be the same as
          # - start and end of FLXPART run
          # - anyway the start of release should coincide with
          # - start of FLXPART integration
          # - end of release can be different to stipulate a
          # - a puff kind of release 

# WRF domains for use with FLXPART
doms=%w{d01}

# Where is WRF output located?
#wrf_out_dir=File.join(ENV["HOME"],"Wrf-fks_coarse/out",wrf_init.join("/"))
wrf_out_dir=File.join(ENV["HOME"],"IceLand")

# Location of template files
cmd_template="./options/COMMAND.template"
rls_template="./options/RELEASES.template"

# Relevant directories
optdir="./options"
outdir="./output"


# --- perhaps there is no need to change lines below here ----- #

flx = FlxInit.new   # the following settings help fill up template files
flx.sdate = flx_start.strftime("%Y%m%d")
flx.stime = flx_start.strftime("%H%M%S")
flx.edate = flx_end.strftime("%Y%m%d")
flx.etime = flx_end.strftime("%k%M%S")
flx.ipin = 0
flx.ipout = 1


flx.wod = wrf_out_dir
first=true
pathfil=File.open("pathnames","w")
doms.each do |dom|

  fout = File.open("AVAILABLE.#{dom}","w")
  4.times { fout.puts }

  time_at_output = first_wrfout_time
  #fname=flx.wrf_out_fil(dom,time_at_output)
  while File.exist?(fname=flx.wrf_out_fil(dom,time_at_output))
    break if time_at_output > flx_end
    times = flx.write_wrf_times(fname, fout)
    time_at_output+= 86400 # advance 1 day
    print flx_end,"\t", time_at_output, "\n"
  end
  flx.write_wrf_paths(pathfil,dom,wrf_out_dir,first)
  first=false
  fout.close
end
pathfil.puts "============================================"
pathfil.puts
pathfil.close

# Write option/[command,releases] files
flx.write_flx_specs(cmd_template,optdir)
flx.write_flx_specs(rls_template,optdir)
