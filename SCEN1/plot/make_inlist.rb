# encoding: utf-8
dfiles=[]
File.open("input_list","w") do |fil|
  dfiles=Dir.glob("../output/grid_conc*").sort
  fil.puts dfiles
end

# Mining the time information
#../output/grid_conc_20110330080000
ftim=File.open("time_list","w") 
dfiles.each do |fil|
  dstring= fil.match(/_2011\d+/).to_s
  dstr_arr=dstring.delete("_").scan(/./)
  yr=dstr_arr[0..3].join.to_i
  mo=dstr_arr[4..5].join.to_i
  dy=dstr_arr[6..7].join.to_i
  hr=dstr_arr[8..9].join.to_i
  utime=Time.utc(yr,mo,dy,hr)
  time=utime.localtime
  #times<<[yr,mo,dy,hr]
  ftim.puts "#{time.year},#{time.mon}, #{time.day}, #{time.hour}"
end


