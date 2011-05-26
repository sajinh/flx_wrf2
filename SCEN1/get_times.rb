require 'erb'

class WRF_Times
  attr_accessor :lyr, :lmo, :ldy, :lhr
  attr_reader :yr, :mn, :dy, :hr
  attr_accessor :wrf_root
  attr_accessor :dom
  attr_reader  :fnam
  attr_accessor  :edate, :etime
  attr_accessor  :sdate, :stime
  attr_accessor  :ipin, :ipout
 
  def initialize(in_tmp=false)
    @in_tmp=in_tmp
    @ipin ||=0
    @ipout ||=1
  end

  def in_tmp?
    @in_tmp
  end

  def times_from_multiple_files
    abort "not yet implemented"
  end

  def get_wrf_times
    local_time=Time.local(lyr,lmo,ldy,lhr)
    utc_time=local_time.utc
    @dy= utc_time.day
    @mn= utc_time.mon
    @yr= utc_time.year
    @hr= utc_time.hour
    day="%02i" % dy
    mon="%02i" % mn
    hour="%02i" % hr
    if SPLIT_FILES
      return(times_from_multiple_files)
    end

    @fnam="wrfout_#{dom}_#{yr}-#{mon}-#{day}_#{hour}\:00\:00"
    string= `ncdump -v Times #{wrf_root}/#{@fnam}`
    tarr=[]
    string.split("\n").each do |line|
      tarr << line.chomp if line.match /^\s*\"\d\d\d\d-\d\d-\d\d/
    end
    tarr
  end

  def sub_dir
    if in_tmp?
      "run"
    else
      "out/#{yr}/#{mn}/#{dy}"
    end
  end

  def spc1
    "\s\s\s\s\s\s"
  end

  def spc2
    "\s\s\s\s"
  end

  def write_times
    tarr=get_wrf_times
    outfil=File.open("AVAILABLE.#{dom}","w")
    4.times { outfil.puts }
    istart=true
    tarr.each do |lin|
      lin=lin.delete("\"").delete(",").delete(";").strip
      date, time = lin.split("_")
      date.delete!("-")
      time.delete!(":")
      (@sdate, @stime = [date, time]) if istart
      istart=false
      outfil.puts "#{date} #{time}#{spc1}'#{fnam}'#{spc2}' '"
      @edate, @etime = date, time
    end
    outfil.close
  end

  def write_paths(path_fname,doms)
    pathfil=File.open(path_fname,"w")
    pathfil.puts "./options/"
    pathfil.puts "./output/"
    doms.each do |dom|
      pathfil.puts "#{wrf_root}/#{sub_dir}/"
      pathfil.puts "./AVAILABLE.#{dom}"
    end
    pathfil.puts "============================================"
    pathfil.puts
    pathfil.close
  end

  def write_specs(template_file,outdir="./")
    template = ERB.new(IO.read(template_file))
    outdir=File.expand_path(outdir)
    outfil=File.join(outdir, File.basename(template_file,".template"))
    fout = File.open(outfil,"w")
    fout.puts template.result(self.get_binding)
    fout.flush
    fout.fsync
    fout.close
  end

   def get_binding
    binding
  end

end
