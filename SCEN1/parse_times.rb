require 'erb'
class FlxInit
  attr_accessor :sdate, :edate
  attr_accessor :stime, :etime
  attr_accessor :ipin, :ipout
  attr_accessor :wod

  def wrf_out_fil(dom,date)
    fse, fmi, fhr, fdy, fmo, fyr = date.to_a[0..5]
    fhr = "%02i" % fhr
    fdy = "%02i" % fdy
    fmo = "%02i" % fmo
    fmi = "%02i" % fmi
    fse = "%02i" % fse
    File.join(wod,"wrfout_#{dom}_#{fyr}-#{fmo}-#{fdy}_#{fhr}:#{fmi}:#{fse}")
  end

  def parse_wrf_times(fnam)

    string= `ncdump -v Times #{fnam}`
    tarr=[]
    string.split("\n").each do |line|
      tarr << line.chomp if line.match /^\s*\"\d\d\d\d-\d\d-\d\d/
    end
    tarr
  end

  def spc1
    "\s\s\s\s\s\s"
  end

  def spc2
    "\s\s\s\s"
  end

  def write_wrf_times(fnam,fout)
    tarr=parse_wrf_times(fnam)
    fnam=File.basename(fnam)
    tarr.each do |lin|

      lin=lin.delete("\"").delete(",").delete(";").strip
      date, time = lin.split("_")
      date.delete!("-")
      time.delete!(":")
      fout.puts "#{date} #{time}#{spc1}'#{fnam}'#{spc2}' '"
    end
  end

  def write_wrf_paths(pathfil,dom,wrf_out_dir,first)
    if first
      pathfil.puts "./options/"
      pathfil.puts "./output/"
    end
    pathfil.puts wrf_out_dir+"/"
    pathfil.puts "./AVAILABLE.#{dom}"
  end

  def write_flx_specs(template_file,optdir="./")
    template = ERB.new(IO.read(template_file))
    optdir=File.expand_path(optdir)
    optfil=File.join(optdir, File.basename(template_file,".template"))
    fopt = File.open(optfil,"w")
    fopt.puts template.result(self.get_binding)
    fopt.flush
    fopt.fsync
    fopt.close
  end

   def get_binding
    binding
  end
end
