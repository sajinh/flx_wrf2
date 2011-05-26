def skip(fin,nlines)
  nlines.times.each {|l| fin.gets}
end

header="../output/header"
fhead=File.new(header,"r")
sdate, stime = fhead.gets.split
skip(fhead,2)

hdims=fhead.gets.split
xgrd, ygrd = hdims[2..3]
skip(fhead,1) if hdims.size <= 5
zgrd = fhead.gets.split[0]

fhead.close

File.open("dims.check","w+") {|fil| fil.puts [xgrd, ygrd, zgrd].join(",")}
