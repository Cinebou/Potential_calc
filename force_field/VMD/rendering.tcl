#!/usr/local/bin/vmd


proc rendering_high_res {file_name} {
	
	render Tachyon $file_name\.dat "/home/auti/VMD/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" -aasamples 10 -auto_skylight 0.7 -fullshade %s -format BMP -res 6400 3600 -o $file_name\.bmp
}

proc rendering_low_res {file_name} {

	render snapshot $file_name exlporer %s

}